/******************************************************************************
 * RunImpl structural sketch — gdal mdim get-refs, Stage 1
 *
 * NOT finished code. A stage-by-stage skeleton showing where the seams fall.
 * Every non-obvious choice is tagged with the evidence-log finding behind it.
 * Fill in against gdalalg_raster_footprint.cpp (layer building) and
 * gdalalg_mdim_info.cpp / gdalalg_mdim_convert.cpp (mdim open + array resolve).
 ******************************************************************************/

#include "gdalalg_mdim_get_refs.h"
#include "cpl_conv.h"
#include "gdal_priv.h"

#ifndef _
#define _(x) (x)
#endif

/************************************************************************/
/*                        GDALMdimGetRefsAlgorithm()                  */
/************************************************************************/

GDALMdimGetRefsAlgorithm::GDALMdimGetRefsAlgorithm()
    : GDALAlgorithm(NAME, DESCRIPTION, HELP_URL)
{
    AddProgressArg();
    AddOutputFormatArg(&m_outputFormat, /* bStreamAllowed = */ false,
                       /* bGDALGAllowed = */ false)
        .AddMetadataItem(GAAMDI_REQUIRED_CAPABILITIES,
                         {GDAL_DCAP_VECTOR, GDAL_DCAP_CREATE});
    AddOpenOptionsArg(&m_openOptions);
    AddInputFormatsArg(&m_inputFormats)
        .AddMetadataItem(GAAMDI_REQUIRED_CAPABILITIES,
                         {GDAL_DCAP_MULTIDIM_RASTER});
    AddInputDatasetArg(&m_inputDataset, GDAL_OF_MULTIDIM_RASTER)
        .AddAlias("dataset");
    AddOutputDatasetArg(&m_outputDataset, GDAL_OF_VECTOR)
        .SetDatasetInputFlags(GADV_NAME | GADV_OBJECT);
    AddArrayNameArg(&m_array, _("Name of the array, used to restrict the "
                                "output to the specified array."))
        .SetRequired();
}

bool GDALMdimGetRefsAlgorithm::RunImpl(GDALProgressFunc pfnProgress, void *)
{
    // ----------------------------------------------------------------------
    // STAGE A — resolve the input array
    // ----------------------------------------------------------------------
    // A1. Get the input GDALDataset from m_inputDataset (already opened by the
    //     framework as OF_MULTIDIM_RASTER — confirm footprint/convert rely on
    //     the framework open rather than re-opening). GetDatasetRef().

    auto poSrcDS = m_inputDataset.GetDatasetRef();
    CPLAssert(poSrcDS);
    auto poRootGroup = poSrcDS->GetRootGroup();
    CPLDebug("get-refs", "input: %s, root group: %s", poSrcDS->GetDescription(),
             poRootGroup ? "present" : "NULL");
    // A2. GetRootGroup(). Null root => driver lacks mdim support => fail with a
    //     clear CPLError, return false.
    // A3. Resolve m_array against the root group. m_array is REQUIRED at
    //     Stage 1 (single-array contract). Support both a bare name and a
    //     '/'-prefixed full path — the HDF5 probe hit a deeply nested array
    //     (/HDFEOS/SWATHS/MySwath/Data Fields/MyDataField), so full-path
    //     resolution is not optional. OpenMDArrayFromFullname for '/'-prefixed,
    //     OpenMDArray from root otherwise. Null => fail clearly.
    //
    //     Decision: factor A3 into a helper now (resolveArray). It is the one
    //     piece of mdim plumbing the eventual whole-dataset traversal (deferred
    //     past Stage 1) will loop over — and the eventual classic-API widening
    //     will want a sibling of. Keep it in get_refs_common, owned here, not
    //     pushed to gcore.

    // ----------------------------------------------------------------------
    // STAGE B — describe the array (these facts become LAYER METADATA, not
    //           per-feature columns — evidence log Q4 + companion notes)
    // ----------------------------------------------------------------------
    // B1. GetDimensions() -> count, per-dim name + size.
    // B2. GetBlockSize() -> per-dim block extent.
    // B3. GUARD: any block extent == 0 => array is not chunk-enumerable. This
    //     is a VALID declined state (mosaic-VRT synthesised coord arrays report
    //     it — evidence log Q5/bonus), NOT a divide-by-zero to crash on. At
    //     Stage 1 / single-array, declining cleanly == fail with a clear
    //     message ("array X has no natural block size, not chunk-enumerable").
    //     When whole-dataset traversal arrives this becomes skip-with-warning.
    // B4. dtype: GetDataType(). C++ accessor still to confirm in
    //     gdal_multidim.h — the GetNumericDataType() route, NOT GetName()
    //     (empty for numerics). Low risk, not on the enumeration path.
    // B5. Hold these; they are written onto the OGRLayer via SetMetadataItem
    //     in Stage D. The codec chain in particular is array-level — it does
    //     NOT go in every row.

    // ----------------------------------------------------------------------
    // STAGE C — compute the chunk grid (CEIL division)
    // ----------------------------------------------------------------------
    // C1. n_chunks[i] = (dim_size[i] + block[i] - 1) / block[i]
    //     CEIL, confirmed for ZARR + netCDF + HDF5, all size-corroborated
    //     (evidence log Q1 + HDF5 addendum). Floor would silently drop every
    //     trailing partial chunk — on the HDF5 fixture that was 140 of 392.
    // C2. Total feature count = product(n_chunks). Compute up front so the
    //     progress callback (pfnProgress) has a denominator.
    // C3. No in-loop bounds paranoia needed: an out-of-range coordinate makes
    //     GetRawBlockInfo error loudly (evidence log bonus, all 3 drivers). If
    //     the ceil arithmetic is right every generated coord is valid; if it's
    //     wrong GDAL fails loudly rather than corrupting. Trust that.

    // ----------------------------------------------------------------------
    // STAGE D — create the output layer + field schema
    // ----------------------------------------------------------------------
    // D1. Output dataset from m_outputDataset (GADV_OUTPUT semantics — see
    //     footprint for the create-vs-update handling and how m_outputFormat
    //     selects the driver; default "Parquet").
    // D2. CreateLayer. Stage 1 = attribute-only, NO geometry column, NO SRS
    //     (geometry is Stage 2/3). wkbNone geometry type.
    // D3. Field schema (RFC section 4, Stage 1 table):
    //       dim_0 .. dim_n : OFTInteger / OFTInteger64  (chunk coord per dim)
    //       present        : OFTInteger, boolean subtype
    //       path           : OFTString   (nullable)
    //       offset         : OFTInteger64 (nullable)  <-- MUST be 64-bit:
    //                        BRAN had a real offset > 5 GB (evidence log bonus)
    //       size           : OFTInteger64 (nullable)
    //       info           : OFTString   (nullable, joined codec text)
    //     Open question carried from the RFC: whether per-row `info` survives
    //     once it's also hoisted to layer metadata. Keep it for Stage 1;
    //     decide later. Don't pre-resolve the open question in code.
    // D4. SetMetadataItem the array-level facts from Stage B onto the layer:
    //     dim names/sizes, block shape, dtype, and the codec chain (papszInfo
    //     joined / hoisted). This is the authoritative copy of the codec info.

    // ----------------------------------------------------------------------
    // STAGE E — enumerate the chunk grid, one feature per chunk
    // ----------------------------------------------------------------------
    // E1. One reusable GDALMDArrayRawBlockInfo, declared OUTSIDE the loop.
    //     It owns heap memory (pszFilename/papszInfo/pabyInlineData) — call
    //     .clear() at the top of each iteration, or rely on the fact that
    //     GetRawBlockInfo overwrites cleanly (CONFIRM which — clear() is the
    //     safe assumption). Do NOT accumulate filled structs.
    // E2. Build the uint64_t coordinate vector for this chunk (size = ndim),
    //     pass .data() to GetRawBlockInfo(coords, info).
    // E3. bool return false => the array declined this block. At Stage 1
    //     single-array, surface it (clear CPLError, fail) rather than guessing.
    // E4. Classify the three states from the REAL struct fields
    //     (evidence log Commit 1 reconciliation):
    //       present (file-backed): info.pszFilename != nullptr
    //       inline               : info.pszFilename == nullptr &&
    //                              info.pabyInlineData != nullptr
    //       absent (sparse)      : info.pszFilename == nullptr &&
    //                              info.pabyInlineData == nullptr
    //     Classify inline by pabyInlineData != nullptr, NOT by nSize > 0 —
    //     the struct's copy ctor can leave pabyInlineData NULL with nSize
    //     non-zero on alloc failure (documented). Stage 1 doesn't read inline
    //     bytes, but writing the check this way means Stage 1b inherits it.
    //     NEVER use nOffset == 0 as an absence signal — 0 is a legal offset
    //     (Zarr is one-file-per-chunk, every chunk at offset 0 — evidence
    //     log Q2).
    // E5. Create feature, populate:
    //       dim_* : the chunk coordinate
    //       present: true for file-backed AND inline; false for sparse
    //       path/offset/size:
    //         file-backed -> info.pszFilename, info.nOffset, info.nSize
    //         inline      -> path/offset NULL, size = info.nSize (Stage 1
    //                        reports size but NOT the bytes — Stage 1b adds
    //                        an OFTBinary field for pabyInlineData)
    //         sparse      -> path/offset/size all NULL
    //       info  : info.papszInfo joined (CSL helper), or NULL
    //     CreateFeature on the layer.
    // E6. Increment progress, honour pfnProgress; check for user interrupt.

    // ----------------------------------------------------------------------
    // STAGE F — finalize
    // ----------------------------------------------------------------------
    // F1. Layer/dataset flushing: the framework's Finalize() path closes the
    //     output dataset (see the C++ "gdal CLI from C++" note — Finalize()
    //     is what properly closes output datasets). Confirm whether RunImpl
    //     should explicitly flush the layer or leave it to framework Finalize.
    // F2. return true.

    return true;
}
