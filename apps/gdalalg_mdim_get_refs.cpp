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
#include "ogrsf_frmts.h"
#include "get_refs_common.h"

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
                         {GDAL_DCAP_VECTOR, GDAL_DCAP_CREATE})
        .SetRequired();
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
    AddOverwriteArg(&m_overwrite);
}

// Local helper for one-line vector formatting (used in debug + later metadata).
auto FormatVec = [](const std::vector<size_t> &v) -> CPLString
{
    CPLString os;
    for (size_t i = 0; i < v.size(); ++i)
    {
        if (i > 0)
            os += ", ";
        os += CPLSPrintf("%zu", v[i]);
    }
    return os;
};

bool GDALMdimGetRefsAlgorithm::RunImpl(GDALProgressFunc pfnProgress,
                                       void *pProgressData)
{
    // ----------------------------------------------------------------------
    // STAGE A — resolve the input array
    // ----------------------------------------------------------------------
    // A1. Get the input GDALDataset from m_inputDataset (already opened by the
    //     framework as OF_MULTIDIM_RASTER — confirm footprint/convert rely on
    //     the framework open rather than re-opening). GetDatasetRef().

    auto poSrcDS = m_inputDataset.GetDatasetRef();
    CPLAssert(poSrcDS);

    // A2. GetRootGroup(). Null root => driver lacks mdim support => fail with a
    //     clear CPLError, return false.
    auto poRootGroup = poSrcDS->GetRootGroup();
    CPLDebug("MDIM-GET-REFS", "input: %s, root group: %s",
             poSrcDS->GetDescription(), poRootGroup ? "present" : "NULL");

    if (!poRootGroup)
    {
        ReportError(CE_Failure, CPLE_AppDefined,
                    "Dataset %s has no root group (not multidimensional?)",
                    poSrcDS->GetDescription());
        return false;
    }

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
    //     pushed to gcore
    auto poArray = poRootGroup->OpenMDArrayFromFullname(m_array);
    if (!poArray)
    {
        ReportError(CE_Failure, CPLE_AppDefined,
                    "Cannot find array %s in dataset %s. \n"
                    "Use 'gdal mdim info %s' to list available arrays.",
                    m_array.c_str(), poSrcDS->GetDescription(),
                    poSrcDS->GetDescription());
        return false;
    }

    // ----------------------------------------------------------------------
    // STAGE B — describe the array (these facts become LAYER METADATA, not
    //           per-feature columns — evidence log Q4 + companion notes)
    // ----------------------------------------------------------------------
    // B1. GetDimensions() -> count, per-dim name + size.
    const std::vector<std::shared_ptr<GDALDimension>> apoDims =
        poArray->GetDimensions();

    // B2. GetBlockSize() -> per-dim block extent.
    // B3. GUARD: any block extent == 0 => array is not chunk-enumerable. This
    //     is a VALID declined state (mosaic-VRT synthesised coord arrays report
    //     it — evidence log Q5/bonus), NOT a divide-by-zero to crash on. At
    //     Stage 1 / single-array, declining cleanly == fail with a clear
    //     message ("array X has no natural block size, not chunk-enumerable").
    //     When whole-dataset traversal arrives this becomes skip-with-warning.
    const auto anBlockSize = poArray->GetBlockSize();
    std::vector<uint64_t> anBlockSizeU64(anBlockSize.begin(),
                                         anBlockSize.end());
    CPLAssert(anBlockSize.size() == poArray->GetDimensionCount());
    for (size_t i = 0; i < anBlockSizeU64.size(); ++i)
    {
        if (anBlockSize[i] == 0)
        {
            ReportError(CE_Failure, CPLE_AppDefined,
                        "Array %s has no natural block size on dimension %zu; "
                        "not chunk-enumerable",
                        m_array.c_str(), i);
            return false;
        }
    }

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

    const auto &dt = poArray->GetDataType();
    GDALDataType nDataType = dt.GetNumericDataType();
    const char *dt_name = GDALGetDataTypeName(nDataType);

    // Build the dim-size vector (needed by ComputeChunkGrid and for debug)
    std::vector<uint64_t> anDimSize(apoDims.size());
    for (size_t i = 0; i < apoDims.size(); ++i)
        anDimSize[i] = apoDims[i]->GetSize();

    // Stage C: chunk grid via the helper
    std::vector<size_t> n_chunks;
    const size_t nTotalChunks =
        get_refs::ComputeChunkGrid(anDimSize, anBlockSizeU64, n_chunks);

    // Stage C debug: one consolidated line using FormatVec on each vector
    CPLDebug(
        "MDIM-GET-REFS",
        "array %s: dims=[%s], blocks=[%s], chunks=[%s], total=%zu, dtype=%s",
        m_array.c_str(), FormatVec(anDimSize).c_str(),
        FormatVec(anBlockSizeU64).c_str(), FormatVec(n_chunks).c_str(),
        nTotalChunks, dt_name);
    // ----------------------------------------------------------------------
    // STAGE D — create the output dataset + layer + field schema
    // ----------------------------------------------------------------------
    // D1. Output path is what the user typed; the framework parsed and
    //     validated it but did not create the dataset (GADV_NAME | GADV_OBJECT
    //     stores intent, creation is RunImpl's job).
    const std::string osOutputPath = m_outputDataset.GetName();

    // D2. Driver lookup. m_outputFormat is .SetRequired() at the constructor,
    //     so it's guaranteed non-empty here.
    GDALDriver *poDriver =
        GetGDALDriverManager()->GetDriverByName(m_outputFormat.c_str());
    if (!poDriver)
    {
        ReportError(CE_Failure, CPLE_AppDefined,
                    "Cannot find vector driver '%s' for output dataset. "
                    "Use 'gdal --formats' to list available drivers.",
                    m_outputFormat.c_str());
        return false;
    }

    // D3. Create the dataset. 0,0,0 = no raster bands (vector dataset).
    //     unique_ptr so destructor calls GDALClose() and flushes to disk
    //     on every return path, including errors below.
    auto poDstDS = std::unique_ptr<GDALDataset>(
        poDriver->Create(osOutputPath.c_str(), 0, 0, 0, GDT_Unknown, nullptr));
    if (!poDstDS)
    {
        // GDALDriver::Create already emits a CPLError on failure;
        // no need to add another.
        return false;
    }

    // D4. Layer name from the array's basename — last '/'-separated segment.
    //     For "/HDFEOS/SWATHS/MySwath/Data Fields/MyDataField" → "MyDataField".
    //     For a bare-name array → the name itself.
    std::string osLayerName = m_array;
    const auto nLastSlash = osLayerName.rfind('/');
    if (nLastSlash != std::string::npos)
        osLayerName = osLayerName.substr(nLastSlash + 1);

    // D5. CreateLayer — wkbNone, no SRS at Stage 1 (attribute-only).
    OGRLayer *poLayer =
        poDstDS->CreateLayer(osLayerName.c_str(), nullptr, wkbNone, nullptr);
    if (!poLayer)
    {
        ReportError(CE_Failure, CPLE_AppDefined,
                    "Cannot create layer '%s' in output dataset '%s'",
                    osLayerName.c_str(), osOutputPath.c_str());
        return false;
    }

    // D6. Per-dimension fields: dim_0, dim_1, ... as OFTInteger64
    //     (names live in layer metadata, not field names — keeps the schema
    //     identical-shape across arrays, sidesteps sanitization).
    for (size_t i = 0; i < apoDims.size(); ++i)
    {
        OGRFieldDefn oField(CPLSPrintf("dim_%zu", i), OFTInteger64);
        if (poLayer->CreateField(&oField) != OGRERR_NONE)
        {
            ReportError(CE_Failure, CPLE_AppDefined,
                        "Cannot create field 'dim_%zu'", i);
            return false;
        }
    }

    // D7. Generic fields per RFC Stage 1 schema.
    {
        OGRFieldDefn oField("present", OFTInteger);
        oField.SetSubType(OFSTBoolean);
        if (poLayer->CreateField(&oField) != OGRERR_NONE)
        {
            ReportError(CE_Failure, CPLE_AppDefined,
                        "Cannot create field 'present'");
            return false;
        }
    }
    {
        OGRFieldDefn oField("path", OFTString);
        oField.SetNullable(TRUE);
        if (poLayer->CreateField(&oField) != OGRERR_NONE)
        {
            ReportError(CE_Failure, CPLE_AppDefined,
                        "Cannot create field 'path'");
            return false;
        }
    }
    {
        OGRFieldDefn oField("offset", OFTInteger64);
        oField.SetNullable(TRUE);
        if (poLayer->CreateField(&oField) != OGRERR_NONE)
        {
            ReportError(CE_Failure, CPLE_AppDefined,
                        "Cannot create field 'offset'");
            return false;
        }
    }
    {
        OGRFieldDefn oField("size", OFTInteger64);
        oField.SetNullable(TRUE);
        if (poLayer->CreateField(&oField) != OGRERR_NONE)
        {
            ReportError(CE_Failure, CPLE_AppDefined,
                        "Cannot create field 'size'");
            return false;
        }
    }
    {
        OGRFieldDefn oField("info", OFTString);
        oField.SetNullable(TRUE);
        if (poLayer->CreateField(&oField) != OGRERR_NONE)
        {
            ReportError(CE_Failure, CPLE_AppDefined,
                        "Cannot create field 'info'");
            return false;
        }
    }

    // D8. Array-level metadata on the layer.
    poLayer->SetMetadataItem("ARRAY_NAME", m_array.c_str());
    poLayer->SetMetadataItem("DTYPE", dt_name);
    for (size_t i = 0; i < apoDims.size(); ++i)
    {
        poLayer->SetMetadataItem(CPLSPrintf("DIM_%zu_NAME", i),
                                 apoDims[i]->GetName().c_str());
        poLayer->SetMetadataItem(
            CPLSPrintf("DIM_%zu_SIZE", i),
            CPLSPrintf(CPL_FRMT_GUIB,
                       static_cast<GUIntBig>(apoDims[i]->GetSize())));
        poLayer->SetMetadataItem(
            CPLSPrintf("DIM_%zu_BLOCK", i),
            CPLSPrintf(CPL_FRMT_GUIB, static_cast<GUIntBig>(anBlockSize[i])));
        poLayer->SetMetadataItem(CPLSPrintf("DIM_%zu_CHUNKS", i),
                                 CPLSPrintf("%zu", n_chunks[i]));
    }

    CPLDebug("MDIM-GET-REFS",
             "created layer '%s' with %d fields, ready for %zu features",
             osLayerName.c_str(), poLayer->GetLayerDefn()->GetFieldCount(),
             nTotalChunks);

    // Stage E (next): walk the chunk grid, call GetRawBlockInfo per chunk,
    // populate one feature per chunk, CreateFeature on poLayer.

    std::vector<uint64_t> coords;  // reused, resized inside helper
    GDALMDArrayRawBlockInfo info;  // also reused, .clear() per iteration

    // Progress is reported roughly every 1% of total chunks. For small arrays
    // (HDFEOS = 392) this is once-per-4-chunks; for large (BRAN = 94860) it is
    // once-per-~950. Either way ~100 progress callbacks per run, regardless
    // of array size.
    const size_t nProgressInterval = std::max<size_t>(1, nTotalChunks / 100);
    bool bCodecHoisted = false;
    for (size_t iLinear = 0; iLinear < nTotalChunks; ++iLinear)
    {
        info.clear();
        get_refs::LinearToCoords(iLinear, n_chunks, coords);
        if (!poArray->GetRawBlockInfo(coords.data(), info))
        {
            ReportError(CE_Failure, CPLE_AppDefined,
                        "GetRawBlockInfo failed at linear index %zu", iLinear);
            return false;
        }
        if (iLinear < 3 || iLinear == nTotalChunks - 1)
        {
            CPLDebug("MDIM-GET-REFS",
                     "chunk %zu coords=[%s] path=%s offset=" CPL_FRMT_GUIB
                     " size=" CPL_FRMT_GUIB,
                     iLinear, FormatVec(coords).c_str(),
                     info.pszFilename ? info.pszFilename : "(null)",
                     static_cast<GUIntBig>(info.nOffset),
                     static_cast<GUIntBig>(info.nSize));
        }
        OGRFeature *poFeature =
            OGRFeature::CreateFeature(poLayer->GetLayerDefn());

        // Per-dim coordinates: dim_0 .. dim_{n-1}
        for (size_t i = 0; i < coords.size(); ++i)
            poFeature->SetField(static_cast<int>(i),
                                static_cast<GIntBig>(coords[i]));

        // Three-state classification — exactly the Commit 1 reconciliation
        const int iPresentField = static_cast<int>(coords.size());
        const int iPathField = iPresentField + 1;
        const int iOffsetField = iPresentField + 2;
        const int iSizeField = iPresentField + 3;
        const int iInfoField = iPresentField + 4;

        if (info.pszFilename != nullptr)
        {
            // present (file-backed)
            poFeature->SetField(iPresentField, 1);
            poFeature->SetField(iPathField, info.pszFilename);
            poFeature->SetField(iOffsetField,
                                static_cast<GIntBig>(info.nOffset));
            poFeature->SetField(iSizeField, static_cast<GIntBig>(info.nSize));
        }
        else if (info.pabyInlineData != nullptr)
        {
            // inline — Stage 1 reports size but not bytes
            // (note: classify by pabyInlineData, not nSize > 0 — Commit 1)
            poFeature->SetField(iPresentField, 1);
            // path and offset left null (default state)
            poFeature->SetField(iSizeField, static_cast<GIntBig>(info.nSize));
        }
        else
        {
            // absent (sparse)
            poFeature->SetField(iPresentField, 0);
            // path, offset, size all left null
        }

        // info (papszInfo joined) — applies to all three states when non-null
        if (info.papszInfo != nullptr)
        {
            CPLStringList aosInfo(info.papszInfo, /* bAssign = */ false);
            // join key=value pairs into one string for the per-row field
            CPLString osJoined;
            for (int i = 0; i < aosInfo.size(); ++i)
            {
                if (i > 0)
                    osJoined += "; ";
                osJoined += aosInfo[i];
            }
            poFeature->SetField(iInfoField, osJoined.c_str());
        }

        // Codec hoist to layer metadata on the first successful file-backed chunk.
        // Per the design: the codec chain is array-level, not per-row authoritative.
        // The per-row info field stays (open question in the RFC; keep for Stage 1).
        if (!bCodecHoisted && info.papszInfo != nullptr &&
            info.pszFilename != nullptr)
        {
            for (int i = 0; info.papszInfo[i] != nullptr; ++i)
            {
                char *pszKey = nullptr;
                const char *pszValue =
                    CPLParseNameValue(info.papszInfo[i], &pszKey);
                if (pszKey && pszValue)
                    poLayer->SetMetadataItem(
                        CPLString().Printf("CODEC_%s", pszKey).c_str(),
                        pszValue);
                CPLFree(pszKey);
            }
            bCodecHoisted = true;
        }

        if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
        {
            ReportError(CE_Failure, CPLE_AppDefined,
                        "Cannot write feature for chunk %zu", iLinear);
            OGRFeature::DestroyFeature(poFeature);
            return false;
        }
        OGRFeature::DestroyFeature(poFeature);

        // Throttled progress report. pfnProgress may be null if no callback was
        // provided. Return-false from pfnProgress means the user (or environment)
        // has requested cancellation; treat as a clean failure.
        if (pfnProgress && (iLinear % nProgressInterval == 0))
        {
            const double dfFraction = static_cast<double>(iLinear) /
                                      static_cast<double>(nTotalChunks);
            if (!pfnProgress(dfFraction, nullptr, pProgressData))
            {
                ReportError(CE_Failure, CPLE_UserInterrupt,
                            "User interrupted at chunk %zu of %zu", iLinear,
                            nTotalChunks);
                return false;
            }
        }
    }
    // Final progress tick — completes the bar at 100%.
    if (pfnProgress)
        pfnProgress(1.0, nullptr, pProgressData);

    return true;
}
