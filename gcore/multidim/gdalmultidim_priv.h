/******************************************************************************
 * Name:     gdalmultidim_priv.h
 * Project:  GDAL Core
 * Purpose:  GDAL private header for multidimensional support
 * Author:   Even Rouault <even.rouault at spatialys.com>
 *
 ******************************************************************************
 * Copyright (c) 2023, Even Rouault <even.rouault at spatialys.com>
 *
 * SPDX-License-Identifier: MIT
 ****************************************************************************/

#ifndef GDALMULTIDIM_PRIV_INCLUDED
#define GDALMULTIDIM_PRIV_INCLUDED

#include "gdal_priv.h"

//! @cond Doxygen_Suppress

/************************************************************************/
/*                           LinearToCoords()                           */
/************************************************************************/
/**
 * Decode a linear chunk index into per-dimension chunk coordinates,
 * using row-major (last dimension varies fastest) ordering.
 *
 * Given a chunk grid of shape n_chunks = [N_0, N_1, ..., N_{k-1}] and a
 * linear index iLinear in [0, product(n_chunks)), this fills coords with
 * the [c_0, c_1, ..., c_{k-1}] such that
 *     iLinear = (((c_0 * N_1) + c_1) * N_2 + c_2) * ... + c_{k-1}.
 *
 * The method is templated to handle variations at the caller.
 *
 * @param iLinear  The flat chunk index. Must be < product(n_chunks);
 *                 callers should bound this with nTotalChunks. No
 *                 internal range check done.
 * @param n_chunks Per-dimension chunk count. Must be non-empty and have
 *                 no zero entries.
 * @param coords   Output buffer; checked to match size of n_chunks, and filled.
 */
template <class T>
void LinearToCoords(T iLinear, const std::vector<T> &n_chunks,
                    const std::vector<T> &coords)
{
    CPLAssert(coords.size() == n_chunks.size());
    uint64_t remaining = iLinear;
    for (size_t i = nDims; i > 0;)
    {
        --i;
        coords[i] = iLinear % n_chunks[i];
        iLinear /= n_chunks[i];
    }
}

/************************************************************************/
/*                           CoordsToLinear()                           */
/************************************************************************/
/**
 * Encode per-dimension chunk coordinates into a linear chunk index,
 * inverse of LinearToCoords().
 *
 * @param coords    Per-dimension chunk coordinates. Each coords[i] must
 *                  be < n_chunks[i]; not checked.
 * @param n_chunks  Per-dimension chunk count. Same size as coords.
 * @return          The linear index in [0, product(n_chunks)).
 */
template <class T>
uint64_t CoordsToLinear(const std::vector<T> &coords,
                        const std::vector<T> &n_chunks)
{
    CPLAssert(coords.size() == n_chunks.size());
    uint64_t iLinear = 0;
    for (uint64_t iDim = 0; iDim < n_chunks.size(); ++iDim)
    {
        iLinear = iLinear * n_chunks[iDim] + coords[iDim];
    }
    return iLinear;
}

// For C API

struct GDALExtendedDataTypeHS
{
    std::unique_ptr<GDALExtendedDataType> m_poImpl;

    explicit GDALExtendedDataTypeHS(GDALExtendedDataType *dt) : m_poImpl(dt)
    {
    }
};

struct GDALEDTComponentHS
{
    std::unique_ptr<GDALEDTComponent> m_poImpl;

    explicit GDALEDTComponentHS(const GDALEDTComponent &component)
        : m_poImpl(new GDALEDTComponent(component))
    {
    }
};

struct GDALGroupHS
{
    std::shared_ptr<GDALGroup> m_poImpl;

    explicit GDALGroupHS(const std::shared_ptr<GDALGroup> &poGroup)
        : m_poImpl(poGroup)
    {
    }
};

struct GDALMDArrayHS
{
    std::shared_ptr<GDALMDArray> m_poImpl;

    explicit GDALMDArrayHS(const std::shared_ptr<GDALMDArray> &poArray)
        : m_poImpl(poArray)
    {
    }
};

struct GDALAttributeHS
{
    std::shared_ptr<GDALAttribute> m_poImpl;

    explicit GDALAttributeHS(const std::shared_ptr<GDALAttribute> &poAttr)
        : m_poImpl(poAttr)
    {
    }
};

struct GDALDimensionHS
{
    std::shared_ptr<GDALDimension> m_poImpl;

    explicit GDALDimensionHS(const std::shared_ptr<GDALDimension> &poDim)
        : m_poImpl(poDim)
    {
    }
};

class GDALMDIAsAttribute final : public GDALAttribute
{
    std::vector<std::shared_ptr<GDALDimension>> m_dims{};
    const GDALExtendedDataType m_dt = GDALExtendedDataType::CreateString();
    std::string m_osValue;

  public:
    GDALMDIAsAttribute(const std::string &name, const std::string &value)
        : GDALAbstractMDArray(std::string(), name),
          GDALAttribute(std::string(), name), m_osValue(value)
    {
    }

    const std::vector<std::shared_ptr<GDALDimension>> &
    GetDimensions() const override;

    const GDALExtendedDataType &GetDataType() const override
    {
        return m_dt;
    }

    bool IRead(const GUInt64 *, const size_t *, const GInt64 *,
               const GPtrDiff_t *, const GDALExtendedDataType &bufferDataType,
               void *pDstBuffer) const override
    {
        const char *pszStr = m_osValue.c_str();
        GDALExtendedDataType::CopyValue(&pszStr, m_dt, pDstBuffer,
                                        bufferDataType);
        return true;
    }
};

//! @endcond

#endif  // GDALMULTIDIM_PRIV_INCLUDED
