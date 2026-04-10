/******************************************************************************
 *
 * Project:  GDAL
 * Purpose:  Zarr driver
 * Author:   Even Rouault <even dot rouault at spatialys.com>
 *
 ******************************************************************************
 * Copyright (c) 2025, Even Rouault <even dot rouault at spatialys.com>
 *
 * SPDX-License-Identifier: MIT
 ****************************************************************************/

#include "zarr.h"

#include "cpl_vsi.h"
#include "gdal_priv.h"

#include <cinttypes>

/************************************************************************/
/*                        ZarrTIFFDecompressor()                        */
/************************************************************************/

static bool ZarrTIFFDecompressor(const void *input_data, size_t input_size,
                                 void **output_data, size_t *output_size,
                                 CSLConstList /* options */,
                                 void * /* compressor_user_data */)
{
    if (output_data != nullptr && *output_data != nullptr &&
        output_size != nullptr && *output_size != 0)
    {
        const std::string osTmpFilename =
            VSIMemGenerateHiddenFilename("tmp.tif");
        VSIFCloseL(VSIFileFromMemBuffer(
            osTmpFilename.c_str(),
            const_cast<GByte *>(static_cast<const GByte *>(input_data)),
            input_size, /* bTakeOwnership = */ false));
        const char *const apszDrivers[] = {"GTIFF", "LIBERTIFF", nullptr};
        auto poDS = std::unique_ptr<GDALDataset>(GDALDataset::Open(
            osTmpFilename.c_str(), GDAL_OF_RASTER | GDAL_OF_VERBOSE_ERROR,
            apszDrivers, nullptr, nullptr));
        if (!poDS)
        {
            *output_size = 0;
            return false;
        }

        const int nBands = poDS->GetRasterCount();
        if (nBands != 1)
        {
            // This might be supported but I'm not sure which interleaving
            // should be returned !
            CPLError(CE_Failure, CPLE_NotSupported,
                     "ZarrTIFFDecompressor(): more than 1 band not supported");
            *output_size = 0;
            return false;
        }
        const int nXSize = poDS->GetRasterXSize();
        const int nYSize = poDS->GetRasterYSize();
        const GDALDataType eDT = poDS->GetRasterBand(1)->GetRasterDataType();
        const int nDTSize = GDALGetDataTypeSizeBytes(eDT);
        if (static_cast<uint64_t>(nXSize) * nYSize * nDTSize != *output_size)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "ZarrTIFFDecompressor(): %" PRIu64
                     " bytes expected, but %" PRIu64 " would be returned",
                     static_cast<uint64_t>(*output_size),
                     static_cast<uint64_t>(nXSize) * nYSize * nDTSize);
            *output_size = 0;
            return false;
        }

        const bool bOK = poDS->GetRasterBand(1)->RasterIO(
                             GF_Read, 0, 0, nXSize, nYSize, *output_data,
                             nXSize, nYSize, eDT, 0, 0, nullptr) == CE_None;
#ifdef CPL_MSB
        if (bOK && nDTSize > 1)
        {
            // Very likely we are expected to return in LSB order
            GDALSwapWordsEx(*output_data, nDTSize,
                            static_cast<size_t>(nXSize) * nYSize, nDTSize);
        }
#endif
        return bOK;
    }

    CPLError(CE_Failure, CPLE_AppDefined, "Invalid use of API");
    return false;
}

/************************************************************************/
/*                      ZarrGetTIFFDecompressor()                       */
/************************************************************************/

const CPLCompressor *ZarrGetTIFFDecompressor()
{
    static const CPLCompressor gTIFFDecompressor = {
        /* nStructVersion = */ 1,
        /* pszId = */ "imagecodecs_tiff",
        CCT_COMPRESSOR,
        /* papszMetadata = */ nullptr,
        ZarrTIFFDecompressor,
        /* user_data = */ nullptr};

    return &gTIFFDecompressor;
}

/************************************************************************/
/*                        ZarrJPEGDecompressor()                        */
/************************************************************************/

static bool ZarrJPEGDecompressor(const void *input_data, size_t input_size,
                                 void **output_data, size_t *output_size,
                                 CSLConstList /* options */,
                                 void * /* compressor_user_data */)
{
    if (output_data != nullptr && *output_data != nullptr &&
        output_size != nullptr && *output_size != 0)
    {
        const std::string osTmpFilename =
            VSIMemGenerateHiddenFilename("tmp.jpg");
        VSIFCloseL(VSIFileFromMemBuffer(
            osTmpFilename.c_str(),
            const_cast<GByte *>(static_cast<const GByte *>(input_data)),
            input_size, /* bTakeOwnership = */ false));
        const char *const apszDrivers[] = {"JPEG", nullptr};
        auto poDS = std::unique_ptr<GDALDataset>(GDALDataset::Open(
            osTmpFilename.c_str(), GDAL_OF_RASTER | GDAL_OF_VERBOSE_ERROR,
            apszDrivers, nullptr, nullptr));
        if (!poDS)
        {
            *output_size = 0;
            return false;
        }

        const int nBands = poDS->GetRasterCount();
        if (nBands != 1 && nBands != 3)
        {
            CPLError(
                CE_Failure, CPLE_NotSupported,
                "ZarrJPEGDecompressor(): only 1 or 3 bands supported, got %d",
                nBands);
            *output_size = 0;
            return false;
        }
        const int nXSize = poDS->GetRasterXSize();
        const int nYSize = poDS->GetRasterYSize();
        const GDALDataType eDT = poDS->GetRasterBand(1)->GetRasterDataType();
        const int nDTSize = GDALGetDataTypeSizeBytes(eDT);
        const uint64_t nExpectedSize =
            static_cast<uint64_t>(nXSize) * nYSize * nBands * nDTSize;
        if (nExpectedSize != *output_size)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "ZarrJPEGDecompressor(): %" PRIu64
                     " bytes expected, but %" PRIu64 " would be returned",
                     static_cast<uint64_t>(*output_size), nExpectedSize);
            *output_size = 0;
            return false;
        }

        // Pixel-interleaved layout: RGBRGBRGB...
        // nPixelSpace = bytes between adjacent pixels in same row    = nBands * nDTSize
        // nLineSpace  = bytes between adjacent rows                  = nXSize * nBands * nDTSize
        // nBandSpace  = bytes between same pixel in different bands  = nDTSize
        const GSpacing nPixelSpace = static_cast<GSpacing>(nBands) * nDTSize;
        const GSpacing nLineSpace =
            static_cast<GSpacing>(nXSize) * nBands * nDTSize;
        const GSpacing nBandSpace = static_cast<GSpacing>(nDTSize);

        const bool bOK =
            poDS->RasterIO(
                GF_Read, 0, 0, nXSize, nYSize,  // source window
                *output_data,                   // dest buffer
                nXSize, nYSize, eDT,            // dest size and type
                nBands, nullptr,  // band count, band map (null = 1..nBands)
                nPixelSpace, nLineSpace, nBandSpace, nullptr) == CE_None;
#ifdef CPL_MSB
        if (bOK && nDTSize > 1)
        {
            // Very likely we are expected to return in LSB order
            GDALSwapWordsEx(*output_data, nDTSize,
                            static_cast<size_t>(nXSize) * nYSize, nDTSize);
        }
#endif
        return bOK;
    }

    CPLError(CE_Failure, CPLE_AppDefined, "Invalid use of API");
    return false;
}

/************************************************************************/
/*                      ZarrGetJPEGDecompressor()                       */
/************************************************************************/

const CPLCompressor *ZarrGetJPEGDecompressor()
{
    static const CPLCompressor gJPEGDecompressor = {
        /* nStructVersion = */ 1,
        /* pszId = */ "imagecodecs_jpeg",
        CCT_COMPRESSOR,
        /* papszMetadata = */ nullptr,
        ZarrJPEGDecompressor,
        /* user_data = */ nullptr};

    return &gJPEGDecompressor;
}
