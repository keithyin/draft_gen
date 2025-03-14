// File Description
/// \file BamRecordTag.h
/// \brief Defines the BamRecordTag enum.
//
// Author: Derek Barnett

#ifndef BAMRECORDTAG_H
#define BAMRECORDTAG_H

namespace PacBio {
namespace BAM {

enum class BamRecordTag
{
    ALT_LABEL_QV,
    ALT_LABEL_TAG,
    BARCODE_QUALITY,
    BARCODES,
    CONTEXT_FLAGS,
    DELETION_QV,
    DELETION_TAG,
    HOLE_NUMBER,
    INSERTION_QV,
    IPD,
    LABEL_QV,
    LONG_CIGAR,
    MERGE_QV,
    NUM_PASSES,
    PKMEAN,
    PKMEAN_2,
    PKMID,
    PKMID_2,
    PRE_PULSE_FRAMES,
    PULSE_CALL,
    PULSE_CALL_WIDTH,
    PULSE_EXCLUSION,
    PULSE_MERGE_QV,
    PULSE_WIDTH,
    QUERY_END,
    QUERY_START,
    READ_ACCURACY,
    READ_GROUP,
    SCRAP_REGION_TYPE,
    SCRAP_ZMW_TYPE,
    SNR,
    BASE_LEVEL_CRS,
    BASE_LEVEL_DW_MEANS,
    BASE_LEVEL_DW_VARS,

    START_FRAME,
    SUBSTITUTION_QV,
    SUBSTITUTION_TAG,

    //
    // not tags per se, but faking these here to simplify data fetching
    //
    QUAL,
    SEQ,
    CAPTURE_RATE
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BAMRECORDTAG_H
