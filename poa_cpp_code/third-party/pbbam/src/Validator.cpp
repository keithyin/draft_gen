// File Description
/// \file Validator.cpp
/// \brief Implements the Validator class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/Validator.h"

#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/core/ignore_unused.hpp>

#include "ValidationErrors.h"
#include "Version.h"
#include "pbbam/BamFile.h"
#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/EntireFileQuery.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/ReadGroupInfo.h"

namespace PacBio {
namespace BAM {
namespace internal {

struct ilexcompare_wrapper
{
    bool operator()(const std::string& lhs, const std::string& rhs) const
    {
        return boost::ilexicographical_compare(lhs, rhs);
    }
};

// clang-format off
static const std::set<std::string, ilexcompare_wrapper> AcceptedSortOrders
{
    "unknown",
    "unsorted",
    "queryname",
    "coordinate"
};

static const std::set<std::string> AcceptedReadTypes
{
    "POLYMERASE",
    "HQREGION",
    "SUBREAD",
    "SMC",
    "SCRAP",
    "UNKNOWN"
};
// clang-format on

static void ValidateReadGroup(const ReadGroupInfo& rg, std::unique_ptr<ValidationErrors>& errors)
{
    const std::string id = rg.Id();

    // has required fields
    if (id.empty()) errors->AddReadGroupError(id, "missing ID");
    if (rg.MovieName().empty()) errors->AddReadGroupError(id, "missing movie name (PU tag)");
    // 3.0.2 adds required RG:PM - do not check for now, we'll add version-aware
    // validation down the road

    // description tag has required components
    if (rg.ReadType().empty()) errors->AddReadGroupError(id, "missing READTYPE in description");
    if (rg.BindingKit().empty()) errors->AddReadGroupError(id, "missing BINDINGKIT in description");
    if (rg.SequencingKit().empty())
        errors->AddReadGroupError(id, "missing SEQUENCINGKIT in description");
    if (rg.BasecallerVersion().empty())
        errors->AddReadGroupError(id, "missing BASECALLERVERSION in description");
    if (rg.FrameRateHz().empty())
        errors->AddReadGroupError(id, "missing FRAMERATEHZ in description");

    // stored ID matches expected ID (as calculated from movie & type)
    if (!id.empty()) {
        const auto expectedId = MakeReadGroupId(rg.MovieName(), rg.ReadType());
        if (expectedId != id) {
            const std::string msg{"stored ID: " + id + " does not match computed ID: " +
                                  expectedId};
            errors->AddReadGroupError(id, std::move(msg));
        }
    }

    // valid read type
    if (!rg.ReadType().empty()) {
        if (internal::AcceptedReadTypes.find(rg.ReadType()) == internal::AcceptedReadTypes.cend())
            errors->AddReadGroupError(id, "read type: " + rg.ReadType() + " is unknown");
    }

    // valid read chemistry (binding, sequencing, chemistry)
    if (!rg.BindingKit().empty() && !rg.SequencingKit().empty() &&
        !rg.BasecallerVersion().empty()) {
        try {
            auto chem = rg.SequencingChemistry();
            boost::ignore_unused(chem);
        } catch (std::exception& e) {
            errors->AddReadGroupError(id, e.what());
        }
    }

    // frame rate convertable to floating point
    if (!rg.FrameRateHz().empty()) {
        try {
            const float frameRate = std::stof(rg.FrameRateHz());
            boost::ignore_unused(frameRate);
        } catch (std::exception& e) {
            errors->AddReadGroupError(id, e.what());
        }
    }
}

static void ValidateHeader(const BamHeader& header, const std::string& filename,
                           std::unique_ptr<ValidationErrors>& errors)
{
    const std::string& fn = filename;

    // SAM/BAM version
    try {
        Version v(header.Version());
        boost::ignore_unused(v);
    } catch (std::exception& e) {
        errors->AddFileError(fn, std::string{"SAM version (@HD:VN) failed: "} + e.what());
    }

    // sort order
    const std::string sortOrder = header.SortOrder();
    if (AcceptedSortOrders.find(sortOrder) == AcceptedSortOrders.end())
        errors->AddFileError(fn, std::string{"unknown sort order: "} + sortOrder);

    // PacBio version
    try {
        const Version v{header.PacBioBamVersion()};
        const Version minimum{3, 0, 1};
        if (v < minimum) {

            std::string msg{"PacBioBAM version (@HD:pb) "};
            msg += v.ToString();
            msg += " is older than the minimum supported version (" + minimum.ToString() + ")";
            errors->AddFileError(fn, std::move(msg));
        }
    } catch (std::exception& e) {
        errors->AddFileError(
            fn, std::string{"PacBioBAM version (@HD:pb) failed to parse: "} + e.what());
    }

    // sequences?

    // read groups
    for (const ReadGroupInfo& rg : header.ReadGroups())
        ValidateReadGroup(rg, errors);
}

static void ValidateMetadata(const BamFile& file, std::unique_ptr<ValidationErrors>& errors)
{
    // filename
    const std::string fn{file.Filename()};
    if (fn == "-") {
        errors->AddFileError(fn,
                             "validation not is available for streamed BAM. Please "
                             "write to a file and run validation on it.");
        errors->ThrowErrors();  // quit early
    }
    if (boost::algorithm::ends_with(fn, ".bam") || boost::algorithm::ends_with(fn, ".bam.tmp")) {
        errors->AddFileError(fn, "non-standard file extension");
    }

    // EOF
    if (!file.HasEOF()) errors->AddFileError(fn, "missing end-of-file marker");

    // has PBI
    if (!file.PacBioIndexExists()) errors->AddFileError(fn, "missing PBI file");

    // header
    ValidateHeader(file.Header(), file.Filename(), errors);
}

void ValidateMappedRecord(const BamRecord& b, std::unique_ptr<ValidationErrors>& errors)
{
    const std::string name{b.FullName()};
    if (b.ReferenceStart() < 0) errors->AddRecordError(name, "mapped record position is invalid");
    if (b.ReferenceId() < 0) errors->AddRecordError(name, "mapped record reference ID is invalid");

    // what else??
}

void ValidateRecordCore(const BamRecord& b, std::unique_ptr<ValidationErrors>& errors)
{
    if (!IsCcsOrTranscript(b.Type())) {
        const auto qStart = b.QueryStart();
        const auto qEnd = b.QueryEnd();
        if (qStart >= qEnd) {
            errors->AddRecordError(b.FullName(), "queryStart (qs) should be < queryEnd (qe)");
        }
    }
}

void ValidateRecordReadGroup(const BamRecord& b, std::unique_ptr<ValidationErrors>& errors)
{
    try {
        auto rg = b.ReadGroup();
        boost::ignore_unused(rg);
    } catch (std::exception& e) {
        errors->AddRecordError(b.FullName(), e.what());
    }
}

void ValidateRecordRequiredTags(const BamRecord& b, std::unique_ptr<ValidationErrors>& errors)
{
    const auto name = b.FullName();
    const auto isCcsOrTranscript = IsCcsOrTranscript(b.Type());
    if (!isCcsOrTranscript) {
        // qe/qs
        const bool hasQueryStart = b.HasQueryStart();
        const bool hasQueryEnd = b.HasQueryEnd();
        if (hasQueryStart && hasQueryEnd) {
            const auto qStart = b.QueryStart();
            const auto qEnd = b.QueryEnd();
            if (qStart >= qEnd)
                errors->AddRecordError(name, "queryStart (qs) should be < queryEnd (qe)");
        } else {
            if (!hasQueryStart) errors->AddRecordError(name, "missing tag: qs (queryStart)");
            if (!hasQueryEnd) errors->AddRecordError(name, "missing tag: qe (queryEnd)");
        }
    }

    // zm
    if (!b.HasHoleNumber()) errors->AddRecordError(name, "missing tag: zm (ZMW hole number)");

    // np
    if (!b.HasNumPasses())
        errors->AddRecordError(name, "missing tag: np (num passes)");
    else {
        const auto numPasses = b.NumPasses();
        if (!isCcsOrTranscript && numPasses != 1)
            errors->AddRecordError(name, "np (numPasses) tag for non-CCS records should be 1");
    }

    // rq
    if (!b.HasReadAccuracy()) errors->AddRecordError(name, "missing tag: rq (read accuracy)");

    // sn
    if (!b.HasSignalToNoise())
        errors->AddRecordError(name, "missing tag: sn (signal-to-noise ratio)");
}

void ValidateRecordTagLengths(const BamRecord& b, std::unique_ptr<ValidationErrors>& errors)
{
    const auto name = b.FullName();
    const size_t expectedLength =
        (IsCcsOrTranscript(b.Type()) ? b.Sequence().size() : (b.QueryEnd() - b.QueryStart()));

    // check "per-base"-type data lengths are compatible
    if (b.Sequence().size() != expectedLength)
        errors->AddRecordError(name, "sequence length does not match expected length");

    if (b.HasDeletionQV()) {
        if (b.DeletionQV().size() != expectedLength)
            errors->AddTagLengthError(name, "DeletionQV", "dq", b.DeletionQV().size(),
                                      expectedLength);
    }
    if (b.HasDeletionTag()) {
        if (b.DeletionTag().size() != expectedLength)
            errors->AddTagLengthError(name, "DeletionTag", "dt", b.DeletionTag().size(),
                                      expectedLength);
    }
    if (b.HasInsertionQV()) {
        if (b.InsertionQV().size() != expectedLength)
            errors->AddTagLengthError(name, "InsertionQV", "iq", b.InsertionQV().size(),
                                      expectedLength);
    }
    if (b.HasMergeQV()) {
        if (b.MergeQV().size() != expectedLength)
            errors->AddTagLengthError(name, "MergeQV", "mq", b.MergeQV().size(), expectedLength);
    }
    if (b.HasSubstitutionQV()) {
        if (b.SubstitutionQV().size() != expectedLength)
            errors->AddTagLengthError(name, "SubstitutionQV", "sq", b.SubstitutionQV().size(),
                                      expectedLength);
    }
    if (b.HasSubstitutionTag()) {
        if (b.SubstitutionTag().size() != expectedLength)
            errors->AddTagLengthError(name, "SubstitutionTag", "st", b.SubstitutionTag().size(),
                                      expectedLength);
    }
    if (b.HasIPD()) {
        if (b.IPD().size() != expectedLength)
            errors->AddTagLengthError(name, "IPD", "ip", b.IPD().size(), expectedLength);
    }

    // NOTE: disabling "internal" tag checks for now, only checking "standard"
    //       PacBioBAM tags

    //    if (b.HasAltLabelQV()) {
    //        if (b.AltLabelQV().size() != expectedLength)
    //            errors->AddTagLengthError(name, "AltLabelQV", "pv", b.AltLabelQV().size(), expectedLength);
    //    }
    //    if (b.HasAltLabelTag()) {
    //        if (b.AltLabelTag().size() != expectedLength)
    //            errors->AddTagLengthError(name, "AltLabelTag", "pt", b.AltLabelTag().size(), expectedLength);
    //    }
    //    if (b.HasLabelQV()) {
    //        if (b.LabelQV().size() != expectedLength)
    //            errors->AddTagLengthError(name, "LabelQV", "pq", b.LabelQV().size(), expectedLength);
    //    }
    //    if (b.HasPkmean()) {
    //        if (b.Pkmean().size() != expectedLength)
    //            errors->AddTagLengthError(name, "Pkmean", "pa", b.Pkmean().size(), expectedLength);
    //    }
    //    if (b.HasPkmean2()) {
    //        if (b.Pkmean2().size() != expectedLength)
    //            errors->AddTagLengthError(name, "Pkmean2", "ps", b.Pkmean2().size(), expectedLength);
    //    }
    //    if (b.HasPkmid()) {
    //        if (b.Pkmid().size() != expectedLength)
    //            errors->AddTagLengthError(name, "Pkmid", "pm", b.Pkmid().size(), expectedLength);
    //    }
    //    if (b.HasPkmid2()) {
    //        if (b.Pkmid2().size() != expectedLength)
    //            errors->AddTagLengthError(name, "Pkmid2", "pi", b.Pkmid2().size(), expectedLength);
    //    }
    //    if (b.HasPrePulseFrames()) {
    //        if (b.PrePulseFrames().size() != expectedLength)
    //            errors->AddTagLengthError(name, "PrePulseFrames", "pd", b.PrePulseFrames().size(), expectedLength);
    //    }
    //    if (b.HasPulseCall()) {
    //        if (b.PulseCall().size() != expectedLength)
    //            errors->AddTagLengthError(name, "PulseCall", "pc", b.PulseCall().size(), expectedLength);
    //    }
    //    if (b.HasPulseCallWidth()) {
    //        if (b.PulseCallWidth().size() != expectedLength)
    //            errors->AddTagLengthError(name, "PulseCallWidth", "px", b.PulseCallWidth().size(), expectedLength);
    //    }
    //    if (b.HasPulseMergeQV()) {
    //        if (b.PulseMergeQV().size() != expectedLength)
    //            errors->AddTagLengthError(name, "PulseMergeQV", "pg", b.PulseMergeQV().size(), expectedLength);
    //    }
    //    if (b.HasPulseWidth()) {
    //        if (b.PulseWidth().size() != expectedLength)
    //            errors->AddTagLengthError(name, "PulseWidth", "pw", b.PulseWidth().size(), expectedLength);
    //    }
}

void ValidateUnmappedRecord(const BamRecord& b, std::unique_ptr<ValidationErrors>& errors)
{
    const std::string name{b.FullName()};
    if (b.ReferenceStart() != -1) errors->AddRecordError(name, "unmapped record has a position");
    if (b.ReferenceId() != -1) errors->AddRecordError(name, "unmapped record has a reference ID");
}

static void ValidateRecord(const BamRecord& b, std::unique_ptr<ValidationErrors>& errors)
{
    ValidateRecordCore(b, errors);
    ValidateRecordReadGroup(b, errors);
    ValidateRecordRequiredTags(b, errors);
    ValidateRecordTagLengths(b, errors);
    if (b.IsMapped())
        ValidateMappedRecord(b, errors);
    else
        ValidateUnmappedRecord(b, errors);
}

}  // namespace internal

using internal::ValidationErrors;

void Validator::Validate(const BamHeader& header, const size_t maxErrors)
{
    auto errors = std::make_unique<ValidationErrors>(maxErrors);
    internal::ValidateHeader(header, "unknown", errors);
    if (!errors->IsEmpty()) errors->ThrowErrors();
}

void Validator::Validate(const ReadGroupInfo& rg, const size_t maxErrors)
{
    auto errors = std::make_unique<ValidationErrors>(maxErrors);
    internal::ValidateReadGroup(rg, errors);
    if (!errors->IsEmpty()) errors->ThrowErrors();
}

void Validator::Validate(const BamRecord& b, const size_t maxErrors)
{
    auto errors = std::make_unique<ValidationErrors>(maxErrors);
    internal::ValidateRecord(b, errors);
    if (!errors->IsEmpty()) errors->ThrowErrors();
}

void Validator::ValidateEntireFile(const BamFile& file, const size_t maxErrors)
{
    auto errors = std::make_unique<ValidationErrors>(maxErrors);
    internal::ValidateMetadata(file, errors);

    EntireFileQuery query(file);
    for (const BamRecord& record : query)
        internal::ValidateRecord(record, errors);

    if (!errors->IsEmpty()) errors->ThrowErrors();
}

void Validator::ValidateFileMetadata(const BamFile& file, const size_t maxErrors)
{
    auto errors = std::make_unique<ValidationErrors>(maxErrors);
    internal::ValidateMetadata(file, errors);
    if (!errors->IsEmpty()) errors->ThrowErrors();
}

}  // namespace BAM
}  // namespace PacBio
