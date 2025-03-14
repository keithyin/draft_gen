// File Description
/// \file BamHeader.cpp
/// \brief Implements the BamHeader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BamHeader.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <set>
#include <sstream>

#include <htslib/hts.h>

#include "StringUtils.h"
#include "Version.h"

namespace PacBio {
namespace BAM {
namespace internal {

static const std::string BamHeaderPrefixHD{"@HD"};
static const std::string BamHeaderPrefixSQ{"@SQ"};
static const std::string BamHeaderPrefixRG{"@RG"};
static const std::string BamHeaderPrefixPG{"@PG"};
static const std::string BamHeaderPrefixCO{"@CO"};

static const std::string BamHeaderTokenVN{"VN"};
static const std::string BamHeaderTokenSO{"SO"};
static const std::string BamHeaderTokenpb{"BN"};

static inline bool CheckSortOrder(const std::string& lhs, const std::string& rhs)
{
    return lhs == rhs;
}

static inline bool CheckPbVersion(const std::string& lhs, const std::string& rhs)
{
    return (Version{lhs} >= Version::Minimum && Version{rhs} >= Version::Minimum);
}

static inline bool CheckSequences(const std::string& sortOrder,
                                  const std::vector<SequenceInfo>& lhs,
                                  const std::vector<SequenceInfo>& rhs)
{
    return ((sortOrder == "coordinate") ? lhs == rhs : true);
}

static void EnsureCanMerge(const BamHeader& lhs, const BamHeader& rhs)
{
    // check compatibility
    const auto sortOrderOk = CheckSortOrder(lhs.SortOrder(), rhs.SortOrder());
    const auto pbVersionOk = CheckPbVersion(lhs.PacBioBamVersion(), rhs.PacBioBamVersion());
    const auto sequencesOk = CheckSequences(lhs.SortOrder(), lhs.Sequences(), rhs.Sequences());
    if (sortOrderOk && pbVersionOk && sequencesOk) return;

    // if any checks failed, format error message & throw
    std::ostringstream e;
    e << "could not merge BAM headers:\n";

    if (!sortOrderOk) {
        e << "  mismatched sort orders (@HD:SO) : (" << lhs.SortOrder() << ", " << rhs.SortOrder()
          << ")\n";
    }

    if (!pbVersionOk) {
        e << "  incompatible PacBio BAM versions (@HD:pb) : (" << lhs.PacBioBamVersion() << ", "
          << rhs.PacBioBamVersion() << ")\n";
    }

    if (!sequencesOk) e << "  mismatched sequence lists (@SQ entries)\n";

    throw std::runtime_error{e.str()};
}

}  // namespace internal

BamHeader::BamHeader(const std::string& samHeaderText)
    : d_{std::make_shared<internal::BamHeaderPrivate>()}
{
    std::istringstream s{samHeaderText};
    std::string line;
    std::string firstToken;
    while (std::getline(s, line)) {

        // skip if line is not long enough to contain true values
        if (line.length() < 5) continue;

        // determine token at beginning of line
        firstToken = line.substr(0, 3);

        if (firstToken == internal::BamHeaderPrefixHD) {

            // pop off '@HD\t', then split HD lines into tokens
            const auto tokens = internal::Split(line.substr(4), '\t');
            for (const auto& token : tokens) {
                const auto tokenTag = token.substr(0, 2);
                const auto tokenValue = token.substr(3);

                std::string tokenTagUpper(tokenTag);
                std::transform(tokenTag.begin(), tokenTag.end(), tokenTagUpper.begin(), toupper);

                // set header contents
                if (tokenTagUpper == internal::BamHeaderTokenVN)
                    Version(tokenValue);
                else if (tokenTagUpper == internal::BamHeaderTokenSO)
                    SortOrder(tokenValue);
                else if (tokenTagUpper == internal::BamHeaderTokenpb)
                    PacBioBamVersion(tokenValue);
            }

            // check for required tags
            if (Version().empty()) Version(std::string{hts_version()});
        }

        else if (firstToken == internal::BamHeaderPrefixSQ)
            AddSequence(SequenceInfo::FromSam(line));

        else if (firstToken == internal::BamHeaderPrefixRG)
            AddReadGroup(ReadGroupInfo::FromSam(line));

        else if (firstToken == internal::BamHeaderPrefixPG)
            AddProgram(ProgramInfo::FromSam(line));

        else if (firstToken == internal::BamHeaderPrefixCO)
            AddComment(line.substr(4));
    }
}

BamHeader& BamHeader::operator+=(const BamHeader& other)
{
    internal::EnsureCanMerge(*this, other);

    // merge read groups
    for (const auto& rg : other.ReadGroups()) {
        if (!HasReadGroup(rg.Id())) AddReadGroup(rg);
    }

    // merge programs
    for (const auto& pg : other.Programs()) {
        if (!HasProgram(pg.Id())) AddProgram(pg);
    }

    // merge comments
    for (const auto& comment : other.Comments())
        AddComment(comment);

    return *this;
}

BamHeader& BamHeader::AddSequence(SequenceInfo sequence)
{
    const std::string name = sequence.Name();
    d_->sequences_.push_back(std::move(sequence));
    d_->sequenceIdLookup_[name] = d_->sequences_.size() - 1;
    return *this;
}

BamHeader& BamHeader::ClearSequences()
{
    d_->sequenceIdLookup_.clear();
    d_->sequences_.clear();
    return *this;
}

BamHeader BamHeader::DeepCopy() const
{
    BamHeader result;
    result.d_->version_ = d_->version_;
    result.d_->pacbioBamVersion_ = d_->pacbioBamVersion_;
    result.d_->sortOrder_ = d_->sortOrder_;
    result.d_->headerLineCustom_ = d_->headerLineCustom_;
    result.d_->readGroups_ = d_->readGroups_;
    result.d_->programs_ = d_->programs_;
    result.d_->comments_ = d_->comments_;
    result.d_->sequences_ = d_->sequences_;
    result.d_->sequenceIdLookup_ = d_->sequenceIdLookup_;
    return result;
}

BamHeader& BamHeader::PacBioBamVersion(const std::string& version)
{
    d_->pacbioBamVersion_ = version;
    const internal::Version fileVersion{version};
    if (fileVersion < internal::Version::Minimum) {
        throw std::runtime_error{"invalid PacBio BAM version number (" + fileVersion.ToString() +
                                 ") is older than the minimum supported version (" +
                                 internal::Version::Minimum.ToString() + ")"};
    }
    return *this;
}

ProgramInfo BamHeader::Program(const std::string& id) const
{
    const auto iter = d_->programs_.find(id);
    if (iter == d_->programs_.cend()) throw std::runtime_error{"program ID not found"};
    return iter->second;
}

std::vector<std::string> BamHeader::ProgramIds() const
{
    std::vector<std::string> result;
    result.reserve(d_->programs_.size());
    for (const auto& pg : d_->programs_)
        result.push_back(pg.first);
    return result;
}

std::vector<ProgramInfo> BamHeader::Programs() const
{
    std::vector<ProgramInfo> result;
    result.reserve(d_->programs_.size());
    for (const auto& pg : d_->programs_)
        result.push_back(pg.second);
    return result;
}

BamHeader& BamHeader::Programs(std::vector<ProgramInfo> programs)
{
    d_->programs_.clear();
    for (const auto& pg : programs)
        d_->programs_[pg.Id()] = std::move(pg);
    return *this;
}

ReadGroupInfo BamHeader::ReadGroup(const std::string& id) const
{
    const auto iter = d_->readGroups_.find(id);
    if (iter == d_->readGroups_.cend()) throw std::runtime_error{"read group ID not found"};
    return iter->second;
}

std::vector<std::string> BamHeader::ReadGroupIds() const
{
    std::vector<std::string> result;
    result.reserve(d_->readGroups_.size());
    for (const auto& rg : d_->readGroups_)
        result.push_back(rg.first);
    return result;
}

std::vector<ReadGroupInfo> BamHeader::ReadGroups() const
{
    std::vector<ReadGroupInfo> result;
    result.reserve(d_->readGroups_.size());
    for (const auto& rg : d_->readGroups_)
        result.push_back(rg.second);
    return result;
}

BamHeader& BamHeader::ReadGroups(std::vector<ReadGroupInfo> readGroups)
{
    d_->readGroups_.clear();
    for (auto&& rg : readGroups)
        d_->readGroups_[rg.Id()] = std::move(rg);
    return *this;
}

SequenceInfo BamHeader::Sequence(const std::string& name) const
{
    // TODO: SequenceId(name) throws if not found, should we do so here as well?

    const auto iter = d_->sequenceIdLookup_.find(name);
    if (iter == d_->sequenceIdLookup_.cend()) return SequenceInfo();
    const auto index = iter->second;
    assert(index >= 0 && static_cast<size_t>(index) < d_->sequences_.size());
    return d_->sequences_.at(index);
}

int32_t BamHeader::SequenceId(const std::string& name) const
{
    const auto iter = d_->sequenceIdLookup_.find(name);
    if (iter == d_->sequenceIdLookup_.cend()) throw std::runtime_error{"sequence not found"};
    return iter->second;
}

std::vector<std::string> BamHeader::SequenceNames() const
{
    std::vector<std::string> result;
    result.reserve(d_->sequences_.size());
    for (const auto& seq : d_->sequences_)
        result.push_back(seq.Name());
    return result;
}

BamHeader& BamHeader::Sequences(std::vector<SequenceInfo> sequences)
{
    d_->sequences_.clear();
    for (auto&& seq : sequences)
        AddSequence(std::move(seq));
    return *this;
}

std::string BamHeader::ToSam() const
{
    // init stream
    std::ostringstream out;

    // @HD
    const auto outputVersion = (d_->version_.empty() ? std::string{hts_version()} : d_->version_);
    const auto outputSortOrder = (d_->sortOrder_.empty() ? std::string{"unknown"} : d_->sortOrder_);
    const auto outputPbBamVersion =
        (d_->pacbioBamVersion_.empty() ? internal::Version::Current.ToString()
                                       : d_->pacbioBamVersion_);

    out << internal::BamHeaderPrefixHD
        << internal::MakeSamTag(internal::BamHeaderTokenVN, outputVersion)
        << internal::MakeSamTag(internal::BamHeaderTokenSO, outputSortOrder)
        << internal::MakeSamTag(internal::BamHeaderTokenpb, outputPbBamVersion) << '\n';

    // @SQ
    for (const auto& seq : d_->sequences_)
        out << seq.ToSam() << '\n';

    // @RG
    for (const auto& rgIter : d_->readGroups_)
        out << rgIter.second.ToSam() << '\n';

    // @PG
    for (const auto& progIter : d_->programs_)
        out << progIter.second.ToSam() << '\n';

    // @CO
    for (const auto& comment : d_->comments_)
        out << internal::BamHeaderPrefixCO << '\t' << comment << '\n';

    // return result
    return out.str();
}

}  // namespace BAM
}  // namespace PacBio
