// File Description
/// \file BamRecord.cpp
/// \brief Implements the BamRecord class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BamRecord.h"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>

#include <htslib/sam.h>
#include <boost/numeric/conversion/cast.hpp>

#include "BamRecordTags.h"
#include "MemoryUtils.h"
#include "Pulse2BaseCache.h"
#include "SequenceUtils.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/ZmwTypeMap.h"
#include "pbbam/virtual/VirtualRegionTypeMap.h"

namespace PacBio {
namespace BAM {
namespace internal {

// record type names
static const std::string recordTypeName_ZMW{"ZMW"};
static const std::string recordTypeName_Polymerase{"POLYMERASE"};
static const std::string recordTypeName_HqRegion{"HQREGION"};
static const std::string recordTypeName_Subread{"SUBREAD"};
static const std::string recordTypeName_CCS{"SMC"};
static const std::string recordTypeName_Scrap{"SCRAP"};
static const std::string recordTypeName_Transcript{"TRANSCRIPT"};
static const std::string recordTypeName_Unknown{"UNKNOWN"};

static int32_t HoleNumberFromName(const std::string& fullName)
{
    const auto mainTokens = Split(fullName, '/');
    if (mainTokens.at(0) == "transcript") {
        if (mainTokens.size() != 2) throw std::runtime_error{"malformed transcript record name"};
        return std::stoi(mainTokens.at(1));
    } else {
        if (mainTokens.size() != 3) throw std::runtime_error("malformed record name");
        return std::stoi(mainTokens.at(1));
    }
}

static Position QueryEndFromName(const std::string& fullName)
{
    const auto mainTokens = Split(fullName, '/');
    if (mainTokens.size() != 3) throw std::runtime_error{"malformed record name"};

    const auto queryTokens = Split(mainTokens.at(2), '_');
    if (queryTokens.size() != 2) throw std::runtime_error{"malformed record name"};

    return stoi(queryTokens.at(1));
}

static Position QueryStartFromName(const std::string& fullName)
{
    const auto mainTokens = Split(fullName, '/');
    if (mainTokens.size() != 3) throw std::runtime_error{"malformed record name"};

    const auto queryTokens = Split(mainTokens.at(2), '_');
    if (queryTokens.size() != 2) throw std::runtime_error{"malformed record name"};

    return stoi(queryTokens.at(0));
}

static inline std::string Label(const BamRecordTag tag) { return BamRecordTags::LabelFor(tag); }

static BamRecordImpl* CreateOrEdit(const BamRecordTag tag, const Tag& value, BamRecordImpl* impl)
{
    if (impl->HasTag(tag))
        impl->EditTag(tag, value);
    else
        impl->AddTag(tag, value);
    return impl;
}

static std::pair<int32_t, int32_t> AlignedOffsets(const BamRecord& record, const int seqLength)
{
    int32_t startOffset = 0;
    int32_t endOffset = seqLength;

    const auto b = internal::BamRecordMemory::GetRawData(record);
    uint32_t* cigarData = bam_get_cigar(b.get());
    const size_t numCigarOps = b->core.n_cigar;
    if (numCigarOps > 0) {

        // start offset
        for (size_t i = 0; i < numCigarOps; ++i) {
            const auto type = static_cast<CigarOperationType>(bam_cigar_op(cigarData[i]));
            if (type == CigarOperationType::HARD_CLIP) {
                if (startOffset != 0 && startOffset != seqLength) {
                    startOffset = -1;
                    break;
                }
            } else if (type == CigarOperationType::SOFT_CLIP)
                startOffset += bam_cigar_oplen(cigarData[i]);
            else
                break;
        }

        // end offset
        for (int i = numCigarOps - 1; i >= 0; --i) {
            const auto type = static_cast<CigarOperationType>(bam_cigar_op(cigarData[i]));
            if (type == CigarOperationType::HARD_CLIP) {
                if (endOffset != 0 && endOffset != seqLength) {
                    endOffset = -1;
                    break;
                }
            } else if (type == CigarOperationType::SOFT_CLIP)
                endOffset -= bam_cigar_oplen(cigarData[i]);
            else
                break;
        }

        if (endOffset == 0) endOffset = seqLength;
    }
    return {startOffset, endOffset};
}

template <typename T>
T Clip(const T& input, const size_t pos, const size_t len)
{
    if (input.empty()) return {};
    return T{input.cbegin() + pos, input.cbegin() + pos + len};
}

template <typename T>
T ClipPulse(const T& input, internal::Pulse2BaseCache* p2bCache, const size_t pos, const size_t len)
{
    assert(p2bCache);
    if (input.empty()) return {};

    // find start
    size_t start = p2bCache->FindFirst();
    size_t basesSeen = 0;
    while (basesSeen < pos) {
        start = p2bCache->FindNext(start);
        ++basesSeen;
    }

    // find end
    size_t end = start;
    size_t seen = 1;
    while (seen < len) {
        end = p2bCache->FindNext(end);
        ++seen;
    }

    // return clipped
    return {input.cbegin() + start, input.cbegin() + end + 1};
}

template <class InputIt, class Size, class OutputIt>
OutputIt Move_N(InputIt first, Size count, OutputIt result)
{
    return std::move(first, first + count, result);
}

template <typename F, typename N>
static void ClipAndGapify(const BamRecordImpl& impl, const bool aligned, const bool exciseSoftClips,
                          F* seq, N paddingNullValue, N deletionNullValue)
{
    assert(seq);

    const bool clipOrGapRequested = aligned || exciseSoftClips;
    if (impl.IsMapped() && clipOrGapRequested) {
        // determine final container length
        auto incrementsOutputLength = [](const CigarOperationType type, const bool isAligned,
                                         const bool exciseSoftClipsFromAln) {
            if (type == CigarOperationType::HARD_CLIP ||
                type == CigarOperationType::REFERENCE_SKIP) {
                return false;
            } else if (type == CigarOperationType::SOFT_CLIP && exciseSoftClipsFromAln) {
                return false;
            } else if (!isAligned && (type == CigarOperationType::DELETION ||
                                      type == CigarOperationType::PADDING)) {
                return false;
            } else
                return true;
        };

        size_t outputLength = 0;
        const auto cigar = impl.CigarData();
        for (const CigarOperation& op : cigar) {
            if (incrementsOutputLength(op.Type(), aligned, exciseSoftClips))
                outputLength += op.Length();
        }

        // move original data to temp, prep output container size
        F originalSeq = std::move(*seq);
        seq->resize(outputLength);

        // apply CIGAR ops
        size_t srcIndex = 0;
        size_t dstIndex = 0;
        for (const CigarOperation& op : cigar) {
            const auto opType = op.Type();
            const auto opLength = op.Length();

            // nothing to do for hard-clipped & ref-skipped positions
            if (opType == CigarOperationType::HARD_CLIP ||
                opType == CigarOperationType::REFERENCE_SKIP) {
                continue;
            }

            // maybe skip soft-clipped positions
            else if (opType == CigarOperationType::SOFT_CLIP) {
                if (exciseSoftClips)
                    srcIndex += opLength;
                else {
                    Move_N(originalSeq.begin() + srcIndex, opLength, seq->begin() + dstIndex);
                    srcIndex += opLength;
                    dstIndex += opLength;
                }
            }

            // maybe add deletion/padding values
            else if (aligned && opType == CigarOperationType::DELETION) {
                for (size_t i = 0; i < opLength; ++i)
                    (*seq)[dstIndex++] = deletionNullValue;
            } else if (aligned && opType == CigarOperationType::PADDING) {
                for (size_t i = 0; i < opLength; ++i)
                    (*seq)[dstIndex++] = paddingNullValue;
            }

            // all other CIGAR ops
            else {
                Move_N(originalSeq.begin() + srcIndex, opLength, seq->begin() + dstIndex);
                srcIndex += opLength;
                dstIndex += opLength;
            }
        }
    }
}

static inline void ClipAndGapifyBases(const BamRecordImpl& impl, const bool aligned,
                                      const bool exciseSoftClips, std::string* seq)
{
    ClipAndGapify<std::string, char>(impl, aligned, exciseSoftClips, seq, '*', '-');
}

static inline void ClipAndGapifyFrames(const BamRecordImpl& impl, const bool aligned,
                                       const bool exciseSoftClips, Frames* frames)
{
    assert(frames);
    std::vector<uint16_t> data{std::move(frames->Data())};
    ClipAndGapify<std::vector<uint16_t>, uint16_t>(impl, aligned, exciseSoftClips, &data, 0, 0);
    frames->Data(data);
}

static inline void ClipAndGapifyPhotons(const BamRecordImpl& impl, const bool aligned,
                                        const bool exciseSoftClips, std::vector<float>* data)
{
    ClipAndGapify<std::vector<float>, float>(impl, aligned, exciseSoftClips, data, 0.0, 0.0);
}

static inline void ClipAndGapifyQualities(const BamRecordImpl& impl, const bool aligned,
                                          const bool exciseSoftClips, QualityValues* quals)
{
    ClipAndGapify<QualityValues, QualityValue>(impl, aligned, exciseSoftClips, quals,
                                               QualityValue(0), QualityValue(0));
}

static inline void ClipAndGapifyUInts(const BamRecordImpl& impl, const bool aligned,
                                      const bool exciseSoftClips, std::vector<uint32_t>* data)
{
    ClipAndGapify<std::vector<uint32_t>, uint32_t>(impl, aligned, exciseSoftClips, data, 0, 0);
}

static inline void ClipAndGapifyUInt8s(const BamRecordImpl& impl, const bool aligned,
                                       const bool exciseSoftClips, std::vector<uint8_t>* data)
{
    ClipAndGapify<std::vector<uint8_t>, uint8_t>(impl, aligned, exciseSoftClips, data, 0, 0);
}

static RecordType NameToType(const std::string& name)
{
    if (name == recordTypeName_Subread) return RecordType::SUBREAD;
    if (name == recordTypeName_ZMW || name == recordTypeName_Polymerase) return RecordType::ZMW;
    if (name == recordTypeName_HqRegion) return RecordType::HQREGION;
    if (name == recordTypeName_CCS) return RecordType::CCS;
    if (name == recordTypeName_Scrap) return RecordType::SCRAP;
    if (name == recordTypeName_Transcript) return RecordType::TRANSCRIPT;
    return RecordType::UNKNOWN;
}

static void OrientBasesAsRequested(std::string* bases, Orientation current, Orientation requested,
                                   bool isReverseStrand, bool isPulse)
{
    assert(bases);
    if (current != requested && isReverseStrand) {
        if (isPulse)
            internal::ReverseComplementCaseSens(*bases);
        else
            internal::ReverseComplement(*bases);
    }
}

template <typename Container>
inline void OrientTagDataAsRequested(Container* data, Orientation current, Orientation requested,
                                     bool isReverseStrand)
{
    assert(data);
    if (current != requested && isReverseStrand) std::reverse(data->begin(), data->end());
}

static inline bool ConsumesQuery(const CigarOperationType type)
{
    return (bam_cigar_type(static_cast<int>(type)) & 0x1) != 0;
}

static inline bool ConsumesReference(const CigarOperationType type)
{
    return (bam_cigar_type(static_cast<int>(type)) & 0x2) != 0;
}

}  // namespace internal

const float BamRecord::photonFactor = 10.0;

BamRecord::BamRecord()
    : alignedStart_{PacBio::BAM::UnmappedPosition}, alignedEnd_{PacBio::BAM::UnmappedPosition}
{
}

BamRecord::BamRecord(BamHeader header)
    : header_{std::move(header)}
    , alignedStart_{PacBio::BAM::UnmappedPosition}
    , alignedEnd_{PacBio::BAM::UnmappedPosition}
{
}

BamRecord::BamRecord(BamRecordImpl impl)
    : impl_{std::move(impl)}
    , alignedStart_{PacBio::BAM::UnmappedPosition}
    , alignedEnd_{PacBio::BAM::UnmappedPosition}
{
}

BamRecord::BamRecord(const BamRecord& other)
    : impl_{other.impl_}
    , header_{other.header_}
    , alignedStart_{other.alignedStart_}
    , alignedEnd_{other.alignedEnd_}
{
}

BamRecord::BamRecord(BamRecord&& other)
    : impl_{std::move(other.impl_)}
    , header_{std::move(other.header_)}
    , alignedStart_{std::move(other.alignedStart_)}
    , alignedEnd_{std::move(other.alignedEnd_)}
    , p2bCache_{std::move(other.p2bCache_)}
{
}

BamRecord& BamRecord::operator=(const BamRecord& other)
{
    if (this != &other) {
        impl_ = other.impl_;
        header_ = other.header_;
        alignedStart_ = other.alignedStart_;
        alignedEnd_ = other.alignedEnd_;
        p2bCache_.reset();  // just reset, for now at least
    }
    return *this;
}

BamRecord& BamRecord::operator=(BamRecord&& other)
{
    if (this != &other) {
        impl_ = std::move(other.impl_);
        header_ = std::move(other.header_);
        alignedStart_ = std::move(other.alignedStart_);
        alignedEnd_ = std::move(other.alignedEnd_);
        p2bCache_ = std::move(other.p2bCache_);
    }
    return *this;
}

BamRecord::~BamRecord() {}

Position BamRecord::AlignedEnd() const
{
    if (alignedEnd_ == PacBio::BAM::UnmappedPosition) CalculateAlignedPositions();
    return alignedEnd_;
}

Position BamRecord::AlignedStart() const
{
    if (alignedStart_ == PacBio::BAM::UnmappedPosition) CalculateAlignedPositions();
    return alignedStart_;
}

Strand BamRecord::AlignedStrand() const
{
    return impl_.IsReverseStrand() ? Strand::REVERSE : Strand::FORWARD;
}

QualityValues BamRecord::AltLabelQV(Orientation orientation, bool aligned, bool exciseSoftClips,
                                    PulseBehavior pulseBehavior) const
{
    return FetchQualities(BamRecordTag::ALT_LABEL_QV, orientation, aligned, exciseSoftClips,
                          pulseBehavior);
}

BamRecord& BamRecord::AltLabelQV(const QualityValues& altLabelQVs)
{
    internal::CreateOrEdit(BamRecordTag::ALT_LABEL_QV, altLabelQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::AltLabelTag(Orientation orientation, bool aligned, bool exciseSoftClips,
                                   PulseBehavior pulseBehavior) const
{
    return FetchBases(BamRecordTag::ALT_LABEL_TAG, orientation, aligned, exciseSoftClips,
                      pulseBehavior);
}

BamRecord& BamRecord::AltLabelTag(const std::string& tags)
{
    internal::CreateOrEdit(BamRecordTag::ALT_LABEL_TAG, tags, &impl_);
    return *this;
}

int16_t BamRecord::BarcodeForward() const { return Barcodes().first; }

int16_t BamRecord::BarcodeReverse() const { return Barcodes().second; }

uint8_t BamRecord::BarcodeQuality() const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::BARCODE_QUALITY);
    const auto bq = impl_.TagValue(tagName);
    if (bq.IsNull())
        return 0;  // ?? "missing" value for tags ?? should we consider boost::optional<T> for these kind of guys ??
    return bq.ToUInt8();
}

BamRecord& BamRecord::BarcodeQuality(const uint8_t quality)
{
    internal::CreateOrEdit(BamRecordTag::BARCODE_QUALITY, quality, &impl_);
    return *this;
}

std::pair<int16_t, int16_t> BamRecord::Barcodes() const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::BARCODES);
    const Tag bc = impl_.TagValue(tagName);
    if (bc.IsNull()) throw std::runtime_error{"barcode tag (bc) was requested but is missing"};

    // NOTE: barcodes are still stored, per the spec, as uint16, even though
    // we're now using them as int16_t in the API (bug 31511)
    //
    if (!bc.IsUInt16Array())
        throw std::runtime_error{
            "barcode tag (bc) is malformed: should be a uint16_t array of size==2."};
    const auto bcArray = bc.ToUInt16Array();
    if (bcArray.size() != 2)
        throw std::runtime_error{
            "barcode tag (bc) is malformed: should be a uint16_t array of size==2."};

    return {boost::numeric_cast<int16_t>(bcArray[0]), boost::numeric_cast<int16_t>(bcArray[1])};
}

BamRecord& BamRecord::Barcodes(const std::pair<int16_t, int16_t>& barcodeIds)
{
    const std::vector<uint16_t> data{boost::numeric_cast<uint16_t>(barcodeIds.first),
                                     boost::numeric_cast<uint16_t>(barcodeIds.second)};
    internal::CreateOrEdit(BamRecordTag::BARCODES, data, &impl_);
    return *this;
}

void BamRecord::CalculateAlignedPositions() const
{
    // reset
    ResetCachedPositions();

    // skip if unmapped, or has no queryStart/End
    if (!impl_.IsMapped()) return;

    // get the query start/end
    const auto seqLength = static_cast<int>(impl_.SequenceLength());
    const bool isCcsOrTranscript = IsCcsOrTranscript(Type());
    const Position qStart = isCcsOrTranscript ? 0 : QueryStart();
    const Position qEnd = isCcsOrTranscript ? seqLength : QueryEnd();

    if (qStart == PacBio::BAM::UnmappedPosition || qEnd == PacBio::BAM::UnmappedPosition) return;

    // determine clipped end ranges
    const auto alignedOffsets = internal::AlignedOffsets(*this, seqLength);
    const auto startOffset = alignedOffsets.first;
    const auto endOffset = alignedOffsets.second;
    if (endOffset == -1 || startOffset == -1) return;  // TODO: handle error more??

    // store aligned positions (polymerase read coordinates)
    if (impl_.IsReverseStrand()) {
        alignedStart_ = qStart + (seqLength - endOffset);
        alignedEnd_ = qEnd - startOffset;
    } else {
        alignedStart_ = qStart + startOffset;
        alignedEnd_ = qEnd - (seqLength - endOffset);
    }
}

void BamRecord::CalculatePulse2BaseCache() const
{
    // skip already calculated
    if (p2bCache_) return;

    // else try to calculate p2b cache.
    if (!HasPulseCall())
        throw std::runtime_error{"BamRecord cannot calculate pulse2base mapping without 'pc' tag."};
    const auto pulseCalls =
        FetchBases(BamRecordTag::PULSE_CALL, Orientation::NATIVE, false, false, PulseBehavior::ALL);
    p2bCache_ = std::make_unique<internal::Pulse2BaseCache>(pulseCalls);
}

Cigar BamRecord::CigarData(bool exciseAllClips) const
{
    auto isClippingOp = [](const CigarOperation& op) {
        const auto type = op.Type();
        return type == CigarOperationType::SOFT_CLIP || type == CigarOperationType::HARD_CLIP;
    };

    auto cigar = impl_.CigarData();
    if (exciseAllClips) {
        cigar.erase(std::remove_if(cigar.begin(), cigar.end(), isClippingOp), cigar.end());
    }
    return cigar;
}

BamRecord& BamRecord::Clip(const ClipType clipType, const Position start, const Position end)
{
    switch (clipType) {
        case ClipType::CLIP_NONE:
            return *this;
        case ClipType::CLIP_TO_QUERY:
            return ClipToQuery(start, end);
        case ClipType::CLIP_TO_REFERENCE:
            return ClipToReference(start, end);
        default:
            throw std::runtime_error{"unsupported clip type requested"};
    }
}

void BamRecord::ClipTags(const size_t clipFrom, const size_t clipLength)
{
    const auto ipdCodec = ReadGroup().IpdCodec();
    const auto pwCodec = ReadGroup().PulseWidthCodec();

    // update BAM tags
    TagCollection tags = impl_.Tags();
    if (HasDeletionQV())
        tags[internal::Label(BamRecordTag::DELETION_QV)] =
            internal::Clip(DeletionQV(Orientation::NATIVE), clipFrom, clipLength).Fastq();
    if (HasInsertionQV())
        tags[internal::Label(BamRecordTag::INSERTION_QV)] =
            internal::Clip(InsertionQV(Orientation::NATIVE), clipFrom, clipLength).Fastq();
    if (HasMergeQV())
        tags[internal::Label(BamRecordTag::MERGE_QV)] =
            internal::Clip(MergeQV(Orientation::NATIVE), clipFrom, clipLength).Fastq();
    if (HasSubstitutionQV())
        tags[internal::Label(BamRecordTag::SUBSTITUTION_QV)] =
            internal::Clip(SubstitutionQV(Orientation::NATIVE), clipFrom, clipLength).Fastq();
    if (HasIPD()) {
        if (ipdCodec == FrameCodec::RAW)
            tags[internal::Label(BamRecordTag::IPD)] =
                internal::Clip(IPD(Orientation::NATIVE).Data(), clipFrom, clipLength);
        else if (ipdCodec == FrameCodec::V1)
            tags[internal::Label(BamRecordTag::IPD)] =
                internal::Clip(IPD(Orientation::NATIVE).Encode(), clipFrom, clipLength);
    }
    if (HasPulseWidth()) {
        if (pwCodec == FrameCodec::RAW)
            tags[internal::Label(BamRecordTag::PULSE_WIDTH)] =
                internal::Clip(PulseWidth(Orientation::NATIVE).Data(), clipFrom, clipLength);
        else if (pwCodec == FrameCodec::V1)
            tags[internal::Label(BamRecordTag::PULSE_WIDTH)] =
                internal::Clip(PulseWidth(Orientation::NATIVE).Encode(), clipFrom, clipLength);
    }
    if (HasDeletionTag())
        tags[internal::Label(BamRecordTag::DELETION_TAG)] =
            internal::Clip(DeletionTag(Orientation::NATIVE), clipFrom, clipLength);
    if (HasSubstitutionTag())
        tags[internal::Label(BamRecordTag::SUBSTITUTION_TAG)] =
            internal::Clip(SubstitutionTag(Orientation::NATIVE), clipFrom, clipLength);

    // internal BAM tags
    if (HasPulseCall()) {

        // ensure p2bCache initialized
        CalculatePulse2BaseCache();
        internal::Pulse2BaseCache* p2bCache = p2bCache_.get();

        if (HasAltLabelQV())
            tags[internal::Label(BamRecordTag::ALT_LABEL_QV)] =
                internal::ClipPulse(AltLabelQV(Orientation::NATIVE), p2bCache, clipFrom, clipLength)
                    .Fastq();
        if (HasLabelQV())
            tags[internal::Label(BamRecordTag::LABEL_QV)] =
                internal::ClipPulse(LabelQV(Orientation::NATIVE), p2bCache, clipFrom, clipLength)
                    .Fastq();
        if (HasPulseMergeQV())
            tags[internal::Label(BamRecordTag::PULSE_MERGE_QV)] =
                internal::ClipPulse(PulseMergeQV(Orientation::NATIVE), p2bCache, clipFrom,
                                    clipLength)
                    .Fastq();
        if (HasAltLabelTag())
            tags[internal::Label(BamRecordTag::ALT_LABEL_TAG)] = internal::ClipPulse(
                AltLabelTag(Orientation::NATIVE), p2bCache, clipFrom, clipLength);
        if (HasPulseCall())
            tags[internal::Label(BamRecordTag::PULSE_CALL)] =
                internal::ClipPulse(PulseCall(Orientation::NATIVE), p2bCache, clipFrom, clipLength);
        if (HasPkmean())
            tags[internal::Label(BamRecordTag::PKMEAN)] = EncodePhotons(
                internal::ClipPulse(Pkmean(Orientation::NATIVE), p2bCache, clipFrom, clipLength));
        if (HasPkmid())
            tags[internal::Label(BamRecordTag::PKMID)] = EncodePhotons(
                internal::ClipPulse(Pkmid(Orientation::NATIVE), p2bCache, clipFrom, clipLength));
        if (HasPkmean2())
            tags[internal::Label(BamRecordTag::PKMEAN_2)] = EncodePhotons(
                internal::ClipPulse(Pkmean2(Orientation::NATIVE), p2bCache, clipFrom, clipLength));
        if (HasPkmid2())
            tags[internal::Label(BamRecordTag::PKMID_2)] = EncodePhotons(
                internal::ClipPulse(Pkmid2(Orientation::NATIVE), p2bCache, clipFrom, clipLength));
        if (HasPrePulseFrames())
            tags[internal::Label(BamRecordTag::PRE_PULSE_FRAMES)] = internal::ClipPulse(
                PrePulseFrames(Orientation::NATIVE).Data(), p2bCache, clipFrom, clipLength);
        if (HasPulseCallWidth())
            tags[internal::Label(BamRecordTag::PULSE_CALL_WIDTH)] = internal::ClipPulse(
                PulseCallWidth(Orientation::NATIVE).Data(), p2bCache, clipFrom, clipLength);
        if (HasStartFrame())
            tags[internal::Label(BamRecordTag::START_FRAME)] = internal::ClipPulse(
                StartFrame(Orientation::NATIVE), p2bCache, clipFrom, clipLength);
    }

    impl_.Tags(tags);
}

void BamRecord::ClipFields(const size_t clipFrom, const size_t clipLength)
{
    const bool isForwardStrand = (AlignedStrand() == Strand::FORWARD);

    // clip seq, quals
    std::string sequence{internal::Clip(Sequence(Orientation::NATIVE), clipFrom, clipLength)};
    QualityValues qualities{internal::Clip(Qualities(Orientation::NATIVE), clipFrom, clipLength)};
    if (!isForwardStrand) {
        internal::ReverseComplement(sequence);
        internal::Reverse(qualities);
    }
    impl_.SetSequenceAndQualities(sequence, qualities.Fastq());

    ClipTags(clipFrom, clipLength);
}

BamRecord& BamRecord::ClipToQuery(const Position start, const Position end)
{
    // cache original coords, skip out if clip not needed
    const auto seqLength = static_cast<int>(impl_.SequenceLength());
    const bool isCcsOrTranscript = IsCcsOrTranscript(Type());
    const Position origQStart = isCcsOrTranscript ? 0 : QueryStart();
    const Position origQEnd = isCcsOrTranscript ? seqLength : QueryEnd();
    if (start <= origQStart && end >= origQEnd) return *this;

    // determine new offsets into data
    const size_t startOffset = start - origQStart;
    const size_t endOffset = origQEnd - end;

    // maybe update CIGAR & aligned position
    if (IsMapped()) {

        // fetch a 'working copy' of CIGAR data
        Cigar cigar = impl_.CigarData();

        // clip leading CIGAR ops
        size_t referencePositionOffset = 0;
        size_t remaining = startOffset;
        while (remaining > 0 && !cigar.empty()) {
            CigarOperation& firstOp = cigar.front();
            const auto firstOpLength = firstOp.Length();
            const bool consumesQuery = internal::ConsumesQuery(firstOp.Type());
            const bool consumesRef = internal::ConsumesReference(firstOp.Type());

            // if (!consumesQuery)
            //    just pop (e.g. deletion) ?
            // else {
            //    check bounds, like clip to reference ?
            // }

            // CIGAR op ends at or before clip
            if (firstOpLength <= remaining) {
                cigar.erase(cigar.begin());
                if (consumesQuery) remaining -= firstOpLength;
                if (consumesRef) referencePositionOffset += firstOpLength;
            }

            // CIGAR op straddles clip
            else {
                firstOp.Length(firstOpLength - remaining);
                if (consumesRef) referencePositionOffset += remaining;
                remaining = 0;
            }
        }

        // clip trailing CIGAR ops
        remaining = endOffset;
        while (remaining > 0 && !cigar.empty()) {
            CigarOperation& lastOp = cigar.back();
            const auto lastOpLength = lastOp.Length();
            const bool consumesQuery = internal::ConsumesQuery(lastOp.Type());

            // CIGAR op ends at or after clip
            if (lastOpLength <= remaining) {
                cigar.pop_back();
                if (consumesQuery) remaining -= lastOpLength;
            }

            // CIGAR op straddles clip
            else {
                lastOp.Length(lastOpLength - remaining);
                remaining = 0;
            }
        }

        // update CIGAR & position
        impl_.CigarData(cigar);
        impl_.Position(impl_.Position() + referencePositionOffset);
    }

    // clip SEQ, QUAL, & tags
    const size_t clipFrom = startOffset;
    const size_t clipLength = (end - start);
    ClipFields(clipFrom, clipLength);

    // update query start/end
    // TODO: update name to reflect new QS/QE ???
    internal::CreateOrEdit(BamRecordTag::QUERY_START, start, &impl_);
    internal::CreateOrEdit(BamRecordTag::QUERY_END, end, &impl_);
    //    UpdateName();

    // reset any cached aligned start/end
    ResetCachedPositions();
    return *this;
}

BamRecord& BamRecord::ClipToReference(const Position start, const Position end)
{
    // skip if not mapped, clipping to reference doesn't make sense
    // or should we even consider throwing here?
    if (!IsMapped()) return *this;

    const bool isForwardStrand = (AlignedStrand() == Strand::FORWARD);
    return (isForwardStrand ? ClipToReferenceForward(start, end)
                            : ClipToReferenceReverse(start, end));
}

BamRecord& BamRecord::ClipToReferenceForward(const PacBio::BAM::Position start,
                                             const PacBio::BAM::Position end)
{
    assert(IsMapped());
    assert(AlignedStrand() == Strand::FORWARD);

    // cache original coords
    const int seqLength = static_cast<int>(impl_.SequenceLength());
    const bool isCcsOrTranscript = IsCcsOrTranscript(Type());
    const Position origQStart = isCcsOrTranscript ? 0 : QueryStart();
    const Position origQEnd = isCcsOrTranscript ? seqLength : QueryEnd();
    const Position origTStart = ReferenceStart();
    const Position origTEnd = ReferenceEnd();
    assert(AlignedStart() >= origQStart);
    assert(AlignedEnd() <= origQEnd);

    // skip if already within requested clip range
    if (start <= origTStart && end >= origTEnd) return *this;

    const Position newTStart = std::max(origTStart, start);
    const Position newTEnd = std::min(origTEnd, end);

    // fetch a 'working copy' of CIGAR data
    Cigar cigar = impl_.CigarData();

    // we're going to skip query sequence outside aligned region
    size_t queryPosRemovedFront = 0;
    size_t queryPosRemovedBack = 0;

    // ------------------------
    // clip leading CIGAR ops
    // ------------------------

    size_t remaining = newTStart - origTStart;
    while (remaining > 0 && !cigar.empty()) {
        CigarOperation& firstOp = cigar.front();
        const auto firstOpLength = firstOp.Length();
        const bool consumesQuery = internal::ConsumesQuery(firstOp.Type());
        const bool consumesRef = internal::ConsumesReference(firstOp.Type());

        if (!consumesRef) {

            // e.g. softclip - just pop it completely
            cigar.erase(cigar.begin());
            if (consumesQuery) queryPosRemovedFront += firstOpLength;

        } else {
            assert(consumesRef);

            // CIGAR ends at or before clip
            if (firstOpLength <= remaining) {
                cigar.erase(cigar.begin());
                if (consumesQuery) queryPosRemovedFront += firstOpLength;
                if (consumesRef) remaining -= firstOpLength;
            }

            // CIGAR straddles clip
            else {
                assert(firstOpLength > remaining);
                firstOp.Length(firstOpLength - remaining);
                if (consumesQuery) queryPosRemovedFront += remaining;
                remaining = 0;
            }
        }
    }

    // -------------------------
    // clip trailing CIGAR ops
    // -------------------------

    remaining = origTEnd - newTEnd;
    while (remaining > 0 && !cigar.empty()) {
        CigarOperation& lastOp = cigar.back();
        const auto lastOpLength = lastOp.Length();
        const bool consumesQuery = internal::ConsumesQuery(lastOp.Type());
        const bool consumesRef = internal::ConsumesReference(lastOp.Type());

        if (!consumesRef) {

            // e.g. softclip - just pop it completely
            cigar.pop_back();
            if (consumesQuery) queryPosRemovedBack += lastOpLength;

        } else {
            assert(consumesRef);

            // CIGAR ends at or after clip
            if (lastOpLength <= remaining) {
                cigar.pop_back();
                if (consumesQuery) queryPosRemovedBack += lastOpLength;
                if (consumesRef) remaining -= lastOpLength;
            }

            // CIGAR straddles clip
            else {
                assert(lastOpLength > remaining);
                lastOp.Length(lastOpLength - remaining);
                if (consumesQuery) queryPosRemovedBack += remaining;
                remaining = 0;
            }
        }
    }

    // update CIGAR and position
    impl_.CigarData(cigar);
    impl_.Position(newTStart);

    // clip SEQ, QUAL, tags
    const Position qStart = origQStart + queryPosRemovedFront;
    const Position qEnd = origQEnd - queryPosRemovedBack;
    const size_t clipFrom = queryPosRemovedFront;
    const size_t clipLength = qEnd - qStart;
    ClipFields(clipFrom, clipLength);

    // update query start/end
    internal::CreateOrEdit(BamRecordTag::QUERY_START, qStart, &impl_);
    internal::CreateOrEdit(BamRecordTag::QUERY_END, qEnd, &impl_);
    //    UpdateName();

    // reset any cached aligned start/end
    ResetCachedPositions();
    return *this;
}

BamRecord& BamRecord::ClipToReferenceReverse(const PacBio::BAM::Position start,
                                             const PacBio::BAM::Position end)
{
    assert(IsMapped());
    assert(AlignedStrand() == Strand::REVERSE);

    // cache original coords
    const int seqLength = static_cast<int>(impl_.SequenceLength());
    const bool isCcsOrTranscript = IsCcsOrTranscript(Type());
    const Position origQStart = isCcsOrTranscript ? 0 : QueryStart();
    const Position origQEnd = isCcsOrTranscript ? seqLength : QueryEnd();
    const Position origTStart = ReferenceStart();
    const Position origTEnd = ReferenceEnd();

    // skip if already within requested clip range
    if (start <= origTStart && end >= origTEnd) return *this;
    assert(AlignedStart() >= origQStart);
    assert(AlignedEnd() <= origQEnd);

    const Position newTStart = std::max(origTStart, start);
    const Position newTEnd = std::min(origTEnd, end);

    Cigar cigar = impl_.CigarData();

    size_t queryPosRemovedFront = 0;
    size_t queryPosRemovedBack = 0;

    // update CIGAR - clip front ops, then clip back ops
    size_t remaining = newTStart - origTStart;
    while (remaining > 0 && !cigar.empty()) {
        CigarOperation& firstOp = cigar.front();
        const auto firstOpType = firstOp.Type();
        const auto firstOpLength = firstOp.Length();
        const bool consumesQuery = internal::ConsumesQuery(firstOpType);
        const bool consumesRef = internal::ConsumesReference(firstOpType);

        if (!consumesRef) {

            // e.g. softclip - just pop it completely
            cigar.erase(cigar.begin());
            if (consumesQuery) queryPosRemovedBack += firstOpLength;

        } else {
            assert(consumesRef);

            // CIGAR ends at or before clip
            if (firstOpLength <= remaining) {
                cigar.erase(cigar.begin());
                if (consumesQuery) queryPosRemovedBack += firstOpLength;
                if (consumesRef) remaining -= firstOpLength;
            }

            // CIGAR straddles clip
            else {
                assert(firstOpLength > remaining);
                firstOp.Length(firstOpLength - remaining);
                if (consumesQuery) queryPosRemovedBack += remaining;
                remaining = 0;
            }
        }
    }

    remaining = origTEnd - newTEnd;
    while (remaining > 0 && !cigar.empty()) {
        CigarOperation& lastOp = cigar.back();
        const auto lastOpType = lastOp.Type();
        const auto lastOpLength = lastOp.Length();
        const bool consumesQuery = internal::ConsumesQuery(lastOpType);
        const bool consumesRef = internal::ConsumesReference(lastOpType);

        if (!consumesRef) {

            // e.g. softclip - just pop it completely
            cigar.pop_back();
            if (consumesQuery) queryPosRemovedFront += lastOpLength;

        } else {
            assert(consumesRef);

            // CIGAR ends at or before clip
            if (lastOpLength <= remaining) {
                cigar.pop_back();
                if (consumesQuery) queryPosRemovedFront += lastOpLength;
                if (consumesRef) remaining -= lastOpLength;
            }

            // CIGAR straddles clip
            else {
                assert(lastOpLength > remaining);
                lastOp.Length(lastOpLength - remaining);
                if (consumesQuery) queryPosRemovedFront += remaining;
                remaining = 0;
            }
        }
    }
    impl_.CigarData(cigar);

    // update aligned reference position
    impl_.Position(newTStart);

    // clip SEQ, QUAL, tags
    const Position qStart = origQStart + queryPosRemovedFront;
    const Position qEnd = origQEnd - queryPosRemovedBack;
    const size_t clipFrom = queryPosRemovedFront;
    const size_t clipLength = qEnd - qStart;
    ClipFields(clipFrom, clipLength);

    // update query start/end
    internal::CreateOrEdit(BamRecordTag::QUERY_START, qStart, &impl_);
    internal::CreateOrEdit(BamRecordTag::QUERY_END, qEnd, &impl_);
    //    UpdateName();

    // reset any cached aligned start/end
    ResetCachedPositions();
    return *this;
}

QualityValues BamRecord::DeletionQV(Orientation orientation, bool aligned,
                                    bool exciseSoftClips) const
{
    return FetchQualities(BamRecordTag::DELETION_QV, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::DeletionQV(const QualityValues& deletionQVs)
{
    internal::CreateOrEdit(BamRecordTag::DELETION_QV, deletionQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::DeletionTag(Orientation orientation, bool aligned,
                                   bool exciseSoftClips) const
{
    return FetchBases(BamRecordTag::DELETION_TAG, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::DeletionTag(const std::string& tags)
{
    internal::CreateOrEdit(BamRecordTag::DELETION_TAG, tags, &impl_);
    return *this;
}

std::vector<uint16_t> BamRecord::EncodePhotons(const std::vector<float>& data)
{
    std::vector<uint16_t> encoded;
    encoded.reserve(data.size());
    for (const auto& d : data)
        encoded.emplace_back(d * photonFactor);
    return encoded;
}

std::string BamRecord::FetchBasesRaw(const BamRecordTag tag) const
{
    const Tag seqTag = impl_.TagValue(tag);
    return seqTag.ToString();
}

std::string BamRecord::FetchBases(const BamRecordTag tag, const Orientation orientation,
                                  const bool aligned, const bool exciseSoftClips,
                                  const PulseBehavior pulseBehavior) const
{
    const bool isBamSeq = (tag == BamRecordTag::SEQ);
    const bool isPulse = internal::BamRecordTags::IsPulse(tag);

    // fetch raw
    std::string bases;
    Orientation current;
    if (isBamSeq) {  // SEQ stored in genomic orientation
        bases = impl_.Sequence();
        current = Orientation::GENOMIC;
    } else {  // all tags stored in native orientation
        bases = FetchBasesRaw(tag);
        current = Orientation::NATIVE;
    }

    // maybe strip 'squashed' pulse loci
    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        CalculatePulse2BaseCache();
        bases = p2bCache_->RemoveSquashedPulses(bases);
    }

    // if we need to touch CIGAR
    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY)
            throw std::runtime_error{
                "Cannot return data at all pulses when gapping and/or soft-clipping are requested. "
                "Use PulseBehavior::BASECALLS_ONLY instead."};

        // force into genomic orientation
        internal::OrientBasesAsRequested(&bases, current, Orientation::GENOMIC,
                                         impl_.IsReverseStrand(), isPulse);
        current = Orientation::GENOMIC;

        // clip & gapify as requested
        internal::ClipAndGapifyBases(impl_, aligned, exciseSoftClips, &bases);
    }

    // return in the orientation requested
    internal::OrientBasesAsRequested(&bases, current, orientation, impl_.IsReverseStrand(),
                                     isPulse);
    return bases;
}

Frames BamRecord::FetchFramesRaw(const BamRecordTag tag) const
{
    const Tag frameTag = impl_.TagValue(tag);
    if (frameTag.IsNull()) return {};  // throw ?

    // lossy frame codes
    if (frameTag.IsUInt8Array()) { // 这里代码认为 u8 是有损的, u16是无损的
        const auto codes = frameTag.ToUInt8Array();
        return Frames::Decode(codes);
    }

    // lossless frame data
    else {
        assert(frameTag.IsUInt16Array());
        return Frames{frameTag.ToUInt16Array()};
    }
}

Frames BamRecord::FetchFrames(const BamRecordTag tag, const Orientation orientation,
                              const bool aligned, const bool exciseSoftClips,
                              const PulseBehavior pulseBehavior) const
{
    const bool isPulse = internal::BamRecordTags::IsPulse(tag);

    // fetch raw
    Frames frames = FetchFramesRaw(tag);
    Orientation current = Orientation::NATIVE;

    // maybe strip 'squashed' pulse loci
    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        CalculatePulse2BaseCache();
        frames.DataRaw() = p2bCache_->RemoveSquashedPulses(frames.Data());
    }

    // if we need to touch the CIGAR
    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY)
            throw std::runtime_error{
                "Cannot return data at all pulses when gapping and/or soft-clipping are requested. "
                "Use PulseBehavior::BASECALLS_ONLY instead."};

        // force into genomic orientation
        internal::OrientTagDataAsRequested(&frames, current, Orientation::GENOMIC,
                                           impl_.IsReverseStrand());
        current = Orientation::GENOMIC;

        // clip & gapify as requested
        internal::ClipAndGapifyFrames(impl_, aligned, exciseSoftClips, &frames);
    }

    // return in the orientation requested
    internal::OrientTagDataAsRequested(&frames, current, orientation, impl_.IsReverseStrand());
    return frames;
}

std::vector<float> BamRecord::FetchPhotonsRaw(const BamRecordTag tag) const
{
    const Tag frameTag = impl_.TagValue(tag);
    if (frameTag.IsNull()) return {};
    if (!frameTag.IsUInt16Array())
        throw std::runtime_error{"Photons are not a uint16_t array, tag " +
                                 internal::BamRecordTags::LabelFor(tag)};

    const auto data = frameTag.ToUInt16Array();
    std::vector<float> photons;
    photons.reserve(data.size());
    for (const auto& d : data)
        photons.emplace_back(d / photonFactor);
    return photons;
}

std::vector<float> BamRecord::FetchPhotons(const BamRecordTag tag, const Orientation orientation,
                                           const bool aligned, const bool exciseSoftClips,
                                           const PulseBehavior pulseBehavior) const
{
    const bool isPulse = internal::BamRecordTags::IsPulse(tag);

    // fetch raw
    auto data = FetchPhotonsRaw(tag);
    Orientation current = Orientation::NATIVE;

    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        // strip 'squashed' pulse loci
        CalculatePulse2BaseCache();
        data = p2bCache_->RemoveSquashedPulses(data);
    }

    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY)
            throw std::runtime_error{
                "Cannot return data at all pulses when gapping and/or soft-clipping are requested. "
                "Use PulseBehavior::BASECALLS_ONLY instead."};

        // force into genomic orientation
        internal::OrientTagDataAsRequested(&data, current, Orientation::GENOMIC,
                                           impl_.IsReverseStrand());
        current = Orientation::GENOMIC;

        // clip & gapify as requested
        internal::ClipAndGapifyPhotons(impl_, aligned, exciseSoftClips, &data);
    }

    // return in the orientation requested
    internal::OrientTagDataAsRequested(&data, current, orientation, impl_.IsReverseStrand());
    return data;
}

QualityValues BamRecord::FetchQualitiesRaw(const BamRecordTag tag) const
{
    const Tag qvsTag = impl_.TagValue(tag);
    return QualityValues::FromFastq(qvsTag.ToString());
}

QualityValues BamRecord::FetchQualities(const BamRecordTag tag, const Orientation orientation,
                                        const bool aligned, const bool exciseSoftClips,
                                        const PulseBehavior pulseBehavior) const
{
    // requested data info
    const bool isBamQual = (tag == BamRecordTag::QUAL);
    const bool isPulse = internal::BamRecordTags::IsPulse(tag);

    // fetch raw
    QualityValues quals;
    Orientation current;
    if (isBamQual) {  // QUAL stored in genomic orientation
        quals = impl_.Qualities();
        current = Orientation::GENOMIC;
    } else {  // all tags stored in native orientation
        quals = FetchQualitiesRaw(tag);
        current = Orientation::NATIVE;
    }

    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        // strip 'squashed' pulse loci
        CalculatePulse2BaseCache();
        quals = p2bCache_->RemoveSquashedPulses(quals);
    }

    // if we need to touch CIGAR
    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY)
            throw std::runtime_error{
                "Cannot return data at all pulses when gapping and/or soft-clipping are requested. "
                "Use PulseBehavior::BASECALLS_ONLY instead."};

        // force into genomic orientation
        internal::OrientTagDataAsRequested(&quals, current, Orientation::GENOMIC,
                                           impl_.IsReverseStrand());
        current = Orientation::GENOMIC;

        // clip & gapify as requested
        internal::ClipAndGapifyQualities(impl_, aligned, exciseSoftClips, &quals);
    }

    // return in the orientation requested
    internal::OrientTagDataAsRequested(&quals, current, orientation, impl_.IsReverseStrand());
    return quals;
}

std::vector<uint32_t> BamRecord::FetchUInt32sRaw(const BamRecordTag tag) const
{
    // fetch tag data
    const Tag frameTag = impl_.TagValue(tag);
    if (frameTag.IsNull()) return {};
    if (!frameTag.IsUInt32Array())
        throw std::runtime_error{"Tag data are not a uint32_t array, tag " +
                                 internal::BamRecordTags::LabelFor(tag)};
    return frameTag.ToUInt32Array();
}

std::vector<uint32_t> BamRecord::FetchUInt32s(const BamRecordTag tag, const Orientation orientation,
                                              const bool aligned, const bool exciseSoftClips,
                                              const PulseBehavior pulseBehavior) const
{
    const bool isPulse = internal::BamRecordTags::IsPulse(tag);

    // fetch raw
    auto arr = FetchUInt32sRaw(tag);
    Orientation current = Orientation::NATIVE;

    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        // strip 'squashed' pulse loci
        CalculatePulse2BaseCache();
        arr = p2bCache_->RemoveSquashedPulses(arr);
    }

    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY)
            throw std::runtime_error{
                "Cannot return data at all pulses when gapping and/or soft-clipping are requested. "
                "Use PulseBehavior::BASECALLS_ONLY instead."};

        // force into genomic orientation
        internal::OrientTagDataAsRequested(&arr, current, Orientation::GENOMIC,
                                           impl_.IsReverseStrand());
        current = Orientation::GENOMIC;

        // clip & gapify as requested
        internal::ClipAndGapifyUInts(impl_, aligned, exciseSoftClips, &arr);
    }

    // return in the orientation requested
    internal::OrientTagDataAsRequested(&arr, current, orientation, impl_.IsReverseStrand());
    return arr;
}

std::vector<uint8_t> BamRecord::FetchUInt8sRaw(const BamRecordTag tag) const
{
    // fetch tag data
    const Tag frameTag = impl_.TagValue(tag);
    if (frameTag.IsNull()) return {};
    if (!frameTag.IsUInt8Array())
        throw std::runtime_error{"Tag data are not a uint8_t array, tag " +
                                 internal::BamRecordTags::LabelFor(tag)};
    return frameTag.ToUInt8Array();
}

std::vector<uint8_t> BamRecord::FetchUInt8s(const BamRecordTag tag, const Orientation orientation,
                                            const bool aligned, const bool exciseSoftClips,
                                            const PulseBehavior pulseBehavior) const
{
    const bool isPulse = internal::BamRecordTags::IsPulse(tag);

    // fetch raw
    auto arr = FetchUInt8sRaw(tag);
    Orientation current = Orientation::NATIVE;

    if (isPulse && pulseBehavior == PulseBehavior::BASECALLS_ONLY) {
        // strip 'squashed' pulse loci
        CalculatePulse2BaseCache();
        arr = p2bCache_->RemoveSquashedPulses(arr);
    }

    if (aligned || exciseSoftClips) {

        if (isPulse && pulseBehavior != PulseBehavior::BASECALLS_ONLY)
            throw std::runtime_error{
                "Cannot return data at all pulses when gapping and/or soft-clipping are requested. "
                "Use PulseBehavior::BASECALLS_ONLY instead."};

        // force into genomic orientation
        internal::OrientTagDataAsRequested(&arr, current, Orientation::GENOMIC,
                                           impl_.IsReverseStrand());
        current = Orientation::GENOMIC;

        // clip & gapify as requested
        internal::ClipAndGapifyUInt8s(impl_, aligned, exciseSoftClips, &arr);
    }

    // return in the orientation requested
    internal::OrientTagDataAsRequested(&arr, current, orientation, impl_.IsReverseStrand());
    return arr;
}

std::string BamRecord::FullName() const { return impl_.Name(); }

bool BamRecord::HasAltLabelQV() const { return impl_.HasTag(BamRecordTag::ALT_LABEL_QV); }

bool BamRecord::HasAltLabelTag() const { return impl_.HasTag(BamRecordTag::ALT_LABEL_TAG); }

bool BamRecord::HasBarcodes() const { return impl_.HasTag(BamRecordTag::BARCODES); }

bool BamRecord::HasBarcodeQuality() const { return impl_.HasTag(BamRecordTag::BARCODE_QUALITY); }

bool BamRecord::HasLabelQV() const { return impl_.HasTag(BamRecordTag::LABEL_QV); }

bool BamRecord::HasDeletionQV() const { return impl_.HasTag(BamRecordTag::DELETION_QV); }

bool BamRecord::HasDeletionTag() const { return impl_.HasTag(BamRecordTag::DELETION_TAG); }

bool BamRecord::HasHoleNumber() const
{
    return impl_.HasTag(BamRecordTag::HOLE_NUMBER) &&
           !impl_.TagValue(BamRecordTag::HOLE_NUMBER).IsNull();
}

bool BamRecord::HasInsertionQV() const { return impl_.HasTag(BamRecordTag::INSERTION_QV); }

bool BamRecord::HasNumPasses() const { return impl_.HasTag(BamRecordTag::NUM_PASSES); }

bool BamRecord::HasPreBaseFrames() const { return HasIPD(); }

bool BamRecord::HasIPD() const { return impl_.HasTag(BamRecordTag::IPD); }

bool BamRecord::HasCR() const { return impl_.HasTag(BamRecordTag::CAPTURE_RATE); }

bool BamRecord::HasLocalContextFlags() const { return impl_.HasTag(BamRecordTag::CONTEXT_FLAGS); }

bool BamRecord::HasMergeQV() const { return impl_.HasTag(BamRecordTag::MERGE_QV); }

bool BamRecord::HasPulseMergeQV() const { return impl_.HasTag(BamRecordTag::PULSE_MERGE_QV); }

bool BamRecord::HasPkmean() const { return impl_.HasTag(BamRecordTag::PKMEAN); }

bool BamRecord::HasPkmean2() const { return impl_.HasTag(BamRecordTag::PKMEAN_2); }

bool BamRecord::HasPkmid() const { return impl_.HasTag(BamRecordTag::PKMID); }

bool BamRecord::HasPkmid2() const { return impl_.HasTag(BamRecordTag::PKMID_2); }

bool BamRecord::HasPrePulseFrames() const { return impl_.HasTag(BamRecordTag::PRE_PULSE_FRAMES); }

bool BamRecord::HasPulseCall() const
{
    return impl_.HasTag(BamRecordTag::PULSE_CALL) &&
           !impl_.TagValue(BamRecordTag::PULSE_CALL).IsNull();
}

bool BamRecord::HasPulseExclusion(void) const
{
    return impl_.HasTag(BamRecordTag::PULSE_EXCLUSION);
}

bool BamRecord::HasPulseCallWidth(void) const
{
    return impl_.HasTag(BamRecordTag::PULSE_CALL_WIDTH);
}

bool BamRecord::HasPulseWidth() const { return impl_.HasTag(BamRecordTag::PULSE_WIDTH); }

bool BamRecord::HasQueryEnd() const { return impl_.HasTag(BamRecordTag::QUERY_END); }

bool BamRecord::HasQueryStart() const { return impl_.HasTag(BamRecordTag::QUERY_START); }

bool BamRecord::HasReadAccuracy() const
{
    return impl_.HasTag(BamRecordTag::READ_ACCURACY) &&
           !impl_.TagValue(BamRecordTag::READ_ACCURACY).IsNull();
}

bool BamRecord::HasScrapRegionType() const
{
    return impl_.HasTag(BamRecordTag::SCRAP_REGION_TYPE) &&
           !impl_.TagValue(BamRecordTag::SCRAP_REGION_TYPE).IsNull();
}

bool BamRecord::HasScrapZmwType() const
{
    return impl_.HasTag(BamRecordTag::SCRAP_ZMW_TYPE) &&
           !impl_.TagValue(BamRecordTag::SCRAP_ZMW_TYPE).IsNull();
}

bool BamRecord::HasStartFrame() const { return impl_.HasTag(BamRecordTag::START_FRAME); }

bool BamRecord::HasSignalToNoise() const { return impl_.HasTag(BamRecordTag::SNR); }

bool BamRecord::HasSubstitutionQV() const { return impl_.HasTag(BamRecordTag::SUBSTITUTION_QV); }

bool BamRecord::HasSubstitutionTag() const { return impl_.HasTag(BamRecordTag::SUBSTITUTION_TAG); }

BamHeader BamRecord::Header() const { return header_; }

int32_t BamRecord::HoleNumber() const
{
    const Tag holeNumber = impl_.TagValue(BamRecordTag::HOLE_NUMBER);
    if (!holeNumber.IsNull()) return holeNumber.ToInt32();

    // missing zm tag - try to pull from name
    return internal::HoleNumberFromName(FullName());
}

BamRecord& BamRecord::HoleNumber(const int32_t holeNumber)
{
    internal::CreateOrEdit(BamRecordTag::HOLE_NUMBER, holeNumber, &impl_);
    return *this;
}

BamRecordImpl& BamRecord::Impl() { return impl_; }

const BamRecordImpl& BamRecord::Impl() const { return impl_; }

QualityValues BamRecord::InsertionQV(Orientation orientation, bool aligned,
                                     bool exciseSoftClips) const
{
    return FetchQualities(BamRecordTag::INSERTION_QV, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::InsertionQV(const QualityValues& insertionQVs)
{
    internal::CreateOrEdit(BamRecordTag::INSERTION_QV, insertionQVs.Fastq(), &impl_);
    return *this;
}

Frames BamRecord::IPD(Orientation orientation, bool aligned, bool exciseSoftClips) const
{
    return FetchFrames(BamRecordTag::IPD, orientation, aligned, exciseSoftClips);
}

Frames BamRecord::CR(Orientation orientation, bool aligned, bool exciseSoftClips) const
{
    return FetchFrames(BamRecordTag::CAPTURE_RATE, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::IPD(const Frames& frames, const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY)
        internal::CreateOrEdit(BamRecordTag::IPD, frames.Encode(), &impl_);
    else
        internal::CreateOrEdit(BamRecordTag::IPD, frames.Data(), &impl_);
    return *this;
}

Frames BamRecord::IPDRaw(Orientation orientation) const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::IPD);
    const Tag frameTag = impl_.TagValue(tagName);
    if (frameTag.IsNull()) return {};

    Frames frames;

    // lossy frame codes
    if (frameTag.IsUInt8Array()) {
        const auto codes = frameTag.ToUInt8Array();
        const std::vector<uint16_t> codes16(codes.begin(), codes.end());
        frames.Data(std::move(codes16));
    }

    // lossless frame data
    else {
        assert(frameTag.IsUInt16Array());
        frames.Data(frameTag.ToUInt16Array());
    }

    // return in requested orientation
    internal::OrientTagDataAsRequested(&frames,
                                       Orientation::NATIVE,  // current
                                       orientation,          // requested
                                       impl_.IsReverseStrand());
    return frames;
}

bool BamRecord::IsMapped() const { return impl_.IsMapped(); }

QualityValues BamRecord::LabelQV(Orientation orientation, bool aligned, bool exciseSoftClips,
                                 PulseBehavior pulseBehavior) const
{
    return FetchQualities(BamRecordTag::LABEL_QV, orientation, aligned, exciseSoftClips,
                          pulseBehavior);
}

BamRecord& BamRecord::LabelQV(const QualityValues& labelQVs)
{
    internal::CreateOrEdit(BamRecordTag::LABEL_QV, labelQVs.Fastq(), &impl_);
    return *this;
}

LocalContextFlags BamRecord::LocalContextFlags() const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::CONTEXT_FLAGS);
    const Tag cxTag = impl_.TagValue(tagName);
    return static_cast<PacBio::BAM::LocalContextFlags>(cxTag.ToUInt8());
}

BamRecord& BamRecord::LocalContextFlags(const PacBio::BAM::LocalContextFlags flags)
{
    internal::CreateOrEdit(BamRecordTag::CONTEXT_FLAGS, static_cast<uint8_t>(flags), &impl_);
    return *this;
}

BamRecord& BamRecord::Map(const int32_t referenceId, const Position refStart, const Strand strand,
                          const Cigar& cigar, const uint8_t mappingQuality)
{
    impl_.Position(refStart);
    impl_.ReferenceId(referenceId);
    impl_.CigarData(cigar);
    impl_.MapQuality(mappingQuality);
    impl_.SetMapped(true);

    if (strand == Strand::FORWARD)
        impl_.SetReverseStrand(false);

    else {
        assert(strand == Strand::REVERSE);
        impl_.SetReverseStrand(true);

        // switch seq & qual
        std::string sequence = impl_.Sequence();
        QualityValues qualities = impl_.Qualities();

        internal::ReverseComplement(sequence);
        internal::Reverse(qualities);

        impl_.SetSequenceAndQualities(sequence, qualities.Fastq());
    }

    // reset any cached aligned start/end
    alignedStart_ = PacBio::BAM::UnmappedPosition;
    alignedEnd_ = PacBio::BAM::UnmappedPosition;

    return *this;
}

uint8_t BamRecord::MapQuality() const { return impl_.MapQuality(); }

QualityValues BamRecord::MergeQV(Orientation orientation, bool aligned, bool exciseSoftClips) const
{
    return FetchQualities(BamRecordTag::MERGE_QV, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::MergeQV(const QualityValues& mergeQVs)
{
    internal::CreateOrEdit(BamRecordTag::MERGE_QV, mergeQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::MovieName() const { return ReadGroup().MovieName(); }

size_t BamRecord::NumDeletedBases() const
{
    auto tEnd = ReferenceEnd();
    auto tStart = ReferenceStart();
    auto numMatchesAndMismatches = NumMatchesAndMismatches();
    auto nM = numMatchesAndMismatches.first;
    auto nMM = numMatchesAndMismatches.second;
    return (tEnd - tStart - nM - nMM);
}

size_t BamRecord::NumInsertedBases() const
{
    auto aEnd = AlignedEnd();
    auto aStart = AlignedStart();
    auto numMatchesAndMismatches = NumMatchesAndMismatches();
    auto nM = numMatchesAndMismatches.first;
    auto nMM = numMatchesAndMismatches.second;
    return (aEnd - aStart - nM - nMM);
}

size_t BamRecord::NumMatches() const { return NumMatchesAndMismatches().first; }

std::pair<size_t, size_t> BamRecord::NumMatchesAndMismatches() const
{
    std::pair<size_t, size_t> result = std::make_pair(0, 0);

    auto b = internal::BamRecordMemory::GetRawData(this);
    uint32_t* cigarData = bam_get_cigar(b.get());
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        const auto type = static_cast<CigarOperationType>(bam_cigar_op(cigarData[i]));
        if (type == CigarOperationType::SEQUENCE_MATCH)
            result.first += bam_cigar_oplen(cigarData[i]);
        else if (type == CigarOperationType::SEQUENCE_MISMATCH)
            result.second += bam_cigar_oplen(cigarData[i]);
    }
    return result;
}

size_t BamRecord::NumMismatches() const { return NumMatchesAndMismatches().second; }

int32_t BamRecord::NumPasses() const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::NUM_PASSES);
    const Tag numPasses = impl_.TagValue(tagName);
    return numPasses.ToInt32();
}

BamRecord& BamRecord::NumPasses(const int32_t numPasses)
{
    internal::CreateOrEdit(BamRecordTag::NUM_PASSES, numPasses, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmean(Orientation orientation, bool aligned, bool exciseSoftClips,
                                     PulseBehavior pulseBehavior) const
{
    return FetchPhotons(BamRecordTag::PKMEAN, orientation, aligned, exciseSoftClips, pulseBehavior);
}

BamRecord& BamRecord::Pkmean(const std::vector<float>& photons)
{
    Pkmean(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmean(const std::vector<uint16_t>& encodedPhotons)
{
    internal::CreateOrEdit(BamRecordTag::PKMEAN, encodedPhotons, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmid(Orientation orientation, bool aligned, bool exciseSoftClips,
                                    PulseBehavior pulseBehavior) const
{
    return FetchPhotons(BamRecordTag::PKMID, orientation, aligned, exciseSoftClips, pulseBehavior);
}

BamRecord& BamRecord::Pkmid(const std::vector<float>& photons)
{
    Pkmid(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmid(const std::vector<uint16_t>& encodedPhotons)
{
    internal::CreateOrEdit(BamRecordTag::PKMID, encodedPhotons, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmean2(Orientation orientation, bool aligned, bool exciseSoftClips,
                                      PulseBehavior pulseBehavior) const
{
    return FetchPhotons(BamRecordTag::PKMEAN_2, orientation, aligned, exciseSoftClips,
                        pulseBehavior);
}

BamRecord& BamRecord::Pkmean2(const std::vector<float>& photons)
{
    Pkmean2(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmean2(const std::vector<uint16_t>& encodedPhotons)
{
    internal::CreateOrEdit(BamRecordTag::PKMEAN_2, encodedPhotons, &impl_);
    return *this;
}

std::vector<float> BamRecord::Pkmid2(Orientation orientation, bool aligned, bool exciseSoftClips,
                                     PulseBehavior pulseBehavior) const
{
    return FetchPhotons(BamRecordTag::PKMID_2, orientation, aligned, exciseSoftClips,
                        pulseBehavior);
}

BamRecord& BamRecord::Pkmid2(const std::vector<float>& photons)
{
    Pkmid2(EncodePhotons(photons));
    return *this;
}

BamRecord& BamRecord::Pkmid2(const std::vector<uint16_t>& encodedPhotons)
{
    internal::CreateOrEdit(BamRecordTag::PKMID_2, encodedPhotons, &impl_);
    return *this;
}

Frames BamRecord::PreBaseFrames(Orientation orientation, bool aligned, bool exciseSoftClips) const
{
    return IPD(orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::PreBaseFrames(const Frames& frames, const FrameEncodingType encoding)
{
    return IPD(frames, encoding);
}

Frames BamRecord::PrePulseFrames(Orientation orientation, bool aligned, bool exciseSoftClips,
                                 PulseBehavior pulseBehavior) const
{
    return FetchFrames(BamRecordTag::PRE_PULSE_FRAMES, orientation, aligned, exciseSoftClips,
                       pulseBehavior);
}

BamRecord& BamRecord::PrePulseFrames(const Frames& frames, const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY) {
        internal::CreateOrEdit(BamRecordTag::PRE_PULSE_FRAMES, frames.Encode(), &impl_);
    } else {
        internal::CreateOrEdit(BamRecordTag::PRE_PULSE_FRAMES, frames.Data(), &impl_);
    }
    return *this;
}

Frames BamRecord::PulseWidthRaw(Orientation orientation, bool /* aligned */,
                                bool /* exciseSoftClips */) const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::PULSE_WIDTH);
    const Tag frameTag = impl_.TagValue(tagName);
    if (frameTag.IsNull()) return {};

    Frames frames;

    // lossy frame codes
    if (frameTag.IsUInt8Array()) {
        const auto codes = frameTag.ToUInt8Array();
        const std::vector<uint16_t> codes16(codes.begin(), codes.end());
        frames.Data(std::move(codes16));
    }

    // lossless frame data
    else {
        assert(frameTag.IsUInt16Array());
        frames.Data(frameTag.ToUInt16Array());
    }

    // return in requested orientation
    internal::OrientTagDataAsRequested(&frames,
                                       Orientation::NATIVE,  // current
                                       orientation,          // requested
                                       impl_.IsReverseStrand());
    return frames;
}

QualityValues BamRecord::PulseMergeQV(Orientation orientation, bool aligned, bool exciseSoftClips,
                                      PulseBehavior pulseBehavior) const
{
    return FetchQualities(BamRecordTag::PULSE_MERGE_QV, orientation, aligned, exciseSoftClips,
                          pulseBehavior);
}

BamRecord& BamRecord::PulseMergeQV(const QualityValues& mergeQVs)
{
    internal::CreateOrEdit(BamRecordTag::PULSE_MERGE_QV, mergeQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::PulseCall(Orientation orientation, bool aligned, bool exciseSoftClips,
                                 PulseBehavior pulseBehavior) const
{
    return FetchBases(BamRecordTag::PULSE_CALL, orientation, aligned, exciseSoftClips,
                      pulseBehavior);
}

BamRecord& BamRecord::PulseCall(const std::string& tags)
{
    internal::CreateOrEdit(BamRecordTag::PULSE_CALL, tags, &impl_);
    return *this;
}

Frames BamRecord::PulseCallWidth(Orientation orientation, bool aligned, bool exciseSoftClips,
                                 PulseBehavior pulseBehavior) const
{
    return FetchFrames(BamRecordTag::PULSE_CALL_WIDTH, orientation, aligned, exciseSoftClips,
                       pulseBehavior);
}

BamRecord& BamRecord::PulseCallWidth(const Frames& frames, const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY) {
        internal::CreateOrEdit(BamRecordTag::PULSE_CALL_WIDTH, frames.Encode(), &impl_);
    } else {
        internal::CreateOrEdit(BamRecordTag::PULSE_CALL_WIDTH, frames.Data(), &impl_);
    }
    return *this;
}

std::vector<PacBio::BAM::PulseExclusionReason> BamRecord::PulseExclusionReason(
    Orientation orientation, bool aligned, bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    std::vector<PacBio::BAM::PulseExclusionReason> reasons;

    const auto reasonNums = FetchUInt8s(BamRecordTag::PULSE_EXCLUSION, orientation, aligned,
                                        exciseSoftClips, pulseBehavior);

    std::transform(
        reasonNums.cbegin(), reasonNums.cend(), std::back_inserter(reasons),
        [](const uint8_t num) { return static_cast<PacBio::BAM::PulseExclusionReason>(num); });

    return reasons;
}

BamRecord& BamRecord::PulseExclusionReason(
    const std::vector<PacBio::BAM::PulseExclusionReason>& reasons)
{
    std::vector<uint8_t> reasonNums;
    std::transform(reasons.cbegin(), reasons.cend(), std::back_inserter(reasonNums),
                   [](const PacBio::BAM::PulseExclusionReason& reason) {
                       return static_cast<uint8_t>(reason);
                   });

    internal::CreateOrEdit(BamRecordTag::PULSE_EXCLUSION, reasonNums, &impl_);
    return *this;
}

Frames BamRecord::PulseWidth(Orientation orientation, bool aligned, bool exciseSoftClips) const
{
    return FetchFrames(BamRecordTag::PULSE_WIDTH, orientation, aligned, exciseSoftClips,
                       PulseBehavior::ALL);
}

BamRecord& BamRecord::PulseWidth(const Frames& frames, const FrameEncodingType encoding)
{
    if (encoding == FrameEncodingType::LOSSY) {
        internal::CreateOrEdit(BamRecordTag::PULSE_WIDTH, frames.Encode(), &impl_);
    } else {
        internal::CreateOrEdit(BamRecordTag::PULSE_WIDTH, frames.Data(), &impl_);
    }
    return *this;
}

QualityValues BamRecord::Qualities(Orientation orientation, bool aligned,
                                   bool exciseSoftClips) const
{
    return FetchQualities(BamRecordTag::QUAL, orientation, aligned, exciseSoftClips);
}

Position BamRecord::QueryEnd() const
{
    // try 'qe' tag
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::QUERY_END);
    const Tag qe = impl_.TagValue(tagName);
    if (!qe.IsNull()) return qe.ToInt32();

    // tag missing, need to check movie name (fallback for non-PB BAMs, but ignore for CCS reads)
    RecordType type;
    try {
        type = Type();
    } catch (std::exception&) {
        return 0;
    }
    if (type == RecordType::CCS) throw std::runtime_error{"no query end for CCS read type"};
    if (type == RecordType::TRANSCRIPT)
        throw std::runtime_error{"no query end for transcript read type"};

    // PacBio BAM, non-CCS/transcript
    try {
        return internal::QueryEndFromName(FullName());
    } catch (std::exception&) {
        // return fallback position
        return 0;
    }
}

BamRecord& BamRecord::QueryEnd(const Position pos)
{
    internal::CreateOrEdit(BamRecordTag::QUERY_END, static_cast<int32_t>(pos), &impl_);
    UpdateName();
    return *this;
}

Position BamRecord::QueryStart() const
{
    // try 'qs' tag
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::QUERY_START);
    const Tag qs = impl_.TagValue(tagName);
    if (!qs.IsNull()) return qs.ToInt32();

    // tag missing, need to check movie name (fallback for non-PB BAMs, but ignore for CCS reads)
    RecordType type;
    try {
        type = Type();
    } catch (std::exception&) {
        return 0;
    }
    if (type == RecordType::CCS) throw std::runtime_error{"no query start for CCS read type"};
    if (type == RecordType::TRANSCRIPT)
        throw std::runtime_error{"no query start for transcript read type"};

    // PacBio BAM, non-CCS/transcript
    try {
        return internal::QueryStartFromName(FullName());
    } catch (std::exception&) {
        // return fallback position
        return 0;
    }
}

BamRecord& BamRecord::QueryStart(const Position pos)
{
    internal::CreateOrEdit(BamRecordTag::QUERY_START, static_cast<int32_t>(pos), &impl_);
    UpdateName();
    return *this;
}

Accuracy BamRecord::ReadAccuracy() const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::READ_ACCURACY);
    const Tag readAccuracy = impl_.TagValue(tagName);
    return {readAccuracy.ToFloat()};
}

BamRecord& BamRecord::ReadAccuracy(const Accuracy& accuracy)
{
    internal::CreateOrEdit(BamRecordTag::READ_ACCURACY, static_cast<float>(accuracy), &impl_);
    return *this;
}

ReadGroupInfo BamRecord::ReadGroup() const { return header_.ReadGroup(ReadGroupId()); }

BamRecord& BamRecord::ReadGroup(const ReadGroupInfo& rg)
{
    internal::CreateOrEdit(BamRecordTag::READ_GROUP, rg.Id(), &impl_);
    UpdateName();
    return *this;
}

std::string BamRecord::ReadGroupId() const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::READ_GROUP);
    const Tag rgTag = impl_.TagValue(tagName);
    if (rgTag.IsNull()) return {};
    return rgTag.ToString();
}

BamRecord& BamRecord::ReadGroupId(const std::string& id)
{
    internal::CreateOrEdit(BamRecordTag::READ_GROUP, id, &impl_);
    UpdateName();
    return *this;
}

int32_t BamRecord::ReadGroupNumericId() const { return ReadGroupInfo::IdToInt(ReadGroupId()); }

Position BamRecord::ReferenceEnd() const
{
    if (!impl_.IsMapped()) return PacBio::BAM::UnmappedPosition;
    const auto htsData = internal::BamRecordMemory::GetRawData(impl_);
    if (!htsData) return PacBio::BAM::UnmappedPosition;
    return bam_endpos(htsData.get());
}

int32_t BamRecord::ReferenceId() const { return impl_.ReferenceId(); }

std::string BamRecord::ReferenceName() const
{
    if (IsMapped())
        return Header().SequenceName(ReferenceId());
    else
        throw std::runtime_error{"unmapped record has no associated reference name"};
}

Position BamRecord::ReferenceStart() const { return impl_.Position(); }

void BamRecord::ResetCachedPositions() const
{
    alignedEnd_ = PacBio::BAM::UnmappedPosition;
    alignedStart_ = PacBio::BAM::UnmappedPosition;
}

void BamRecord::ResetCachedPositions()
{
    alignedEnd_ = PacBio::BAM::UnmappedPosition;
    alignedStart_ = PacBio::BAM::UnmappedPosition;
}

VirtualRegionType BamRecord::ScrapRegionType() const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::SCRAP_REGION_TYPE);
    const Tag srTag = impl_.TagValue(tagName);
    return VirtualRegionTypeMap::ParseChar[srTag.ToUInt8()];
}

BamRecord& BamRecord::ScrapRegionType(const VirtualRegionType type)
{
    internal::CreateOrEdit(BamRecordTag::SCRAP_REGION_TYPE, static_cast<uint8_t>(type), &impl_);
    return *this;
}

BamRecord& BamRecord::ScrapRegionType(const char type)
{
    internal::CreateOrEdit(BamRecordTag::SCRAP_REGION_TYPE, type, &impl_);
    return *this;
}

ZmwType BamRecord::ScrapZmwType() const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::SCRAP_ZMW_TYPE);
    const Tag szTag = impl_.TagValue(tagName);
    return ZmwTypeMap::ParseChar[szTag.ToUInt8()];
}

BamRecord& BamRecord::ScrapZmwType(const ZmwType type)
{
    internal::CreateOrEdit(BamRecordTag::SCRAP_ZMW_TYPE, static_cast<uint8_t>(type), &impl_);
    return *this;
}

BamRecord& BamRecord::ScrapZmwType(const char type)
{
    internal::CreateOrEdit(BamRecordTag::SCRAP_ZMW_TYPE, type, &impl_);
    return *this;
}

std::string BamRecord::Sequence(const Orientation orientation, bool aligned,
                                bool exciseSoftClips) const
{
    return FetchBases(BamRecordTag::SEQ, orientation, aligned, exciseSoftClips);
}

std::vector<float> BamRecord::SignalToNoise() const
{
    const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::SNR);
    const Tag snTag = impl_.TagValue(tagName);
    return snTag.ToFloatArray();
}

std::vector<float> BamRecord::BaseLevelCrs() const {
    if (impl_.HasTag(BamRecordTag::BASE_LEVEL_CRS)) {
        const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::BASE_LEVEL_CRS);
        const Tag snTag = impl_.TagValue(tagName);
        return snTag.ToFloatArray();
    }
    return {-1, -1, -1, -1};
    
}
std::vector<float> BamRecord::BaseLevelDwMeans() const {
    if (impl_.HasTag(BamRecordTag::BASE_LEVEL_DW_MEANS)) {
        const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::BASE_LEVEL_DW_MEANS);
        const Tag snTag = impl_.TagValue(tagName);
        return snTag.ToFloatArray();
    }
    return {-1, -1, -1, -1};
}
std::vector<float> BamRecord::BaseLevelDwVars() const {
    if (impl_.HasTag(BamRecordTag::BASE_LEVEL_DW_VARS)) {
        const auto tagName = internal::BamRecordTags::LabelFor(BamRecordTag::BASE_LEVEL_DW_VARS);
        const Tag snTag = impl_.TagValue(tagName);
        return snTag.ToFloatArray();
    }
    return {-1, -1, -1, -1};
}

BamRecord& BamRecord::SignalToNoise(const std::vector<float>& snr)
{
    internal::CreateOrEdit(BamRecordTag::SNR, snr, &impl_);
    return *this;
}

std::vector<uint32_t> BamRecord::StartFrame(Orientation orientation, bool aligned,
                                            bool exciseSoftClips, PulseBehavior pulseBehavior) const
{
    return FetchUInt32s(BamRecordTag::START_FRAME, orientation, aligned, exciseSoftClips,
                        pulseBehavior);
}

BamRecord& BamRecord::StartFrame(const std::vector<uint32_t>& startFrame)
{
    internal::CreateOrEdit(BamRecordTag::START_FRAME, startFrame, &impl_);
    return *this;
}

QualityValues BamRecord::SubstitutionQV(Orientation orientation, bool aligned,
                                        bool exciseSoftClips) const
{
    return FetchQualities(BamRecordTag::SUBSTITUTION_QV, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::SubstitutionQV(const QualityValues& substitutionQVs)
{
    internal::CreateOrEdit(BamRecordTag::SUBSTITUTION_QV, substitutionQVs.Fastq(), &impl_);
    return *this;
}

std::string BamRecord::SubstitutionTag(Orientation orientation, bool aligned,
                                       bool exciseSoftClips) const
{
    return FetchBases(BamRecordTag::SUBSTITUTION_TAG, orientation, aligned, exciseSoftClips);
}

BamRecord& BamRecord::SubstitutionTag(const std::string& tags)
{
    internal::CreateOrEdit(BamRecordTag::SUBSTITUTION_TAG, tags, &impl_);
    return *this;
}

RecordType BamRecord::Type() const
{
    try {
        const auto typeName = ReadGroup().ReadType();
        return internal::NameToType(typeName);
    } catch (std::exception&) {

        // read group not found, peek at name to see if we're possibly
        // CCS or TRANSCRIPT
        //
        const auto name = FullName();
        if (name.find("transcript") == 0)
            return RecordType::TRANSCRIPT;
        else if (name.find("/ccs") != std::string::npos)
            return RecordType::CCS;
        else
            return RecordType::UNKNOWN;
    }
}

void BamRecord::UpdateName()
{
    std::string newName;
    newName.reserve(100);

    const auto holeNumber = (HasHoleNumber() ? std::to_string(HoleNumber()) : "?");
    if (Type() == RecordType::TRANSCRIPT) {
        newName = "transcript/" + holeNumber;
    } else {
        newName += MovieName();
        newName += "/";
        newName += holeNumber;
        newName += "/";

        if (Type() == RecordType::CCS)
            newName += "ccs";

        else {
            if (HasQueryStart())
                newName += std::to_string(QueryStart());
            else
                newName += "?";

            newName += '_';

            if (HasQueryEnd())
                newName += std::to_string(QueryEnd());
            else
                newName += "?";
        }
    }
    impl_.Name(newName);
}

}  // namespace BAM
}  // namespace PacBio
