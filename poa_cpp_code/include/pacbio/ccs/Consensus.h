// Author: Lance Hepler

#pragma once

#include <pacbio/ccs/ConsensusSettings.h>
#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Polish.h>
#include <pacbio/data/ReadId.h>
#include <pacbio/data/Sequence.h>
#include <pacbio/data/State.h>
#include <pacbio/data/StrandType.h>
#include <pacbio/data/SubreadResultCounter.h>
#include <pacbio/denovo/PoaConsensus.h>
#include <pacbio/denovo/SparsePoa.h>
#include <pacbio/util/Timer.h>
#include <pbbam/Accuracy.h>
#include <pbbam/LocalContextFlags.h>
#include <pbcopper/logging/Logging.h>

#include <algorithm>
#include <boost/optional.hpp>
#include <cmath>
#include <functional>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "./EdAlignment/EdlibGlobalAlign.h"

namespace PacBio {
namespace CCS {

using Accuracy = PacBio::BAM::Accuracy;
using Interval = PacBio::Data::Interval;
using LocalContextFlags = PacBio::BAM::LocalContextFlags;
using QualityValues = PacBio::Consensus::QualityValues;
using ReadId = PacBio::Data::ReadId;
using Read = PacBio::Data::Read;
using MappedRead = PacBio::Data::MappedRead;
using SNR = PacBio::Data::SNR;
using State = PacBio::Data::State;
using StrandType = PacBio::Data::StrandType;
using SubreadResultCounter = PacBio::Data::SubreadResultCounter;
using Timer = PacBio::Util::Timer;
using PoaAlignmentSummary = PacBio::Poa::PoaAlignmentSummary;
using SparsePoa = PacBio::Poa::SparsePoa;
using AlignConfig = PacBio::Align::AlignConfig;
using AlignMode = PacBio::Align::AlignMode;

const int BASE_Q_MAX = 45;
const int READ_Q_MAX = 30;
const int READ_Q_MAX_BASE = 30;

template <typename TId>
struct ReadType
{
    TId Id;
    std::string Seq;
    std::vector<uint8_t> IPD;
    std::vector<uint8_t> PulseWidth;
    LocalContextFlags Flags;
    Accuracy ReadAccuracy;
    SNR SignalToNoise;
    Data::NucleotideLevelFeats Feats;
    std::string Chemistry;
    std::vector<uint8_t> CR;
};

template <typename TId, typename TRead>
struct ChunkType
{
    TId Id;
    std::vector<TRead> Reads;
    boost::optional<std::tuple<int16_t, int16_t, uint8_t>> Barcodes;
};

struct ConsensusType
{
    Consensus::PolishResult polishResult;
    ReadId Id;
    boost::optional<StrandType> Strand;
    std::string Sequence;
    QualityValues QVs;
    size_t NumPasses;
    double PredictedAccuracy;
    double AvgZScore;
    std::vector<double> ZScores;
    std::vector<int32_t> StatusCounts;
    float ElapsedMilliseconds;
    boost::optional<SNR> SignalToNoise;
    boost::optional<std::tuple<int16_t, int16_t, uint8_t>> Barcodes;
};

template <typename TConsensus>
class ResultType : public std::vector<TConsensus>
{
public:
    size_t Success;
    size_t PoorSNR;
    size_t NoSubreads;
    size_t TooLong;
    size_t TooShort;
    size_t TooFewPasses1;
    size_t TooFewPasses2;
    size_t TooManyUnusable;
    size_t NonConvergent;
    size_t PoorQuality;
    size_t ExceptionThrown;
    SubreadResultCounter SubreadCounter;

    ResultType()
        : Success{0}
        , PoorSNR{0}
        , NoSubreads{0}
        , TooLong{0}
        , TooShort{0}
        , TooFewPasses1{0}
        , TooFewPasses2{0}
        , TooManyUnusable{0}
        , NonConvergent{0}
        , PoorQuality{0}
        , ExceptionThrown{0}
        , SubreadCounter{}
    {
    }

    ResultType<TConsensus>& operator+=(const ResultType<TConsensus>& other)
    {
        Success += other.Success;
        PoorSNR += other.PoorSNR;
        NoSubreads += other.NoSubreads;
        TooLong += other.TooLong;
        TooShort += other.TooShort;
        TooManyUnusable += other.TooManyUnusable;
        TooFewPasses1 += other.TooFewPasses1;
        TooFewPasses2 += other.TooFewPasses2;
        NonConvergent += other.NonConvergent;
        PoorQuality += other.PoorQuality;
        ExceptionThrown += other.ExceptionThrown;
        SubreadCounter += other.SubreadCounter;
        return *this;
    }

    size_t Total() const
    {
        return (Success + PoorSNR + NoSubreads + TooShort + TooManyUnusable + TooFewPasses1 +
                TooFewPasses2 + NonConvergent + PoorQuality + ExceptionThrown);
    }
};

namespace {  // anonymous

size_t ComputePassesCnt(BAM::LocalContextFlags flag, size_t cnt, const ConsensusSettings& settings)
{
    if (settings.PassesCntVersion == "v1") {
        if (flag & BAM::ADAPTER_BEFORE && flag & BAM::ADAPTER_AFTER) {
            cnt += 1;
        }
    } else if (settings.PassesCntVersion == "v2") {
        cnt += 1;

    } else {
        PBLOG_FATAL << "invalid PassesCntVersion, exptected: [v1,v2]. but got "
                    << settings.PassesCntVersion;
        exit(EXIT_FAILURE);
    }
    return cnt;
}

template <typename TRead>
size_t ComputeNPassesForNoPolish(std::vector<const TRead*>& reads,
                                 std::vector<SparsePoa::ReadKey>& readKeys,
                                 std::vector<PoaAlignmentSummary>& summaries, float minIdentity,
                                 ResultType<ConsensusType>& result)
{

    const size_t nReads = readKeys.size();
    size_t nPasses = 0;
    for (int i = 0; i < nReads; i++) {
        if (readKeys[i] < 0) continue;
        if (summaries[readKeys[i]].AlignmentIdentity < minIdentity) {
            result.SubreadCounter.PoorIdentity += 1;
            continue;
        }

        result.SubreadCounter.Success += 1;

        BAM::LocalContextFlags flags = reads[i]->Flags;
        if (!(flags & BAM::ADAPTER_BEFORE && flags & BAM::ADAPTER_AFTER)) continue;
        nPasses += 1;
    }
    return nPasses;
}

template <typename T>
float Median(std::vector<T>* vs)
{
    const size_t n = vs->size();
    std::sort(vs->begin(), vs->end());
    if (n % 2 == 1) return static_cast<float>(vs->at(n / 2));
    return 0.5 * (vs->at(n / 2 - 1) + vs->at(n / 2));
}

double err_rate_from_qv(std::vector<int>& qualities)
{
    double err_rate = 0;
    int seq_len = qualities.size();
    for (const int qv : qualities) {
        err_rate += pow(10.0, static_cast<double>(qv) / -10.0);
    }
    err_rate /= seq_len;
    return err_rate;
}

double qv_calibration(std::vector<int>& qualities, int _)
{

    double err_rate = err_rate_from_qv(qualities);
    int seq_len = qualities.size();

    double min_err_rate = 0.1 / seq_len;
    double channelq = -10 * log10(err_rate);
    double max_channelq = -10 * log10(min_err_rate);
    double res_channelq = channelq > max_channelq ? max_channelq : channelq;
    double max_baseq = res_channelq + 20;

    double delta = channelq - max_channelq;
    delta = delta >= 0 ? delta : 0;

    for (int& qv : qualities) {
        qv -= delta;
        qv = qv >= 0 ? qv : 0;
        qv = qv > max_baseq ? max_baseq : qv;
        qv = qv > 50 ? 50 : qv;
    }

    err_rate = err_rate_from_qv(qualities);

    return 1 - err_rate;
}

double qv_calibration2(std::vector<int>& qualities, int _)
{

    for (int& qv : qualities) {
        qv = qv > BASE_Q_MAX ? BASE_Q_MAX : qv;
    }

    double err_rate = err_rate_from_qv(qualities);

    return 1 - err_rate;
}

double qv_calibration3(std::vector<int>& qualities, int nPasses)
{

    int predicted_max_read_q = (nPasses + 1) / 2.0 * 10.0;
    predicted_max_read_q = predicted_max_read_q > READ_Q_MAX ? READ_Q_MAX : predicted_max_read_q;
    double err_rate = err_rate_from_qv(qualities);
    double channelq = -10 * log10(err_rate);

    double delta = channelq - predicted_max_read_q;

    delta = delta >= 0 ? delta : 0;

    for (int& qv : qualities) {
        qv -= delta;
        qv = qv >= 0 ? qv : 0;
        qv = qv > BASE_Q_MAX ? BASE_Q_MAX : qv;
    }

    err_rate = err_rate_from_qv(qualities);

    return 1 - err_rate;
}

double qv_calibration4(std::vector<int>& qualities, int nPasses)
{

    int predicted_max_read_q = nPasses / 2 + READ_Q_MAX_BASE;
    double err_rate = err_rate_from_qv(qualities);
    double channelq = -10 * log10(err_rate);

    double delta = channelq - predicted_max_read_q;

    delta = delta >= 0 ? delta : 0;

    for (int& qv : qualities) {
        qv -= delta;
        qv = qv >= 0 ? qv : 0;
        qv = qv > BASE_Q_MAX ? BASE_Q_MAX : qv;
    }

    err_rate = err_rate_from_qv(qualities);

    return 1 - err_rate;
}

template <typename TRead>
std::vector<const TRead*> FilterReads(const std::vector<TRead>& reads,
                                      const ConsensusSettings& settings,
                                      SubreadResultCounter* resultCounter)
{
    // This is a count of subreads removed for bing too short, or too long.
    std::vector<const TRead*> results;

    if (reads.empty()) return results;

    std::vector<size_t> lengths;
    size_t longest = 0;

    // get the lengths for all full-length subreads
    for (const auto& read : reads) {
        longest = std::max(longest, read.Seq.length());
        if ((read.Flags & BAM::ADAPTER_BEFORE) && (read.Flags & BAM::ADAPTER_AFTER) &&
            read.ReadAccuracy >= settings.MinReadScore)
            lengths.emplace_back(read.Seq.length());
    }

    // nonexistent median is just the greatest observed length
    const float median = lengths.empty() ? static_cast<float>(longest) : Median(&lengths);
    size_t maxLen = std::min(2 * static_cast<size_t>(median), settings.MaxLength);

    // if it's too short, return nothing
    if (median < static_cast<float>(settings.MinLength)) {
        resultCounter->FilteredBySize += reads.size();
        return results;
    }
    results.reserve(reads.size());

    for (const auto& read : reads) {
        // if the median exists, then this filters stuff,
        //   otherwise it's twice the longest read and is always true

        if (read.SignalToNoise.Minimum() < settings.MinSNR) {
            resultCounter->ZMWBelowMinSNR += 1;
            results.emplace_back(nullptr);
            PBLOG_DEBUG << "Subreads.FilteredByMinSNR " << read.Id << " " << read.Seq;

        } else if (read.ReadAccuracy < settings.MinReadScore) {
            resultCounter->BelowMinQual += 1;
            results.emplace_back(nullptr);
            PBLOG_DEBUG << "Subreads.FilteredByMinReadScore " << read.Id << " " << read.Seq;

        } else if (read.Seq.length() < maxLen) {
            results.emplace_back(&read);
        } else {
            resultCounter->FilteredBySize += 1;
            PBLOG_DEBUG << "Subreads.FilteredBySize " << read.Id << " " << read.Seq;
            results.emplace_back(nullptr);
        }
    }

    // TODO(lhepler): incorporate per-subread quality here
    // End-to-end reads take priority, hence the lexicographical sort;
    //   always take the read with the least deviation from the median.
    //   In the case of no median, longer reads are prioritized.
    const auto lexForm = [median](const TRead* read) {
        const float l = static_cast<float>(read->Seq.length());
        const float v = std::min(l / median, median / l);

        if (read->Flags & BAM::ADAPTER_BEFORE && read->Flags & BAM::ADAPTER_AFTER)
            return std::make_tuple(v, 0.0f);

        return std::make_tuple(0.0f, v);
    };

    std::stable_sort(results.begin(), results.end(),
                     [&lexForm](const TRead* lhs, const TRead* rhs) {
                         if (lhs == nullptr)
                             return false;
                         else if (rhs == nullptr)
                             return true;

                         return lexForm(lhs) > lexForm(rhs);
                     });

    return results;
}

std::tuple<size_t, size_t> ComputeStartEndBasedOnReverseOrNot(size_t read_seq_len,
                                                              const PoaAlignmentSummary& summary)
{
    size_t readStart = summary.ExtentOnRead.Left();
    size_t readEnd = summary.ExtentOnRead.Right();
    if (summary.ReverseComplementedRead) {
        size_t newReadS = read_seq_len - readEnd;
        size_t newReadE = read_seq_len - readStart;
        readStart = newReadS;
        readEnd = newReadE;
    }
    return std::make_tuple(readStart, readEnd);
}

/**
 * attention. 
 * if ReverseComplemented, summary.ExtentOnRead.Left() should store the reversed subreads start point.
 * if not ReverseComplemented, summary.ExtentOnRead.Left() should store the forward subreads start point.
 * 
*/
template <typename TRead>
boost::optional<MappedRead> ExtractMappedRead(const TRead& read, const PoaAlignmentSummary& summary,
                                              const size_t poaLength,
                                              const ConsensusSettings& settings,
                                              SubreadResultCounter* resultCounter,
                                              const std::string& poa)
{
    constexpr size_t kStickyEnds = 7;

    size_t readStart = summary.ExtentOnRead.Left();
    size_t readEnd = summary.ExtentOnRead.Right();
    size_t tplStart = summary.ExtentOnConsensus.Left();
    size_t tplEnd = summary.ExtentOnConsensus.Right();

    // if we're ADAPTER_BEFORE and _AFTER and mapped nearly end-to-end,
    //   just make it end to end (but for each side, respectively)
    if (summary.ReverseComplementedRead) {
        if (read.Flags & BAM::ADAPTER_BEFORE && (poaLength - tplEnd) <= kStickyEnds)
            tplEnd = poaLength;
        if (read.Flags & BAM::ADAPTER_AFTER && tplStart <= kStickyEnds) tplStart = 0;
    } else {
        if (read.Flags & BAM::ADAPTER_BEFORE && tplStart <= kStickyEnds) tplStart = 0;
        if (read.Flags & BAM::ADAPTER_AFTER && (poaLength - tplEnd) <= kStickyEnds)
            tplEnd = poaLength;
    }

    if (summary.ReverseComplementedRead) {
        size_t newReadS = read.Seq.size() - readEnd;
        size_t newReadE = read.Seq.size() - readStart;
        readStart = newReadS;
        readEnd = newReadE;
    }

    if (readStart > readEnd || readEnd - readStart < settings.MinLength) {
        resultCounter->FilteredBySize += 1;
        PBLOG_DEBUG << "Skipping.Subread " << read.Id << ", too short (<" << settings.MinLength
                    << ')';
        return boost::none;
    } else if (readEnd - readStart > settings.MaxLength) {
        resultCounter->FilteredBySize += 1;
        PBLOG_DEBUG << "Skipping.Subread " << read.Id << ", too long (>" << settings.MaxLength
                    << ')';
        return boost::none;
    }

    const SNR& snr = read.SignalToNoise;
    const std::string& chem = read.Chemistry;

    PBLOG_DEBUG << "MappedRead Reverse:" << summary.ReverseComplementedRead << " "
                << poa.substr(tplStart, tplEnd - tplStart) << " "
                << read.Seq.substr(readStart, readEnd - readStart) << " "
                << read.Seq.substr(summary.ExtentOnRead.Left(),
                                   summary.ExtentOnRead.Right() - summary.ExtentOnRead.Left())
                << " Identity:" << summary.AlignmentIdentity << " " << readStart << "," << readEnd;
    ;

    MappedRead mappedRead(
        Read(read.Id, read.Seq.substr(readStart, readEnd - readStart),
             std::vector<uint8_t>(read.IPD.begin() + readStart, read.IPD.begin() + readEnd),
             std::vector<uint8_t>(read.PulseWidth.begin() + readStart,
                                  read.PulseWidth.begin() + readEnd),
             snr, read.Feats, chem,
             std::vector<uint8_t>(read.CR.begin() + readStart, read.CR.begin() + readEnd)),
        summary.ReverseComplementedRead ? StrandType::REVERSE : StrandType::FORWARD, tplStart,
        tplEnd, (tplStart == 0) ? true : false, (tplEnd == poaLength) ? true : false);

    return boost::make_optional(mappedRead);
}

#if 0
template<typename TRead>
bool ReadAccuracyDescending(const std::pair<size_t, const TRead*>& a,
                            const std::pair<size_t, const TRead*>& b)
{
    return a.second->ReadAccuracy > b.second->ReadAccuracy;
}
#endif

}  // namespace

/// \returns a std::pair containing a std::string for the consensus, and a size_t
//           describing the number of adapter-to-adapter reads successfully added
template <typename TRead>
std::pair<std::string, size_t> PoaConsensus(const std::vector<const TRead*>& reads,
                                            std::vector<SparsePoa::ReadKey>* readKeys,
                                            std::vector<PoaAlignmentSummary>* summaries,
                                            const size_t maxPoaCov,
                                            const ConsensusSettings& settings)
{
    SparsePoa poa;
    size_t cov = 0;
    size_t nPasses = 0;

#if 0
    // create a vector of indices into the original reads vector,
    //   sorted by the ReadAccuracy in descending order
    std::vector<std::pair<size_t, const TRead*>> sorted;

    for (size_t i = 0; i < reads.size(); ++i)
        sorted.emplace_back(std::make_pair(i, reads[i]));

    std::sort(sorted.begin(), sorted.end(), ReadAccuracyDescending<TRead>);
#endif

    // initialize readKeys and resize
    readKeys->clear();
    // readKeys->resize(sorted.size());

    for (const auto read : reads) {
        SparsePoa::ReadKey key = (read == nullptr) ? -1 : poa.OrientAndAddRead(read->Seq, settings);

        // SparsePoa::ReadKey key = poa.OrientAndAddRead(read.second->Seq);
        // readKeys->at(read.first) = key;
        readKeys->emplace_back(key);
        if (key >= 0) {
            // TODO: replace with ComputePassesCnt
            if (read->Flags & BAM::ADAPTER_BEFORE && read->Flags & BAM::ADAPTER_AFTER) ++nPasses;
            if ((++cov) >= maxPoaCov) break;
        }
    }

    // at least 50% of the reads should cover
    // TODO(lhepler) revisit this minimum coverage equation
    const size_t minCov = (cov < 5) ? 1 : (cov + 1) / 2 - 1;
    return std::make_pair(poa.FindConsensus(minCov, settings, &(*summaries))->Sequence, nPasses);
}


// pass unique_ptr by reference to satisfy finickyness wrt move semantics in <future>
//   but then take ownership here with a local unique_ptr
template <typename TChunk>
ResultType<ConsensusType> Consensus(std::unique_ptr<std::vector<TChunk>>& chunksRef,
                                    const ConsensusSettings& settings)
{
    using namespace PacBio::Consensus;

    auto chunks(std::move(chunksRef));
    ResultType<ConsensusType> result;

    if (!chunks) return result;
    // We should only be dealing with chunks of size 1
    if (chunks->size() != 1) {
        throw std::runtime_error("CCS chunk was of size != 1");
    }
    const auto& chunk = chunks->at(0);

    try {
        Timer timer;

        // Do read level SNR filtering first
        size_t readsBelowMinSNR = 0;
        for (const auto& read : chunk.Reads) {
            if (read.SignalToNoise.Minimum() < settings.MinSNR) readsBelowMinSNR++;
        }
        // Only if all reads are below the MinSNR cutoff is this a PoorSNR
        if (readsBelowMinSNR == chunk.Reads.size()) {
            result.SubreadCounter.ZMWBelowMinSNR += chunk.Reads.size();
            result.PoorSNR += 1;
            return result;
        }

        auto reads = FilterReads(chunk.Reads, settings, &result.SubreadCounter);

        if (reads.empty() ||  // Check if subread are present
            std::accumulate(reads.begin(), reads.end(), 0, std::plus<bool>()) == 0) {
            result.NoSubreads += 1;
            PBLOG_DEBUG << "Skipping.Channel" << chunk.Id << ", no high quality subreads available";
            return result;
        }

        // If it is not possible to exceed the minPasses requirement, we will bail here before
        //   generating the POA, filling the matrices and performing all the other checks
        size_t possiblePasses = 0;
        size_t activeReads = 0;
        for (size_t i = 0; i < reads.size(); ++i) {
            if (reads[i] != nullptr) {
                activeReads += 1;
                PBLOG_DEBUG << "PassesCntVersion " << settings.PassesCntVersion;
                possiblePasses = ComputePassesCnt(reads[i]->Flags, possiblePasses, settings);
            }
        }

        //
        if (possiblePasses < settings.MinPasses) {
            result.TooFewPasses1 += 1;
            result.SubreadCounter.ZMWNotEnoughSubReads += activeReads;
            PBLOG_DEBUG << "Skipping.Channel" << chunk.Id << ", not enough possible passes ("
                        << possiblePasses << '<' << settings.MinPasses << ')';
            return result;
        }

        std::vector<SparsePoa::ReadKey> readKeys;
        std::vector<PoaAlignmentSummary> summaries;
        std::string poaConsensus;
        size_t nPasses = 0;
        std::tie(poaConsensus, nPasses) =
            PoaConsensus(reads, &readKeys, &summaries, settings.MaxPoaCoverage, settings);

        if (poaConsensus.length() < settings.MinLength) {
            result.TooShort += 1;
            result.SubreadCounter.Other += activeReads;
            PBLOG_DEBUG << "Skipping.Channel" << chunk.Id << ", initial consensus too short (<"
                        << settings.MinLength << ')';
            PBLOG_DEBUG << "Reads.Poa.FilteredByMinLength" << " " << chunk.Id << " "
                        << poaConsensus;

        } else if (poaConsensus.length() > settings.MaxLength) {
            result.TooLong += 1;
            result.SubreadCounter.Other += activeReads;
            // PBLOG_DEBUG << "Skipping.Channel" << chunk.Id << ", initial consensus too long (>"
            //             << settings.MaxLength << ')';
            // PBLOG_DEBUG << "Reads.Poa.FilteredByMaxLength" << " " << chunk.Id << " "
            //             << poaConsensus;

        } else {
            if (settings.NoPolish) {
                const size_t len = poaConsensus.length();
                nPasses = ComputeNPassesForNoPolish(reads, readKeys, summaries,
                                                    settings.MinIdentity, result);

                if (nPasses < settings.MinPasses) {
                    // Reassign all the successful reads to the other category
                    result.SubreadCounter.AssignSuccessToOther();
                    result.TooFewPasses2 += 1;
                } else {
                    // generate dummy QVs, will use 25
                    // TODO(lhepler): should we use different values for delQVs, insQVs, and subQVs?
                    QualityValues qvs{std::vector<int>(len, 3), std::vector<int>(len, 3),
                                    std::vector<int>(len, 3), std::vector<int>(len, 3)};
                    result.Success += 1;
                    result.SubreadCounter.Success += activeReads;

                    result.emplace_back(ConsensusType{
                        PolishResult(), chunk.Id, boost::none, poaConsensus, qvs, nPasses, 0.5, 0,
                        std::vector<double>(1), result.SubreadCounter.ReturnCountsAsArray(),
                        timer.ElapsedMilliseconds(),
                        boost::make_optional(chunk.Reads[0].SignalToNoise), chunk.Barcodes});
                }

            } else {
                const auto mkConsensus = [&](const boost::optional<StrandType> strand) {
                    // give this consensus attempt a name we can refer to
                    std::string chunkName(chunk.Id);
                    if (strand && *strand == StrandType::FORWARD) chunkName += " [fwd]";
                    if (strand && *strand == StrandType::REVERSE) chunkName += " [rev]";

                    try {
                        // setup the arrow integrator
                        IntegratorConfig cfg(settings.MinZScore);
                        Integrator ai(poaConsensus, cfg);
                        const size_t nReads = readKeys.size();
                        size_t nPasses = 0, nDropped = 0;

                        // If this ZMW could possibly pass,  add the reads to the integrator
                        for (size_t i = 0; i < nReads; ++i) {
                            // skip unadded reads
                            if (readKeys[i] < 0) continue;
                            // skip reads that are not sufficiently similar
                            if (summaries[readKeys[i]].AlignmentIdentity < settings.MinIdentity) {

                                const PoaAlignmentSummary& summary = summaries[readKeys[i]];

                                PBLOG_DEBUG << "Skipping.Subread " << reads[i]->Id
                                            << ", poor identity";
                                PBLOG_DEBUG << "Subreads.FilteredByMinIdentity " << chunk.Id
                                            << " Reverse:" << summary.ReverseComplementedRead << " "
                                            << poaConsensus << " " << reads[i]->Seq
                                            << " Identity:" << summary.AlignmentIdentity;

                                size_t readS;
                                size_t readE;
                                std::tie(readS, readE) = ComputeStartEndBasedOnReverseOrNot(
                                    reads[i]->Seq.size(), summary);

                                PBLOG_DEBUG
                                    << "Subreads.FilteredByMinIdentityMapped" << " " << chunk.Id
                                    << " Reverse:" << summary.ReverseComplementedRead << " "
                                    << poaConsensus.substr(summary.ExtentOnConsensus.Left(),
                                                           summary.ExtentOnConsensus.Right() -
                                                               summary.ExtentOnConsensus.Left())
                                    << " " << reads[i]->Seq.substr(readS, readE - readS)
                                    << " Identity:" << summary.AlignmentIdentity << " " << readS
                                    << "," << readE;
                                ;

                                if (settings.ReComputePoaIdentityVersion == "v1") {
                                    result.SubreadCounter.PoorIdentity += 1;
                                    continue;

                                } else if (settings.ReComputePoaIdentityVersion == "v2") {
                                    if ((summaries[readKeys[i]].AlignmentIdentity + 0.1) >=
                                        settings.MinIdentity) {

                                        std::string query =
                                            reads[i]->Seq.substr(readS, readE - readS);
                                        if (summary.ReverseComplementedRead) {
                                            query = ::PacBio::Data::ReverseComplement(query);
                                        }

                                        std::string target = poaConsensus.substr(
                                            summary.ExtentOnConsensus.Left(),
                                            summary.ExtentOnConsensus.Right() -
                                                summary.ExtentOnConsensus.Left());
                                        float new_identity =
                                            Geneus::Align::EdlibGlobalAlign(query, target);
                                        PBLOG_DEBUG << "Subreads.ReComputePoaIdentity" << " "
                                                    << chunk.Id << " identity:" << new_identity
                                                    << " " << target << " " << query;

                                        if (new_identity < settings.MinIdentity) {
                                            result.SubreadCounter.PoorIdentity += 1;
                                            continue;
                                        }
                                        PBLOG_DEBUG << "Subreads.RescuedByReComputePoaIdentity"
                                                    << " " << chunk.Id << " " << poaConsensus << " "
                                                    << reads[i]->Seq;

                                    } else {
                                        result.SubreadCounter.PoorIdentity += 1;
                                        continue;
                                    }

                                } else {
                                    PBLOG_FATAL << "invalid ReComputePoaIdentityVersion, expected: "
                                                   "v1,v2. but got["
                                                << settings.ReComputePoaIdentityVersion << "]";
                                    continue;
                                }
                            }

                            if (auto mr = ExtractMappedRead(*reads[i], summaries[readKeys[i]],
                                                            poaConsensus.length(), settings,
                                                            &result.SubreadCounter, poaConsensus)) {
                                // skip reads not belonging to this strand, if we're --byStrand
                                if (strand && mr->Strand != *strand) continue;
                                auto status = ai.AddRead(*mr);
                                // increment the status count
                                result.SubreadCounter.AddResult(status);
                                if (status == State::VALID) {
                                    nPasses = ComputePassesCnt(reads[i]->Flags, nPasses, settings);

                                } else if (status != State::VALID) {
                                    nDropped += 1;
                                    PBLOG_DEBUG << "Skipping.Subread " << mr->Name << ", "
                                                << status;
                                    PBLOG_DEBUG << "Subreads.FilteredByStatus" << " " << chunk.Id
                                                << " " << poaConsensus << " " << reads[i]->Seq;
                                }
                            }
                        }

                        if (nPasses < settings.MinPasses) {
                            // Reassign all the successful reads to the other category
                            result.SubreadCounter.AssignSuccessToOther();
                            result.TooFewPasses2 += 1;
                            PBLOG_DEBUG << "Skipping.Channel" << chunkName
                                        << ", insufficient number of passes (" << nPasses << '<'
                                        << settings.MinPasses << ')';
                            return;
                        }

                        // this is hairy, but also relatively straightforward, so bear with me:
                        //   if we're not doing strand-specific consensus, the total number of
                        //   available reads is just nReads. If we're doing byStrand though,
                        //   then the number of available reads is those that mapped to this
                        //   strand, plus half of those that didn't map to the POA (we assume).
                        const size_t nAvail =
                            (!strand) ? nReads
                                      : (std::count_if(
                                             readKeys.cbegin(), readKeys.cend(),
                                             [&](const SparsePoa::ReadKey key) {
                                                 return key >= 0 &&
                                                        summaries[key].ReverseComplementedRead ==
                                                            (*strand == StrandType::REVERSE);
                                             }) +
                                         std::count_if(
                                             readKeys.cbegin(), readKeys.cend(),
                                             [](const SparsePoa::ReadKey key) { return key < 0; }) /
                                             2);

                        const double fracDropped = static_cast<double>(nDropped) / nAvail;
                        if (fracDropped > settings.MaxDropFraction) {
                            result.TooManyUnusable += 1;
                            result.SubreadCounter.AssignSuccessToOther();
                            PBLOG_DEBUG << "Skipping.Channel" << chunkName
                                        << ", too high a fraction of unusable subreads ("
                                        << fracDropped << '>' << settings.MaxDropFraction << ')';
                            return;
                        }

                        const double zAvg = ai.AvgZScore();
                        const auto zScores = ai.ZScores();

                        // find consensus!!
                        const PolishResult polishResult = Polish(&ai, PolishConfig());

                        if (!polishResult.hasConverged) {
                            result.NonConvergent += 1;
                            result.SubreadCounter.AssignSuccessToOther();
                            PBLOG_DEBUG << "Skipping.Channel" << chunkName
                                        << ", failed to converge";
                            return;
                        }

                        // compute predicted accuracy
                        QualityValues qvs = ConsensusQVs(ai);
                        double predAcc = qv_calibration4(qvs.Qualities, nPasses);

                        if (predAcc < settings.MinPredictedAccuracy) {
                            result.PoorQuality += 1;
                            result.SubreadCounter.AssignSuccessToOther();
                            PBLOG_DEBUG << "Skipping.Channel" << chunkName
                                        << ", failed to meet minimum predicted accuracy ("
                                        << predAcc << '<' << settings.MinPredictedAccuracy << ')';
                            return;
                        }

                        // return resulting sequence!!
                        result.Success += 1;
                        result.emplace_back(ConsensusType{
                            polishResult, chunk.Id, strand, std::string(ai), std::move(qvs),
                            nPasses, predAcc, zAvg, zScores,
                            result.SubreadCounter.ReturnCountsAsArray(),
                            timer.ElapsedMilliseconds(),
                            boost::make_optional(chunk.Reads[0].SignalToNoise), chunk.Barcodes});
                    } catch (const std::exception& e) {
                        result.ExceptionThrown += 1;
                        PBLOG_ERROR << "Skipping.Channel" << chunkName << ", caught exception: '"
                                    << e.what() << "\'";
                    }
                };

                if (settings.ByStrand) {
                    mkConsensus(StrandType::FORWARD);
                    mkConsensus(StrandType::REVERSE);
                } else {
                    mkConsensus(boost::none);
                }
            }
        }
    } catch (const std::exception& e) {
        result.ExceptionThrown += 1;
        PBLOG_ERROR << "Skipping.Channel" << chunk.Id << ", caught exception: '" << e.what()
                    << "\'";
    } catch (...) {
        // This should NEVER happen. Only here as a guard, if this is ever printed someone
        // goofed
        // up by throwing something that didn't derive from std::exception.
        result.ExceptionThrown += 1;
        PBLOG_ERROR << "Skipping.Channel" << chunk.Id << ", caught unknown exception type";
    }

    return result;
}

}  // namespace CCS
}  // namespace PacBio
