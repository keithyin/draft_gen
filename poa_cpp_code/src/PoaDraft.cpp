// Author: Lance Hepler

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

using SparsePoa = PacBio::Poa::SparsePoa;
using PoaAlignmentSummary = PacBio::Poa::PoaAlignmentSummary;
using ConsensusSettings = PacBio::CCS::ConsensusSettings;

extern "C"
{
    struct Subread
    {
        char *seq;
        int flags;
    };

    struct Setting
    {
        float min_identity = 0.82;
        int match_score = 3;
        int mismatch_score = -5;
        int insertion_score = -2;
        int deletion_score = -2;
    };

    struct Result
    {
        char *seq;
        size_t n_passes;
    };

    size_t _ComputeNPassesForNoPolish(const Subread *reads,
                                      size_t num_sbr,
                                      std::vector<SparsePoa::ReadKey> &readKeys,
                                      std::vector<PoaAlignmentSummary> &summaries, float minIdentity)
    {

        const size_t nReads = num_sbr;
        size_t nPasses = 0;
        for (int i = 0; i < nReads; i++)
        {
            if (readKeys[i] < 0)
                continue;
            if (summaries[readKeys[i]].AlignmentIdentity < minIdentity)
            {
                continue;
            }

            int flags = (reads + i)->flags;
            if (!(flags & PacBio::BAM::ADAPTER_BEFORE && flags & PacBio::BAM::ADAPTER_AFTER))
                continue;
            nPasses += 1;
        }
        return nPasses;
    }

    std::string _PoaDraftGenCore(const Subread *reads,
                                 size_t num_sbr,
                                 std::vector<SparsePoa::ReadKey> *readKeys,
                                 std::vector<PoaAlignmentSummary> *summaries,
                                 const size_t maxPoaCov,
                                 const Setting *settings)
    {
        SparsePoa poa;
        size_t cov = 0;

        readKeys->clear();

        ConsensusSettings consensus_settings;
        consensus_settings.MatchScore = settings->match_score;
        consensus_settings.MismatchScore = settings->mismatch_score;
        consensus_settings.InsertionScore = settings->insertion_score;
        consensus_settings.DeletionScore = settings->deletion_score;
        consensus_settings.PoaErrSpanVersion = std::string("v1");
        consensus_settings.PoaIdentityVersion = std::string("v1");

        const size_t nReads = num_sbr;
        size_t nPasses = 0;
        for (int i = 0; i < nReads; i++)
        {
            std::string seq((reads + i)->seq);
            SparsePoa::ReadKey key = poa.OrientAndAddRead(seq, consensus_settings);

            readKeys->emplace_back(key);
            if (key >= 0)
            {
                if ((++cov) >= maxPoaCov)
                    break;
            }
        }

        // at least 50% of the reads should cover
        // TODO(lhepler) revisit this minimum coverage equation
        const size_t minCov = (cov < 5) ? 1 : (cov + 1) / 2 - 1;
        return poa.FindConsensus(minCov, consensus_settings, &(*summaries))->Sequence;
    }

    void FreePoaResult(Result result)
    {
        delete[] result.seq;
    }

    Result PoaDraftGen(const Subread *reads, size_t num_sbr, const Setting *setting)
    {
        std::vector<SparsePoa::ReadKey> readKeys;
        std::vector<PoaAlignmentSummary> summaries;

        std::string consensus_seq = _PoaDraftGenCore(reads, num_sbr, &readKeys, &summaries, 1000000, setting);
        size_t n_passes = _ComputeNPassesForNoPolish(reads, num_sbr, readKeys, summaries, setting->min_identity);

        char *consensus = new char[consensus_seq.size() + 1];
        std::strcpy(consensus, consensus_seq.c_str());
        Result result;
        result.n_passes = n_passes;
        result.seq = consensus;

        return result;
    }
}
