#include <tuple>
#include <string>
#include <vector>
#include "SparsePoa.h"

using namespace Poa;

std::pair<std::string, size_t> PoaConsensus(const std::vector<const TRead *> &reads,
                                            std::vector<SparsePoa::ReadKey> *readKeys,
                                            std::vector<PoaAlignmentSummary> *summaries,
                                            const size_t maxPoaCov)
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

    for (const auto read : reads)
    {
        SparsePoa::ReadKey key = (read == nullptr) ? -1 : poa.OrientAndAddRead(read->Seq);
        // SparsePoa::ReadKey key = poa.OrientAndAddRead(read.second->Seq);
        // readKeys->at(read.first) = key;
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
    return poa.FindConsensus(minCov, &(*summaries))->Sequence;
}