// Author: Lance Hepler

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <pacbio/align/AlignConfig.h>
#include <pacbio/data/Sequence.h>
#include <pacbio/denovo/PoaConsensus.h>
#include <pacbio/denovo/PoaGraph.h>
#include <pacbio/ccs/ConsensusSettings.h>

#include <pacbio/ccs/SparseAlignment.h>
#include <pacbio/denovo/SparsePoa.h>

using PacBio::Align::AlignConfig;
using PacBio::Align::AlignMode;
using PacBio::Data::ReverseComplement;
using PacBio::Poa::DefaultPoaConfig;
using PacBio::Poa::PoaConsensus;
using PacBio::Poa::PoaGraph;
using PacBio::Poa::detail::SdpAnchorVector;

namespace PacBio
{
    namespace Poa
    {

        using Vertex = PoaGraph::Vertex;

        void ConsensusSpan(
            const std::vector<Vertex> &readPath,
            std::map<Vertex, size_t> &cssPosition,
            size_t *cssS,
            size_t *cssE,
            size_t *readS,
            size_t *readE,
            size_t *nErr)
        {
            bool foundStart = false;

            for (size_t readPos = 0; readPos < readPath.size(); readPos++)
            {
                Vertex v = readPath[readPos];
                if (cssPosition.find(v) != cssPosition.end())
                {
                    if (!foundStart)
                    {
                        *cssS = cssPosition[v];
                        *readS = readPos;
                        foundStart = true;
                    }

                    *cssE = cssPosition[v] + 1;
                    *readE = readPos + 1;
                }
                else
                {
                    *nErr += 1;
                }
            }
        }

        /**
         * nErr: only care the error in the intersected span
         */
        void ConsensusSpanV2(
            const std::vector<Vertex> &readPath,
            std::map<Vertex, size_t> &cssPosition,
            size_t *cssS,
            size_t *cssE,
            size_t *readS,
            size_t *readE,
            size_t *nErr)
        {
            bool foundStart = false;
            size_t cumErr = 0;
            for (size_t readPos = 0; readPos < readPath.size(); readPos++)
            {
                Vertex v = readPath[readPos];
                if (cssPosition.find(v) != cssPosition.end())
                {
                    if (!foundStart)
                    {
                        *cssS = cssPosition[v];
                        *readS = readPos;
                        foundStart = true;
                    }

                    *cssE = cssPosition[v] + 1;
                    *readE = readPos + 1;

                    *nErr += cumErr;
                    cumErr = 0;
                }
                else
                {
                    if (foundStart)
                    {
                        cumErr += 1;
                    }
                }
            }
        }

        SdpAnchorVector SdpRangeFinder::FindAnchors(const std::string &consensusSequence,
                                                    const std::string &readSequence) const
        {
            return CCS::SparseAlign(6, consensusSequence, readSequence);
        }

        SparsePoa::SparsePoa()
            : graph_(new PoaGraph()), readPaths_(), reverseComplemented_(), rangeFinder_(new SdpRangeFinder())
        {
        }

        SparsePoa::~SparsePoa()
        {
            delete graph_;
            delete rangeFinder_;
        }

        SparsePoa::ReadKey SparsePoa::AddRead(const std::string &readSequence,
                                              const PacBio::CCS::ConsensusSettings &settings,
                                              const PoaAlignmentOptions & /* alnOptions */,
                                              float minScoreToAdd)
        {
            AlignConfig config = DefaultPoaConfig(AlignMode::LOCAL);
            config.Params.Match = settings.MatchScore;
            config.Params.Mismatch = settings.MismatchScore;
            config.Params.Insert = settings.InsertionScore;
            config.Params.Delete = settings.DeletionScore;

            Path outputPath;
            ReadKey key;
            // std::cout << "SparsePoa::AddRead" << std::endl;
            if (graph_->NumReads() == 0)
            {
                graph_->AddFirstRead(readSequence, &outputPath);
                readPaths_.push_back(outputPath);
                reverseComplemented_.push_back(false);
                key = graph_->NumReads() - 1;
            }
            else
            {
                auto c1 = graph_->TryAddRead(readSequence, config, rangeFinder_);

                if (c1->Score() >= minScoreToAdd)
                {
                    graph_->CommitAdd(c1, &outputPath);
                    readPaths_.push_back(outputPath);
                    reverseComplemented_.push_back(false);
                    key = graph_->NumReads() - 1;
                }

                else
                {
                    key = -1;
                }

                delete c1;
            }
            return key;
        }

        SparsePoa::ReadKey SparsePoa::OrientAndAddRead(const std::string &readSequence,
                                                       const PacBio::CCS::ConsensusSettings &settings,
                                                       const PoaAlignmentOptions & /* alnOptions */,
                                                       float minScoreToAdd)
        {
            AlignConfig config = DefaultPoaConfig(AlignMode::LOCAL);
            config.Params.Match = settings.MatchScore;
            config.Params.Mismatch = settings.MismatchScore;
            config.Params.Insert = settings.InsertionScore;
            config.Params.Delete = settings.DeletionScore;

            Path outputPath;
            ReadKey key;

            if (graph_->NumReads() == 0)
            {
                graph_->AddFirstRead(readSequence, &outputPath);
                readPaths_.push_back(outputPath);
                reverseComplemented_.push_back(false);
                key = graph_->NumReads() - 1;
            }
            else
            {
                auto c1 = graph_->TryAddRead(readSequence, config, rangeFinder_);
                auto c2 = graph_->TryAddRead(ReverseComplement(readSequence), config, rangeFinder_);

                if (c1->Score() >= c2->Score() && c1->Score() >= minScoreToAdd)
                {
                    graph_->CommitAdd(c1, &outputPath);
                    readPaths_.push_back(outputPath);
                    reverseComplemented_.push_back(false);
                    key = graph_->NumReads() - 1;
                }
                else if (c2->Score() >= c1->Score() && c2->Score() >= minScoreToAdd)
                {
                    graph_->CommitAdd(c2, &outputPath);
                    readPaths_.push_back(outputPath);
                    reverseComplemented_.push_back(true);
                    key = graph_->NumReads() - 1;
                }
                else
                {
                    key = -1;
                }

                delete c1;
                delete c2;
            }
            return key;
        }

        std::shared_ptr<const PoaConsensus> SparsePoa::FindConsensus(
            int minCoverage, const PacBio::CCS::ConsensusSettings &settings, std::vector<PoaAlignmentSummary> *summaries) const
        {
            AlignConfig config = DefaultPoaConfig(AlignMode::LOCAL);
            std::shared_ptr<const PoaConsensus> pc(graph_->FindConsensus(config, minCoverage));
            std::string css = pc->Sequence;

            if (summaries != nullptr)
            {
                summaries->clear();

                // digest the consensus path consensus into map(vtx, pos)
                // the fold over the readPaths
                std::map<Vertex, size_t> cssPosition;

                int i = 0;
                for (Vertex v : pc->Path)
                {
                    cssPosition[v] = i;
                    i++;
                }

                for (size_t readId = 0; readId < graph_->NumReads(); readId++)
                {
                    size_t readS = 0, readE = 0;
                    size_t cssS = 0, cssE = 0;
                    bool foundStart = false;
                    size_t nErr = 0;

                    const std::vector<Vertex> &readPath = readPaths_[readId];
                    if (settings.PoaErrSpanVersion == std::string("v1"))
                    {
                        ConsensusSpan(readPath, cssPosition, &cssS, &cssE, &readS, &readE, &nErr);
                    }
                    else if (settings.PoaErrSpanVersion == std::string("v2"))
                    {
                        ConsensusSpanV2(readPath, cssPosition, &cssS, &cssE, &readS, &readE, &nErr);
                    }
                    else
                    {
                        exit(EXIT_FAILURE);
                    }

                    Interval readExtent(readS, readE);
                    Interval cssExtent(cssS, cssE);

                    PoaAlignmentSummary summary;
                    summary.ReverseComplementedRead = reverseComplemented_[readId];
                    summary.ExtentOnRead = readExtent;
                    summary.ExtentOnConsensus = cssExtent;
                    summary.AlignmentIdentity = 0.0f;
                    if (settings.PoaIdentityVersion == std::string("v1"))
                    {
                        summary.AlignmentIdentity = std::max(0.0f, 1.0f - 1.0f * nErr / cssPosition.size());
                    }
                    else if (settings.PoaIdentityVersion == std::string("v2"))
                    {
                        summary.AlignmentIdentity = std::max(0.0f, 1.0f - 1.0f * nErr / std::max((readE - readS), (size_t)1));
                    }
                    else
                    {
                        exit(EXIT_FAILURE);
                    }

                    (*summaries).push_back(summary);
                }
            }

            return pc;
        }

        std::string SparsePoa::ToGraphViz(int flags, const PoaConsensus *pc) const
        {
            return graph_->ToGraphViz(flags, pc);
        }

        void SparsePoa::WriteGraphVizFile(const std::string &filename, int flags,
                                          const PoaConsensus *pc) const
        {
            graph_->WriteGraphVizFile(filename, flags, pc);
        }

        void SparsePoa::WriteGraphCsvFile(const std::string &filename) const
        {
            graph_->WriteGraphCsvFile(filename);
        }

        void SparsePoa::PruneGraph(const int minCoverage) { graph_->PruneGraph(minCoverage); }

        void SparsePoa::repCheck()
        {
            assert(graph_->NumReads() == readPaths_.size());
            assert(graph_->NumReads() == reverseComplemented_.size());
        }

    } // namespace Poa
} // namespace PacBio
