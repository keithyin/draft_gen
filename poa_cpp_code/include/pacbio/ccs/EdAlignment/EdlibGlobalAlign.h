#pragma once
#include "Config.h"
#include "Cigar.h"
#include "edlib.h"
#include <memory>
#include <array>
#include <iosfwd>
#include <string>
#include <string_view>
#include <vector>
#include <utility>
#include <span>
#include "DiffCounts.h"
#include "AlignmentResult.h"
#include "AlignmentTools.h"
using namespace std;
//using namespace SparcSimple;
namespace Geneus {
namespace Align {

///
/// RAII wrapper around edlib's alignment result
///
struct EdlibAlignment
{
    //explicit EdlibAlignment(EdlibAlignResult aln);
    EdlibAlignment(const EdlibAlignment&) = delete;
    EdlibAlignment(EdlibAlignment&&) noexcept = default;
    EdlibAlignment& operator=(const EdlibAlignment&) = delete;
    EdlibAlignment& operator=(EdlibAlignment&&) noexcept = default;
    //~EdlibAlignment() noexcept;
    EdlibAlignment(EdlibAlignResult aln) : Data(std::move(aln)) {}

    ~EdlibAlignment() noexcept { edlibFreeAlignResult(Data); }

    EdlibAlignResult Data;
};

///
/// Align query to target
///
//EdlibAlignment EdlibAlign(const std::string& query, const std::string& target,
//                          const EdlibAlignConfig& config);
//EdlibAlignment EdlibAlign(const char* query, int queryLength, const char* target, int targetLength,
//                          const EdlibAlignConfig& config);





//void EdlibAlignStanlay(vector<StrQueryInfo>* strQueryInfo, const vector<string>& queries, const string& target
//    , const EdlibAlignConfig& config, const Pancake::AlignmentParameters& opt);

float EdlibGlobalAlign(const string& query, const string& target)
{
    //  std::vector<std::unique_ptr<EdlibAlignment>> alignments;
    // vector<Pancake::AlignmentResult> alignments;
    const EdlibAlignConfig config{ -1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0 };

    EdlibAlignResult edlibResult = edlibAlign(query.c_str(), query.size(), target.c_str(), target.size(), config);
        
    if (edlibResult.numLocations == 0)
    {   
         edlibFreeAlignResult(edlibResult);
         return 0.0f;
    }
    
    CCS::AlignmentResult ret;
    Data::Cigar cigar = CCS::EdlibAlignmentToCigar(
            // { edlibResultRes.alignment, static_cast<size_t>(edlibResultRes.alignmentLength) }, ret.diffs);
            edlibResult.alignment, static_cast<size_t>(edlibResult.alignmentLength), ret.diffs);
        
    float similarity = (float)ret.diffs.numEq / (float)(ret.diffs.numEq + ret.diffs.numD + ret.diffs.numI + ret.diffs.numX);


    edlibFreeAlignResult(edlibResult);
    return similarity;
    }

///
/// Convert edlib alignment result to CIGAR
///
//Data::Cigar EdlibAlignmentToCigar(const EdlibAlignment& alignment);
//Data::Cigar EdlibAlignmentToCigar(const unsigned char* alignment, std::int32_t alignmentLength);

}  // namespace Align
}  // namespace PacBio


