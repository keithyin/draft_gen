// Author: Lance Hepler

#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <string>

using namespace std::literals::string_literals;  // for std::operator ""s

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/internal/BaseEncoding.h>
#include <pacbio/exception/StateError.h>

using PacBio::Data::State;
using PacBio::Exception::StateError;

namespace PacBio {
namespace Consensus {
namespace {

static constexpr const char idx2base[4] = {'A', 'C', 'G', 'T'};

template <typename T>
inline UNANIMITY_CONSTEXPR T clip(const T val, const T (&range)[2])
{
    return std::max(range[0], std::min(val, range[1]));
}

inline UNANIMITY_CONSTEXPR uint8_t EncodeBase(const char base)
{
    // Encode just the base without its pulsewidth into the form:
    //
    //   zzzzzzBB
    //
    // where zzzz are just 0-filled padding bits with the following
    // numeric correspondence:
    //
    //   Aa -> 0
    //   Cc -> 1
    //   Gg -> 2
    //   Tt -> 3
    //
    // This 2-bit format is also known as the NCBI2na format
    // Reference:
    //   https://www.ncbi.nlm.nih.gov/IEB/ToolBox/SDKDOCS/BIOSEQ.HTML

    const uint8_t em = Data::detail::ASCIIToNCBI2naImpl(base);
    if (em > 3U) throw StateError(State::ILLEGAL_BASE, "invalid base in read!");
    return em;
}

inline UNANIMITY_CONSTEXPR uint8_t EncodeBase(const char base, const uint8_t raw_pw)
{
    // Encode the base AND its pulsewidth into the form:
    //
    //   zzzzWWBB
    //
    // where zzzz are just padding bits, WW are the two bits for
    // the 2-bit pulsewidth value and BB are two bits in NCBI2na
    // format.

    if (raw_pw < 1U) throw StateError(State::ILLEGAL_PW, "invalid PulseWidth in read!");
    const uint8_t pw = std::min(2, raw_pw - 1);
    const uint8_t em = (pw << 2) | EncodeBase(base);
    if (em > 11U) throw StateError(State::INVALID, "read encoding error!");
    return em;
}


inline uint8_t EncodeBaseV2(const char base, const uint8_t raw_pw, const std::vector<uint8_t>& boundaries)
{
    // Encode the base AND its pulsewidth into the form:
    //
    //   zzzzWWBB
    //
    // where zzzz are just padding bits, WW are the two bits for
    // the 2-bit pulsewidth value and BB are two bits in NCBI2na
    // format.
    auto it = std::lower_bound(boundaries.begin(), boundaries.end(), static_cast<uint8_t>(raw_pw));
    int insertion_position = 2;
    if (it != boundaries.end()) {
        insertion_position = it - boundaries.begin();
    }
    insertion_position = std::min(2, insertion_position);
    const uint8_t pw = static_cast<uint8_t>(insertion_position);
    const uint8_t em = (pw << 2) | EncodeBase(base);
    if (em > 11U) throw StateError(State::INVALID, "read encoding error!");
    return em;
}

inline uint8_t EncodeBaseV3(const char base, const uint8_t raw_pw, const std::vector<uint8_t>& pw_boundaries, const uint8_t raw_cr, const std::vector<uint8_t>& cr_boundaries)
{
    // Encode the base AND its pulsewidth into the form:
    //
    //   zzCCPPBB
    //
    // where zzzz are just padding bits, WW are the two bits for
    // the 2-bit pulsewidth value and BB are two bits in NCBI2na
    // format.
    auto it = std::lower_bound(pw_boundaries.begin(), pw_boundaries.end(), static_cast<uint8_t>(raw_pw));
    int insertion_position = pw_boundaries.size();
    if (it != pw_boundaries.end()) {
        insertion_position = it - pw_boundaries.begin();
    }
    const uint8_t pw = static_cast<uint8_t>(insertion_position);

    it = std::lower_bound(cr_boundaries.begin(), cr_boundaries.end(), raw_cr);
    insertion_position = cr_boundaries.size();
    if (it != cr_boundaries.end()) {
        insertion_position = it - cr_boundaries.begin();
    }
    const uint8_t cr = static_cast<uint8_t>(insertion_position);


    const uint8_t em = (cr << 4) | (pw << 2) | EncodeBase(base);
    if (em > 63U) throw StateError(State::INVALID, "read encoding error!");
    return em;
}

// context order for A=0, C=1, G=2, T=3:
//   AA, CC, GG, TT, NA, NC, NG, NT
inline UNANIMITY_CONSTEXPR uint8_t EncodeContext8(const NCBI2na prev, const NCBI2na curr)
{
    return ((prev.Data() != curr.Data()) << 2) | curr.Data();
}

// context order for A=0, C=1, G=2, T=3:
//   AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT
inline UNANIMITY_CONSTEXPR uint8_t EncodeContext16(const NCBI2na prev, const NCBI2na curr)
{
    return (prev.Data() << 2) | curr.Data();
}

// first: base
// second: pw
inline UNANIMITY_CONSTEXPR std::pair<char, uint8_t> DecodeEmission(const uint8_t em)
{
    if (em > 11U) throw std::runtime_error("encoded emission value is invalid!");
    const uint8_t NCBI2na = em & 3;
    const uint8_t pw = (em >> 2) + 1;
    if (pw > 3U) throw std::runtime_error("invalid generated PulseWidth!");
    return {Data::detail::NCBI2naToASCIIImpl(NCBI2na), pw};
}

// Emission probabilities overloads
// The idea of these emission overloads is to reduce code duplication
// The way models are currently specified is in a hierarchy, from simplest
// to most parametrized.

// * Simple (only-single base context, match and mismatch)
//   1. P6C4
//   2. Snr
inline UNANIMITY_CONSTEXPR double EmissionTableLookup(const double (&emissionTable)[3][1][2],
                                                      const MoveType move, const uint8_t emission,
                                                      const NCBI2na prev, const NCBI2na curr)
{
    assert(move != MoveType::DELETION);

    return emissionTable[static_cast<uint8_t>(move)][0][curr.Data() != emission];
}

// * Di-nucleotide context, unequal neighbor context
//   3. Marginal
//   4. S_P1C1Beta
inline UNANIMITY_CONSTEXPR double EmissionTableLookup(const double (&emissionTable)[3][8][4],
                                                      const MoveType move, const uint8_t emission,
                                                      const NCBI2na prev, const NCBI2na curr)
{
    assert(move != MoveType::DELETION);

    const auto row = EncodeContext8(prev, curr);
    return emissionTable[static_cast<uint8_t>(move)][row][emission];
}

// * Di-nucleotide, full combinatorial context
//   5. PwSnrA
//   6. PwSnr
//   7. S_P1C1v1
//   8. S_P1C1v2
//   9. S_P2C2v5
inline UNANIMITY_CONSTEXPR double EmissionTableLookup(const double (&emissionTable)[3][16][12],
                                                      const MoveType move, const uint8_t emission,
                                                      const NCBI2na prev, const NCBI2na curr)
{
    assert(move != MoveType::DELETION);

    const auto row = EncodeContext16(prev, curr);
    return emissionTable[static_cast<uint8_t>(move)][row][emission];
}

// S_P2C2v6
inline UNANIMITY_CONSTEXPR double EmissionTableLookup(const double (&emissionTable)[3][16][64],
                                                      const MoveType move, const uint8_t emission,
                                                      const NCBI2na prev, const NCBI2na curr)
{
    assert(move != MoveType::DELETION);

    const auto row = EncodeContext16(prev, curr);
    return emissionTable[static_cast<uint8_t>(move)][row][emission];
}

template <size_t EmissionContextNumber, size_t EmissionOutcomeNumber>
inline UNANIMITY_CONSTEXPR double AbstractEmissionPr(
    const double (&emissionTable)[3][EmissionContextNumber][EmissionOutcomeNumber],
    const MoveType move, const uint8_t emission, const AlleleRep& prev, const AlleleRep& curr)
{
    // constrain template
    static_assert(
        EmissionContextNumber == 1 || EmissionContextNumber == 8 || EmissionContextNumber == 16,
        "The emission context has to be 1 (simple), 8 (di-nuc unequal) or 16 (di-nuc full)!");
    static_assert(
        EmissionOutcomeNumber == 2 || EmissionOutcomeNumber == 4 || EmissionOutcomeNumber == 12,
        "The emission outcome has to be 2 (simple), 4 (di-nuc unequal) or 12 (di-nuc full)!");

    assert(move != MoveType::DELETION);

    // Recall that 0 in NCBI4na indicates a gap
    // which is non-sensical for an emission
    assert(prev.IsValid());
    assert(curr.IsValid());

    using namespace Data::detail;

    if (prev.IsPure() && curr.IsPure()) {
        // pure haploid context
        const NCBI2na prevNCBI2na = prev.GetNCBI2na();
        const NCBI2na currNCBI2na = curr.GetNCBI2na();

        return EmissionTableLookup(emissionTable, move, emission, prevNCBI2na, currNCBI2na);
    } else {
        assert(prev.IsAmbig() || curr.IsAmbig());
        // diploid context
        double result = 0;

        for (uint8_t p = 0; p < 4; ++p) {
            const auto pNCBI2na = NCBI2na::FromRaw(p);

            for (uint8_t c = 0; c < 4; ++c) {
                const auto cNCBI2na = NCBI2na::FromRaw(c);

                if (prev.Contains(pNCBI2na) && curr.Contains(cNCBI2na)) {
                    result +=
                        EmissionTableLookup(emissionTable, move, emission, pNCBI2na, cNCBI2na);
                }
            }
        }

        // normalize by the combinatorial space
        // say we have context 'RR', that is, we
        // have two adjacent loci having 50% A and 50% G
        // The resulting mixture model is
        //
        //   PrEm(emission | tpl = AA) * P(tpl = AA) +
        //   PrEm(emission | tpl = AC) * P(tpl = AC) +
        //   PrEm(emission | tpl = CA) * P(tpl = CA) +
        //   PrEm(emission | tpl = CC) * P(tpl = CC)
        //
        // Given that P(tpl = AA) is just equal to the
        // reciprocal of the cardinality of the
        // combinatorial space, then all terms
        // P(tpl = XY) boil down to 1/|space|, i.e. 1/4
        result /= (prev.NumContainedBases() * curr.NumContainedBases());
        return result;
    }
}

// S_P2C2v6
inline UNANIMITY_CONSTEXPR double EmissionTableLookupV2(const double (&emissionTable)[3][16][64],
                                                      const MoveType move, const uint8_t emission,
                                                      const NCBI2na prev, const NCBI2na curr)
{
    assert(move != MoveType::DELETION);

    const auto row = EncodeContext16(prev, curr);
    return emissionTable[static_cast<uint8_t>(move)][row][emission];
}

template <size_t EmissionContextNumber, size_t EmissionOutcomeNumber>
inline UNANIMITY_CONSTEXPR double AbstractEmissionPrV2(
    const double (&emissionTable)[3][EmissionContextNumber][EmissionOutcomeNumber],
    const MoveType move, const uint8_t emission, const AlleleRep& prev, const AlleleRep& curr)
{
    // constrain template
    static_assert(
        EmissionContextNumber == 1 || EmissionContextNumber == 8 || EmissionContextNumber == 16,
        "The emission context has to be 1 (simple), 8 (di-nuc unequal) or 16 (di-nuc full)!");
    static_assert(
        EmissionOutcomeNumber == 2 || EmissionOutcomeNumber == 4 || EmissionOutcomeNumber == 64,
        "The emission outcome has to be 2 (simple), 4 (di-nuc unequal) or 64 (di-nuc full)!");

    assert(move != MoveType::DELETION);

    // Recall that 0 in NCBI4na indicates a gap
    // which is non-sensical for an emission
    assert(prev.IsValid());
    assert(curr.IsValid());

    using namespace Data::detail;

    if (prev.IsPure() && curr.IsPure()) {
        // pure haploid context
        const NCBI2na prevNCBI2na = prev.GetNCBI2na();
        const NCBI2na currNCBI2na = curr.GetNCBI2na();

        return EmissionTableLookup(emissionTable, move, emission, prevNCBI2na, currNCBI2na);
    } else {
        assert(prev.IsAmbig() || curr.IsAmbig());
        // diploid context
        double result = 0;

        for (uint8_t p = 0; p < 4; ++p) {
            const auto pNCBI2na = NCBI2na::FromRaw(p);

            for (uint8_t c = 0; c < 4; ++c) {
                const auto cNCBI2na = NCBI2na::FromRaw(c);

                if (prev.Contains(pNCBI2na) && curr.Contains(cNCBI2na)) {
                    result +=
                        EmissionTableLookup(emissionTable, move, emission, pNCBI2na, cNCBI2na);
                }
            }
        }

        // normalize by the combinatorial space
        // say we have context 'RR', that is, we
        // have two adjacent loci having 50% A and 50% G
        // The resulting mixture model is
        //
        //   PrEm(emission | tpl = AA) * P(tpl = AA) +
        //   PrEm(emission | tpl = AC) * P(tpl = AC) +
        //   PrEm(emission | tpl = CA) * P(tpl = CA) +
        //   PrEm(emission | tpl = CC) * P(tpl = CC)
        //
        // Given that P(tpl = AA) is just equal to the
        // reciprocal of the cardinality of the
        // combinatorial space, then all terms
        // P(tpl = XY) boil down to 1/|space|, i.e. 1/4
        result /= (prev.NumContainedBases() * curr.NumContainedBases());
        return result;
    }
}

// Generic cache expectation interface
template <typename Callable>
inline UNANIMITY_CONSTEXPR double AbstractExpectedLLForEmission(const MoveType move,
                                                                const AlleleRep& prev,
                                                                const AlleleRep& curr,
                                                                const MomentType moment,
                                                                Callable cacheExpectationFetcher)
{
    // Recall that 0 in NCBI4na indicates a gap
    // which is non-sensical for an emission
    assert(prev.IsValid());
    assert(curr.IsValid());

    using namespace Data::detail;

    if (prev.IsPure() && curr.IsPure()) {
        // pure haploid context
        const NCBI2na prevNCBI2na = prev.GetNCBI2na();
        const NCBI2na currNCBI2na = curr.GetNCBI2na();

        return cacheExpectationFetcher(move, prevNCBI2na, currNCBI2na, moment);
    } else {
        assert(prev.IsAmbig() || curr.IsAmbig());
        // diploid context
        double result = 0;

        for (uint8_t p = 0; p < 4; ++p) {
            const auto pNCBI2na = NCBI2na::FromRaw(p);

            for (uint8_t c = 0; c < 4; ++c) {
                const auto cNCBI2na = NCBI2na::FromRaw(c);

                if (prev.Contains(pNCBI2na) && curr.Contains(cNCBI2na)) {
                    result += cacheExpectationFetcher(move, pNCBI2na, cNCBI2na, moment);
                }
            }
        }

        // normalize by the combinatorial space
        result /= (prev.NumContainedBases() * curr.NumContainedBases());
        return result;
    }
}

// Generic population interface
// rowFetcher is a generic function that takes the previous and current base as
// arguments in NCBI2na encoding. The return value *ret* has to be an array-like
// object of size 4, where values ret[0] + ret[1] + ret[2] + ret[3] sum to 1,
// otherwise the resulting weighted transition probabilities are invalid.
template <typename Callable>
inline std::vector<TemplatePosition> AbstractPopulater(const std::string& tpl, Callable rowFetcher)
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    result.reserve(tpl.size());

    // calculate transition probabilities
    auto prev = AlleleRep::FromASCII(tpl[0], false);
    if (!prev.IsValid())
        throw std::invalid_argument("invalid character ('"s + tpl[0] + "', ordinal "s +
                                    std::to_string(static_cast<int>(tpl[0])) +
                                    ") in template at position 0!"s);

    for (size_t i = 1; i < tpl.size(); ++i) {
        const auto curr = AlleleRep::FromASCII(tpl[i], false);
        if (!curr.IsValid())
            throw std::invalid_argument("invalid character ('"s + tpl[i] + "', ordinal "s +
                                        std::to_string(static_cast<int>(tpl[i])) +
                                        ") in template at position "s + std::to_string(i) + '!');

        // 1. Perform a weighted averaging of
        //    the transition probabilities
        std::array<double, 4> weightedTrans{{0, 0, 0, 0}};
        for (uint8_t p = 0; p < 4; ++p) {
            const auto pNCBI2na = NCBI2na::FromRaw(p);

            for (uint8_t c = 0; c < 4; ++c) {
                const auto cNCBI2na = NCBI2na::FromRaw(c);

                if (prev.Contains(pNCBI2na) && curr.Contains(cNCBI2na)) {
                    const auto& params = rowFetcher(pNCBI2na, cNCBI2na);

                    for (uint8_t i = 0; i < 4; ++i) {
                        weightedTrans[i] += params[i];
                    }
                }
            }
        }

        // 2. Renormalize the weighted transition probabilities
        //    by the combinatorial space
        const uint8_t cominatorialSpace = prev.NumContainedBases() * curr.NumContainedBases();
        for (uint8_t i = 0; i < 4; ++i) {
            weightedTrans[i] /= cominatorialSpace;
        }

        // 3. Finally populate Template vector
        result.emplace_back(TemplatePosition{
            tpl[i - 1],
            weightedTrans[0],  // match
            weightedTrans[1],  // branch
            weightedTrans[2],  // stick
            weightedTrans[3]   // deletion
        });
        prev = curr;
    }

    result.emplace_back(TemplatePosition{tpl.back(), 1.0, 0.0, 0.0, 0.0});

    return result;
}

std::array<char, 2> Ctx2Bases(size_t ctx) {
    const uint8_t bp = ctx & 3;
    char second = idx2base[bp];
    char first = idx2base[ctx >> 2];
    return std::array<char, 2>({first, second});
}

template <typename Iter, typename V>
int FindInsertionPosition(Iter iter_begin, Iter iter_end, V v, int length) {  

    auto it = std::lower_bound(iter_begin, iter_end, v);
    int insertion_position = length;
    if (it != iter_end) {
        insertion_position = it - iter_begin;
    }
    return insertion_position;
}

}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
