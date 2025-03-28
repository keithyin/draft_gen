#pragma once

#include "Config.h"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#if __cplusplus >= 201703L
#include <string_view>
#endif

#include <cctype>
#include <cstdint>

namespace Geneus {
namespace Utility {

inline char Complement(const char base)
{
    constexpr std::array<char, 256> LOOKUP_TABLE{
        {/*   0 -   7: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*   8 -  15: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  16 -  23: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  24 -  31: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  32 -  39: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  40 -  47: */ 0,   0,   '*', 0,   0,   '-', 0,   0,
         /*  48 -  55: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  56 -  63: */ 0,   0,   0,   0,   0,   0,   0,   0,

         /*  64 -  71: */ 0,   'T', 'V', 'G', 'H', 0,   0,   'C',
         /*  72 -  79: */ 'D', 0,   0,   'M', 0,   'K', 'N', 0,
         /*  80 -  87: */ 0,   0,   'Y', 'S', 'A', 'A', 'B', 'W',
         /*  88 -  95: */ 0,   'R', 0,   0,   0,   0,   0,   0,

         /*  96 - 103: */ 0,   'T', 'V', 'G', 'H', 0,   0,   'C',
         /* 104 - 111: */ 'D', 0,   0,   'M', 0,   'K', 'N', 0,
         /* 112 - 119: */ 0,   0,   'Y', 'S', 'A', 'A', 'B', 'W',
         /* 120 - 127: */ 0,   'R', 0,   0,   0,   0,   0,   0,

         /* 128 - 135: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 136 - 143: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 144 - 151: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 152 - 159: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 160 - 167: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 168 - 175: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 176 - 183: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 184 - 191: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 192 - 199: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 200 - 207: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 208 - 215: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 216 - 223: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 224 - 231: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 232 - 239: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 240 - 247: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 248 - 255: */ 0,   0,   0,   0,   0,   0,   0,   0}};

    const char result = LOOKUP_TABLE[static_cast<unsigned char>(base)];

    if (result == 0) {
        throw std::invalid_argument(base + std::string{" is an invalid base!"});
    }

    return result;
}

template <typename T>
void Reverse(T& input)
{
    std::reverse(input.begin(), input.end());
}

template <typename T>
T MaybeReverse(T&& input, bool reverse)
{
    if (reverse) {
        std::reverse(input.begin(), input.end());
    }
    return input;
}

template <typename T>
T Reversed(const T& input)
{
    T result = input;
    Reverse(result);
    return result;
}

inline void ReverseComplement(std::string& seq)
{
    std::transform(seq.begin(), seq.end(), seq.begin(), Complement);
    Reverse(seq);
}

inline std::string MaybeReverseComplement(std::string&& seq, bool reverse)
{
    if (reverse) {
        ReverseComplement(seq);
    }
    return std::move(seq);
}

/// Reverse complement a DNA sequence case-sensitive
inline void ReverseComplementCaseSens(std::string& seq)
{
    const std::string original = seq;
    constexpr std::array<char, 128> RC_TABLE = {
        4,  4, 4, 4, 4, 4, 4,  4,  4, 4, 4, 4, 4,  4,  4, 4,  4,  4, 4, 4,   4, 4,   4, 4, 4, 4,
        4,  4, 4, 4, 4, 4, 32, 4,  4, 4, 4, 4, 4,  4,  4, 4,  42, 4, 4, 45,  4, 4,   4, 4, 4, 4,
        4,  4, 4, 4, 4, 4, 4,  4,  4, 4, 4, 4, 4,  84, 4, 71, 4,  4, 4, 67,  4, 4,   4, 4, 4, 4,
        78, 4, 4, 4, 4, 4, 65, 65, 4, 4, 4, 4, 4,  4,  4, 4,  4,  4, 4, 116, 4, 103, 4, 4, 4, 99,
        4,  4, 4, 4, 4, 4, 4,  4,  4, 4, 4, 4, 97, 97, 4, 4,  4,  4, 4, 4,   4, 4,   4, 4};
    std::string reverseCompl(original.length(), 'N');
    for (std::uint32_t i = 0; i < original.length(); ++i) {
        reverseCompl[original.length() - i - 1] = RC_TABLE[static_cast<unsigned char>(original[i])];
    }
    seq = reverseCompl;
}

inline std::string MaybeReverseComplementCaseSens(std::string&& seq, bool reverse)
{
    if (reverse) {
        ReverseComplementCaseSens(seq);
    }
    return std::move(seq);
}

inline std::string ReverseComplemented(const std::string& input)
{
    std::string result = input;
    ReverseComplement(result);
    return result;
}

#if __cplusplus >= 201703L
inline std::string_view ReverseComplement(const std::string_view input, char* output)
{
    const std::size_t strLen = input.length();
    for (std::size_t i = 0; i < strLen; ++i) {
        output[i] = Complement(input[strLen - 1 - i]);
    }
    return {output, strLen};
}
#endif

}  // namespace Utility
}  // namespace PacBio


