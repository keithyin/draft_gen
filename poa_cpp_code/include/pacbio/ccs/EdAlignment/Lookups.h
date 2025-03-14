// Author: Ivan Sovic

#ifndef PANCAKE_LOOKUPS_HPP
#define PANCAKE_LOOKUPS_HPP

#include <array>
#include <cstdint>

namespace Geneus {
namespace CCS {

// clang-format off
constexpr std::array<int8_t, 256> BASE_TO_TWO_BIT = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0 - 15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 16 - 31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 32 - 47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 48 - 63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 64 - 79 (A, C, G)
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 80 - 95 (T)
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 96 - 111
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 112 - 127
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 128 - 143
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 144 - 159
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 160 - 176
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 176 - 191
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 192 - 208
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 208 - 223
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 224 - 239
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4   // 240 - 256
};

constexpr std::array<int8_t, 256> BaseToTwobitComplement = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0 - 15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 16 - 31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 32 - 47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 48 - 63
    4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4,  // 64 - 79 (A, C, G)
    4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 80 - 95 (T)
    4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4,  // 96 - 111
    4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 112 - 127
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 128 - 143
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 144 - 159
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 160 - 176
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 176 - 191
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 192 - 208
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 208 - 223
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 224 - 239
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4   // 240 - 256
};

constexpr std::array<int8_t, 256> TwobitToBase = {
    65, 67, 71, 84, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 0 - 15 (A, C, G, T, N)
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 16 - 31
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 32 - 47
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 48 - 63
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 64 - 79
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 80 - 95
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 96 - 111
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 112 - 127
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 128 - 143
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 144 - 159
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 160 - 176
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 176 - 191
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 192 - 208
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 208 - 223
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 224 - 239
    0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0   // 240 - 256
};

constexpr std::array<int8_t, 256> BaseToBaseComplement = {
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 0 - 15
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 16 - 31
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 32 - 47
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 48 - 63
    78, 84, 78, 71, 78, 78, 78, 67, 78, 78, 78, 78, 78, 78, 78, 78,  // 64 - 79 (A = 65, C = 67, G = 71)
    78, 78, 78, 78, 65, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 80 - 95 (T = 84)
    78, 116, 78, 103, 78, 78, 78, 99, 78, 78, 78, 78, 78, 78, 78, 78,  // 96 - 111 (a = 97, c = 99, g = 103)
    78, 78, 78, 78, 97, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 112 - 127 (t = 116)
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 128 - 143
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 144 - 159
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 160 - 176
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 176 - 191
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 192 - 208
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 208 - 223
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,  // 224 - 239
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78   // 240 - 256
};

constexpr std::array<int8_t, 256> TwobitToTwobitComplement = {
    3, 2, 1, 0, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 0 - 15 (A, C, G, T, N)
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 16 - 31
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 32 - 47
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 48 - 63
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 64 - 79
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 80 - 95
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 96 - 111
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 112 - 127
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 128 - 143
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 144 - 159
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 160 - 176
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 176 - 191
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 192 - 208
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 208 - 223
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 224 - 239
    0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0   // 240 - 256
};

constexpr std::array<int8_t, 256> IsNucleotide = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 0 - 15
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 16 - 31
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 32 - 47
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 48 - 63
    0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  // 64 - 79 (A, C, G)
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 80 - 95 (T)
    0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  // 96 - 111
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 112 - 127
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 128 - 143
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 144 - 159
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 160 - 176
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 176 - 191
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 192 - 208
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 208 - 223
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 224 - 239
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0   // 240 - 256
};

const char ByteToBases[256][5] = {
    {'A', 'A', 'A', 'A', '\0'}, {'A', 'A', 'A', 'C', '\0'}, {'A', 'A', 'A', 'G', '\0'},
    {'A', 'A', 'A', 'T', '\0'}, {'A', 'A', 'C', 'A', '\0'}, {'A', 'A', 'C', 'C', '\0'},
    {'A', 'A', 'C', 'G', '\0'}, {'A', 'A', 'C', 'T', '\0'}, {'A', 'A', 'G', 'A', '\0'},
    {'A', 'A', 'G', 'C', '\0'}, {'A', 'A', 'G', 'G', '\0'}, {'A', 'A', 'G', 'T', '\0'},
    {'A', 'A', 'T', 'A', '\0'}, {'A', 'A', 'T', 'C', '\0'}, {'A', 'A', 'T', 'G', '\0'},
    {'A', 'A', 'T', 'T', '\0'}, {'A', 'C', 'A', 'A', '\0'}, {'A', 'C', 'A', 'C', '\0'},
    {'A', 'C', 'A', 'G', '\0'}, {'A', 'C', 'A', 'T', '\0'}, {'A', 'C', 'C', 'A', '\0'},
    {'A', 'C', 'C', 'C', '\0'}, {'A', 'C', 'C', 'G', '\0'}, {'A', 'C', 'C', 'T', '\0'},
    {'A', 'C', 'G', 'A', '\0'}, {'A', 'C', 'G', 'C', '\0'}, {'A', 'C', 'G', 'G', '\0'},
    {'A', 'C', 'G', 'T', '\0'}, {'A', 'C', 'T', 'A', '\0'}, {'A', 'C', 'T', 'C', '\0'},
    {'A', 'C', 'T', 'G', '\0'}, {'A', 'C', 'T', 'T', '\0'}, {'A', 'G', 'A', 'A', '\0'},
    {'A', 'G', 'A', 'C', '\0'}, {'A', 'G', 'A', 'G', '\0'}, {'A', 'G', 'A', 'T', '\0'},
    {'A', 'G', 'C', 'A', '\0'}, {'A', 'G', 'C', 'C', '\0'}, {'A', 'G', 'C', 'G', '\0'},
    {'A', 'G', 'C', 'T', '\0'}, {'A', 'G', 'G', 'A', '\0'}, {'A', 'G', 'G', 'C', '\0'},
    {'A', 'G', 'G', 'G', '\0'}, {'A', 'G', 'G', 'T', '\0'}, {'A', 'G', 'T', 'A', '\0'},
    {'A', 'G', 'T', 'C', '\0'}, {'A', 'G', 'T', 'G', '\0'}, {'A', 'G', 'T', 'T', '\0'},
    {'A', 'T', 'A', 'A', '\0'}, {'A', 'T', 'A', 'C', '\0'}, {'A', 'T', 'A', 'G', '\0'},
    {'A', 'T', 'A', 'T', '\0'}, {'A', 'T', 'C', 'A', '\0'}, {'A', 'T', 'C', 'C', '\0'},
    {'A', 'T', 'C', 'G', '\0'}, {'A', 'T', 'C', 'T', '\0'}, {'A', 'T', 'G', 'A', '\0'},
    {'A', 'T', 'G', 'C', '\0'}, {'A', 'T', 'G', 'G', '\0'}, {'A', 'T', 'G', 'T', '\0'},
    {'A', 'T', 'T', 'A', '\0'}, {'A', 'T', 'T', 'C', '\0'}, {'A', 'T', 'T', 'G', '\0'},
    {'A', 'T', 'T', 'T', '\0'}, {'C', 'A', 'A', 'A', '\0'}, {'C', 'A', 'A', 'C', '\0'},
    {'C', 'A', 'A', 'G', '\0'}, {'C', 'A', 'A', 'T', '\0'}, {'C', 'A', 'C', 'A', '\0'},
    {'C', 'A', 'C', 'C', '\0'}, {'C', 'A', 'C', 'G', '\0'}, {'C', 'A', 'C', 'T', '\0'},
    {'C', 'A', 'G', 'A', '\0'}, {'C', 'A', 'G', 'C', '\0'}, {'C', 'A', 'G', 'G', '\0'},
    {'C', 'A', 'G', 'T', '\0'}, {'C', 'A', 'T', 'A', '\0'}, {'C', 'A', 'T', 'C', '\0'},
    {'C', 'A', 'T', 'G', '\0'}, {'C', 'A', 'T', 'T', '\0'}, {'C', 'C', 'A', 'A', '\0'},
    {'C', 'C', 'A', 'C', '\0'}, {'C', 'C', 'A', 'G', '\0'}, {'C', 'C', 'A', 'T', '\0'},
    {'C', 'C', 'C', 'A', '\0'}, {'C', 'C', 'C', 'C', '\0'}, {'C', 'C', 'C', 'G', '\0'},
    {'C', 'C', 'C', 'T', '\0'}, {'C', 'C', 'G', 'A', '\0'}, {'C', 'C', 'G', 'C', '\0'},
    {'C', 'C', 'G', 'G', '\0'}, {'C', 'C', 'G', 'T', '\0'}, {'C', 'C', 'T', 'A', '\0'},
    {'C', 'C', 'T', 'C', '\0'}, {'C', 'C', 'T', 'G', '\0'}, {'C', 'C', 'T', 'T', '\0'},
    {'C', 'G', 'A', 'A', '\0'}, {'C', 'G', 'A', 'C', '\0'}, {'C', 'G', 'A', 'G', '\0'},
    {'C', 'G', 'A', 'T', '\0'}, {'C', 'G', 'C', 'A', '\0'}, {'C', 'G', 'C', 'C', '\0'},
    {'C', 'G', 'C', 'G', '\0'}, {'C', 'G', 'C', 'T', '\0'}, {'C', 'G', 'G', 'A', '\0'},
    {'C', 'G', 'G', 'C', '\0'}, {'C', 'G', 'G', 'G', '\0'}, {'C', 'G', 'G', 'T', '\0'},
    {'C', 'G', 'T', 'A', '\0'}, {'C', 'G', 'T', 'C', '\0'}, {'C', 'G', 'T', 'G', '\0'},
    {'C', 'G', 'T', 'T', '\0'}, {'C', 'T', 'A', 'A', '\0'}, {'C', 'T', 'A', 'C', '\0'},
    {'C', 'T', 'A', 'G', '\0'}, {'C', 'T', 'A', 'T', '\0'}, {'C', 'T', 'C', 'A', '\0'},
    {'C', 'T', 'C', 'C', '\0'}, {'C', 'T', 'C', 'G', '\0'}, {'C', 'T', 'C', 'T', '\0'},
    {'C', 'T', 'G', 'A', '\0'}, {'C', 'T', 'G', 'C', '\0'}, {'C', 'T', 'G', 'G', '\0'},
    {'C', 'T', 'G', 'T', '\0'}, {'C', 'T', 'T', 'A', '\0'}, {'C', 'T', 'T', 'C', '\0'},
    {'C', 'T', 'T', 'G', '\0'}, {'C', 'T', 'T', 'T', '\0'}, {'G', 'A', 'A', 'A', '\0'},
    {'G', 'A', 'A', 'C', '\0'}, {'G', 'A', 'A', 'G', '\0'}, {'G', 'A', 'A', 'T', '\0'},
    {'G', 'A', 'C', 'A', '\0'}, {'G', 'A', 'C', 'C', '\0'}, {'G', 'A', 'C', 'G', '\0'},
    {'G', 'A', 'C', 'T', '\0'}, {'G', 'A', 'G', 'A', '\0'}, {'G', 'A', 'G', 'C', '\0'},
    {'G', 'A', 'G', 'G', '\0'}, {'G', 'A', 'G', 'T', '\0'}, {'G', 'A', 'T', 'A', '\0'},
    {'G', 'A', 'T', 'C', '\0'}, {'G', 'A', 'T', 'G', '\0'}, {'G', 'A', 'T', 'T', '\0'},
    {'G', 'C', 'A', 'A', '\0'}, {'G', 'C', 'A', 'C', '\0'}, {'G', 'C', 'A', 'G', '\0'},
    {'G', 'C', 'A', 'T', '\0'}, {'G', 'C', 'C', 'A', '\0'}, {'G', 'C', 'C', 'C', '\0'},
    {'G', 'C', 'C', 'G', '\0'}, {'G', 'C', 'C', 'T', '\0'}, {'G', 'C', 'G', 'A', '\0'},
    {'G', 'C', 'G', 'C', '\0'}, {'G', 'C', 'G', 'G', '\0'}, {'G', 'C', 'G', 'T', '\0'},
    {'G', 'C', 'T', 'A', '\0'}, {'G', 'C', 'T', 'C', '\0'}, {'G', 'C', 'T', 'G', '\0'},
    {'G', 'C', 'T', 'T', '\0'}, {'G', 'G', 'A', 'A', '\0'}, {'G', 'G', 'A', 'C', '\0'},
    {'G', 'G', 'A', 'G', '\0'}, {'G', 'G', 'A', 'T', '\0'}, {'G', 'G', 'C', 'A', '\0'},
    {'G', 'G', 'C', 'C', '\0'}, {'G', 'G', 'C', 'G', '\0'}, {'G', 'G', 'C', 'T', '\0'},
    {'G', 'G', 'G', 'A', '\0'}, {'G', 'G', 'G', 'C', '\0'}, {'G', 'G', 'G', 'G', '\0'},
    {'G', 'G', 'G', 'T', '\0'}, {'G', 'G', 'T', 'A', '\0'}, {'G', 'G', 'T', 'C', '\0'},
    {'G', 'G', 'T', 'G', '\0'}, {'G', 'G', 'T', 'T', '\0'}, {'G', 'T', 'A', 'A', '\0'},
    {'G', 'T', 'A', 'C', '\0'}, {'G', 'T', 'A', 'G', '\0'}, {'G', 'T', 'A', 'T', '\0'},
    {'G', 'T', 'C', 'A', '\0'}, {'G', 'T', 'C', 'C', '\0'}, {'G', 'T', 'C', 'G', '\0'},
    {'G', 'T', 'C', 'T', '\0'}, {'G', 'T', 'G', 'A', '\0'}, {'G', 'T', 'G', 'C', '\0'},
    {'G', 'T', 'G', 'G', '\0'}, {'G', 'T', 'G', 'T', '\0'}, {'G', 'T', 'T', 'A', '\0'},
    {'G', 'T', 'T', 'C', '\0'}, {'G', 'T', 'T', 'G', '\0'}, {'G', 'T', 'T', 'T', '\0'},
    {'T', 'A', 'A', 'A', '\0'}, {'T', 'A', 'A', 'C', '\0'}, {'T', 'A', 'A', 'G', '\0'},
    {'T', 'A', 'A', 'T', '\0'}, {'T', 'A', 'C', 'A', '\0'}, {'T', 'A', 'C', 'C', '\0'},
    {'T', 'A', 'C', 'G', '\0'}, {'T', 'A', 'C', 'T', '\0'}, {'T', 'A', 'G', 'A', '\0'},
    {'T', 'A', 'G', 'C', '\0'}, {'T', 'A', 'G', 'G', '\0'}, {'T', 'A', 'G', 'T', '\0'},
    {'T', 'A', 'T', 'A', '\0'}, {'T', 'A', 'T', 'C', '\0'}, {'T', 'A', 'T', 'G', '\0'},
    {'T', 'A', 'T', 'T', '\0'}, {'T', 'C', 'A', 'A', '\0'}, {'T', 'C', 'A', 'C', '\0'},
    {'T', 'C', 'A', 'G', '\0'}, {'T', 'C', 'A', 'T', '\0'}, {'T', 'C', 'C', 'A', '\0'},
    {'T', 'C', 'C', 'C', '\0'}, {'T', 'C', 'C', 'G', '\0'}, {'T', 'C', 'C', 'T', '\0'},
    {'T', 'C', 'G', 'A', '\0'}, {'T', 'C', 'G', 'C', '\0'}, {'T', 'C', 'G', 'G', '\0'},
    {'T', 'C', 'G', 'T', '\0'}, {'T', 'C', 'T', 'A', '\0'}, {'T', 'C', 'T', 'C', '\0'},
    {'T', 'C', 'T', 'G', '\0'}, {'T', 'C', 'T', 'T', '\0'}, {'T', 'G', 'A', 'A', '\0'},
    {'T', 'G', 'A', 'C', '\0'}, {'T', 'G', 'A', 'G', '\0'}, {'T', 'G', 'A', 'T', '\0'},
    {'T', 'G', 'C', 'A', '\0'}, {'T', 'G', 'C', 'C', '\0'}, {'T', 'G', 'C', 'G', '\0'},
    {'T', 'G', 'C', 'T', '\0'}, {'T', 'G', 'G', 'A', '\0'}, {'T', 'G', 'G', 'C', '\0'},
    {'T', 'G', 'G', 'G', '\0'}, {'T', 'G', 'G', 'T', '\0'}, {'T', 'G', 'T', 'A', '\0'},
    {'T', 'G', 'T', 'C', '\0'}, {'T', 'G', 'T', 'G', '\0'}, {'T', 'G', 'T', 'T', '\0'},
    {'T', 'T', 'A', 'A', '\0'}, {'T', 'T', 'A', 'C', '\0'}, {'T', 'T', 'A', 'G', '\0'},
    {'T', 'T', 'A', 'T', '\0'}, {'T', 'T', 'C', 'A', '\0'}, {'T', 'T', 'C', 'C', '\0'},
    {'T', 'T', 'C', 'G', '\0'}, {'T', 'T', 'C', 'T', '\0'}, {'T', 'T', 'G', 'A', '\0'},
    {'T', 'T', 'G', 'C', '\0'}, {'T', 'T', 'G', 'G', '\0'}, {'T', 'T', 'G', 'T', '\0'},
    {'T', 'T', 'T', 'A', '\0'}, {'T', 'T', 'T', 'C', '\0'}, {'T', 'T', 'T', 'G', '\0'},
    {'T', 'T', 'T', 'T'}};

constexpr std::array<int32_t, 256> CigarCharToNum = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0 - 16
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 16 - 32
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 32 - 48
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, // 48 - 64, "=": 7
    0, 0, 0, 0, 2, 0, 0, 0, 5, 1, 0, 0, 0, 0, 3, 0, // 64 - 80, "D": 2, "H": 5, "I": 1, "M": 0, "N": 3
    6, 0, 0, 4, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, // 80 - 96, "P": 6, "S": 4, "X": 8
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 96 - 112
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 112 - 128
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 128 - 144
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 144 - 160
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 160 - 176
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 176 - 192
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 192 - 208
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 208 - 224
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 224 - 240
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 // 240 - 256
};

constexpr std::array<char, 256> CigarNumToChar = {
    'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'
};

// This is intentionally an enum and not an enum class, so we
// can convert the ops to numeric values.
enum SamCigarOperations {
    CIGAR_OP_M = 0,
    CIGAR_OP_I = 1,
    CIGAR_OP_D = 2,
    CIGAR_OP_N = 3,
    CIGAR_OP_S = 4,
    CIGAR_OP_H = 5,
    CIGAR_OP_P = 6,
    CIGAR_OP_EQ = 7,
    CIGAR_OP_X = 8,
};
// clang-format on

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_LOOKUPS_HPP
