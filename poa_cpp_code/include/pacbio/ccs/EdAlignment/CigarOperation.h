#pragma once

#include "Config.h"

#include <stdexcept>
#include <type_traits>

#include <cassert>
#include <cstdint>

//#include <string_view>

using namespace std;

namespace Geneus {
namespace Data {

/// \brief Describes a CIGAR operation.
///
/// Bracketed character is the corresponding SAM/BAM character code.
///
/// \warning ALIGNMENT_MATCH ('M') is included in this enum to maintain
///          consistency with htslib. However, as of PacBio BAM spec version
///          3.0b7, this CIGAR operation \b forbidden. Any attempt to read or
///          write a record containing this operation will trigger a
///          std::runtime_error. SEQUENCE_MATCH('=') or SEQUENCE_MISMATCH('X')
///          should be used instead.
///
enum class CigarOperationType
{
    ALIGNMENT_MATCH = 0,  ///< alignment match (can be a sequence match or mismatch) [M]
    INSERTION,            ///< insertion to the reference [I]
    DELETION,             ///< deletion from the reference [D]
    REFERENCE_SKIP,       ///< skipped region from the reference [N]
    SOFT_CLIP,            ///< soft clipping (clipped sequences present in SEQ) [S]
    HARD_CLIP = 5,        ///< hard clipping (clipped sequences NOT present in SEQ) [H]
    PADDING,              ///< padding (silent deletion from padded reference) [P]
    SEQUENCE_MATCH,       ///< sequence match [=]
    SEQUENCE_MISMATCH,    ///< sequence mismatch [X]
    UNKNOWN_OP = 15,      ///< unknown/invalid CIGAR operator
};

/// \brief The CigarOperation class represents a single CIGAR operation
///        (consisting of a type & length).
///
class CigarOperation
{
public:
    //static void EnableAutoValidation();
    static void EnableAutoValidation() { AutoValidateCigar = true; }

    //static void DisableAutoValidation();
    static void DisableAutoValidation() { AutoValidateCigar = false; }

public:
    /// \name Operation Type Conversion Methods
    /// \{

    /// Convert between CigarOperationType enum & SAM/BAM character code.
    ///
    /// \param[in] type CigarOperationType value
    /// \returns SAM/BAM character code
    //static char TypeToChar(CigarOperationType type);
    static char TypeToChar(const CigarOperationType type)
    {
        string LOOKUP_TABLE{ "MIDNSHP=XB??????" };
        return LOOKUP_TABLE[static_cast<int32_t>(type)];
    }

    /// Convert between CigarOperationType enum & SAM/BAM character code.
    ///
    /// \param[in] c SAM/BAM character code
    /// \returns CigarOperationType value
    PB_CUDA_HOST PB_CUDA_DEVICE constexpr static CigarOperationType CharToType(const char c)
    {
        switch (c) {
            case 'S':
                return CigarOperationType::SOFT_CLIP;
            case '=':
                return CigarOperationType::SEQUENCE_MATCH;
            case 'X':
                return CigarOperationType::SEQUENCE_MISMATCH;
            case 'I':
                return CigarOperationType::INSERTION;
            case 'D':
                return CigarOperationType::DELETION;
            case 'N':
                return CigarOperationType::REFERENCE_SKIP;
            case 'H':
                return CigarOperationType::HARD_CLIP;
            case 'P':
                return CigarOperationType::PADDING;
            case 'M':
                return CigarOperationType::ALIGNMENT_MATCH;
            default:
                return CigarOperationType::UNKNOWN_OP;
        }
    }

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    CigarOperation() = default;

    PB_CUDA_HOST PB_CUDA_DEVICE
#ifdef __cpp_lib_is_constant_evaluated
        constexpr
#endif
        CigarOperation(const char c, const uint32_t length)
        : CigarOperation{CigarOperation::CharToType(c), length}
    {}

    PB_CUDA_HOST PB_CUDA_DEVICE
#ifdef __cpp_lib_is_constant_evaluated
        constexpr
#endif
        CigarOperation(const CigarOperationType op, const uint32_t length)
        : data_{(length << 4) | static_cast<uint32_t>(op)}
    {
#ifndef __CUDA_ARCH__  // host
        if (
#ifdef __cpp_lib_is_constant_evaluated
            !is_constant_evaluated() &&
#endif
            AutoValidateCigar && (Type() == CigarOperationType::ALIGNMENT_MATCH)) {
            throw std::runtime_error{
                "[pbcopper] CIGAR operation ERROR: 'M' is not allowed in PacBio BAM files. Use "
                "'X/=' instead."};
        }
#endif
    }

    /// \}

public:
    /// \returns operation type as SAM/BAM char code
    //char Char() const;
    char Char() const { return CigarOperation::TypeToChar(Type()); }
    /// \returns operation length
    PB_CUDA_HOST PB_CUDA_DEVICE constexpr uint32_t Length() const noexcept
    {
        return data_ >> 4;
    }

    /// \returns operation type as CigarOperationType enum value
    PB_CUDA_HOST PB_CUDA_DEVICE constexpr CigarOperationType Type() const noexcept
    {
        return static_cast<CigarOperationType>(data_ & 0b1111U);
    }

    /// \}

public:
    /// \name Attributes
    /// \{

    /// Sets this operation type.
    ///
    /// \param[in] opChar SAM/BAM character code
    /// \returns reference to this operation
    constexpr CigarOperation& Char(const char opChar) noexcept
    {
        return Type(CigarOperation::CharToType(opChar));
    }

    /// Sets this operation length.
    ///
    /// \param[in] length
    /// \returns reference to this operation
    constexpr CigarOperation& Length(const std::uint32_t length) noexcept
    {
        data_ = (data_ & 0b1111U) | (length << 4);
        return *this;
    }

    /// Sets this operation type.
    ///
    /// \param[in] opType CigarOperationType value
    /// \returns reference to this operation
    constexpr CigarOperation& Type(const CigarOperationType opType) noexcept
    {
        data_ = (data_ & ~0b1111U) | static_cast<uint32_t>(opType);
        return *this;
    }

    /// \}

public:
    /// \name Comparison Operators
    /// \{

    /// \returns true if both CIGAR operation type & length match
    constexpr bool operator==(const CigarOperation& rhs) const noexcept
    {
        return data_ == rhs.data_;
    }

    /// \returns true if either CIGAR operation type or length differ
    constexpr bool operator!=(const CigarOperation& rhs) const noexcept { return !(*this == rhs); }

    /// \}

private:
    static bool AutoValidateCigar;// = true;

private:
    // we use the same representation as the SAM/BAM spec
    // 4 least significant bits: the operation
    // 28 most significant bits: the length

    uint32_t data_ = static_cast<uint32_t>(CigarOperationType::UNKNOWN_OP);
};
bool CigarOperation::AutoValidateCigar = true;
//bool ConsumesQuery(CigarOperationType type) noexcept;
bool ConsumesQuery(const CigarOperationType type) noexcept
{
    //                                X=PHSNDIM
    constexpr uint32_t LOOKUP_TABLE{ 0b110010011 };
    const auto val = static_cast<uint32_t>(type);
    assert(val <= 8);

    return (LOOKUP_TABLE >> val) & 0b1U;
}
//bool ConsumesReference(CigarOperationType type) noexcept;

bool ConsumesReference(const CigarOperationType type) noexcept
{
    //                                X=PHSNDIM
    constexpr uint32_t LOOKUP_TABLE{ 0b110001101 };
    const auto val = static_cast<uint32_t>(type);
    assert(val <= 8);

    return (LOOKUP_TABLE >> val) & 0b1U;
}

}  // namespace Data
}  // namespace PacBio


