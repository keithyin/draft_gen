// Author: Derek Barnett

#ifndef PBCOPPER_QGRAM_INDEXHIT_H
#define PBCOPPER_QGRAM_INDEXHIT_H

#include <cstdint>

#include <pbcopper/PbcopperConfig.h>

namespace PacBio {
namespace QGram {

///
/// \brief The IndexHit class represents a q-gram hit from a query onto input
///        sequence.
///
/// Each instance reports which sequence was hit (useful in the case of
/// multi-sequence indices) and at which position.
///
class IndexHit
{
public:
    ///
    /// \brief IndexHit
    ///
    IndexHit(void) : id_(0), position_(0) {}

    ///
    /// \brief IndexHit
    /// \param[in] id   index number of input sequence
    /// \param[in] pos  position in input sequence
    ///
    IndexHit(const uint32_t id, const uint64_t pos) : id_(id), position_(pos) {}

    IndexHit(const IndexHit&) = default;
    IndexHit(IndexHit&&) = default;
    IndexHit& operator=(const IndexHit&) = default;
    IndexHit& operator=(IndexHit&&) = default;

public:
    ///
    /// \brief Id
    /// \return id (index number) of Index input sequence containing this hit
    ///
    uint32_t Id(void) const { return id_; }

    ///
    /// \brief Position
    /// \return position of hit
    ///
    uint64_t Position(void) const { return position_; }

private:
    uint32_t id_;
    uint64_t position_;
};

///
/// \brief operator ==
/// \param[in] lhs
/// \param[in] rhs
/// \return true if the hits share the same sequence id & position
///
inline bool operator==(const IndexHit& lhs, const IndexHit& rhs)
{
    return lhs.Id() == rhs.Id() && lhs.Position() == rhs.Position();
}

}  // namespace QGram
}  // namespace PacBio

#endif  // PBCOPPER_QGRAM_INDEXHIT_H
