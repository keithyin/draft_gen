// File Description
/// \file ReadGroupInfo.inl
/// \brief Inline implementations for the ReadGroupInfo class.
//
// Author: Derek Barnett

#include <stdexcept>
#include "pbbam/ReadGroupInfo.h"

namespace PacBio {
namespace BAM {

inline size_t ReadGroupInfo::BarcodeCount() const
{
    if (!hasBarcodeData_)
        throw std::runtime_error{"barcode count requested but barcode data is missing"};
    return barcodeCount_;
}

inline ReadGroupInfo& ReadGroupInfo::BarcodeData(std::string barcodeFile,
                                                 std::string barcodeHash,
                                                 size_t barcodeCount,
                                                 BarcodeModeType barcodeMode,
                                                 BarcodeQualityType barcodeQuality)
{
    barcodeFile_ = std::move(barcodeFile);
    barcodeHash_ = std::move(barcodeHash);
    barcodeCount_ = barcodeCount;
    barcodeMode_ = barcodeMode;
    barcodeQuality_ = barcodeQuality;
    hasBarcodeData_ = true;
    return *this;
}

inline std::string ReadGroupInfo::BarcodeFile() const
{
    if (!hasBarcodeData_)
        throw std::runtime_error{"barcode file requested but barcode data is missing"};
    return barcodeFile_;
}

inline std::string ReadGroupInfo::BarcodeHash() const
{
    if (!hasBarcodeData_)
        throw std::runtime_error{"barcode hash requested but barcode data is missing"};
    return barcodeHash_;
}

inline BarcodeModeType ReadGroupInfo::BarcodeMode() const
{
    if (!hasBarcodeData_)
        throw std::runtime_error{"barcode mode requested but barcode data is missing"};
    return barcodeMode_;
}

inline BarcodeQualityType ReadGroupInfo::BarcodeQuality() const
{
    if (!hasBarcodeData_)
        throw std::runtime_error{"barcode quality requested but barcode data is missing"};
    return barcodeQuality_;
}

inline std::string ReadGroupInfo::BasecallerVersion() const
{ return basecallerVersion_; }

inline ReadGroupInfo& ReadGroupInfo::BasecallerVersion(std::string versionNumber)
{
    if (basecallerVersion_ != versionNumber) { 
        basecallerVersion_ = std::move(versionNumber);
        sequencingChemistry_.clear(); // reset cached chemistry name
    }
    return *this; 
}

inline std::string ReadGroupInfo::BaseFeatureTag(BaseFeature feature) const
{
    const auto iter = features_.find(feature);
    if (iter == features_.end())
        return {};
    return iter->second;
}

inline ReadGroupInfo& ReadGroupInfo::BaseFeatureTag(BaseFeature feature,
                                                    std::string tag)
{ features_[feature] = std::move(tag); return *this; }

inline std::string ReadGroupInfo::BindingKit() const
{ return bindingKit_; }

inline std::string ReadGroupInfo::Model() const {
    return model_;
}

inline ReadGroupInfo& ReadGroupInfo::BindingKit(std::string kitNumber)
{
    if (bindingKit_ != kitNumber) { 
        bindingKit_ = std::move(kitNumber);
        sequencingChemistry_.clear(); // reset cached chemistry name
    }
    return *this; 
}

inline ReadGroupInfo& ReadGroupInfo::ClearBarcodeData()
{
    barcodeFile_.clear();
    barcodeHash_.clear();
    hasBarcodeData_ = false;
    return *this;
}

inline ReadGroupInfo& ReadGroupInfo::ClearBaseFeatures()
{
    features_.clear();
    return *this;
}

inline bool ReadGroupInfo::Control() const
{ return control_; }

inline ReadGroupInfo& ReadGroupInfo::Control(bool ctrl)
{ control_ = ctrl; return *this; }

inline std::map<std::string, std::string> ReadGroupInfo::CustomTags() const
{ return custom_; }

inline ReadGroupInfo& ReadGroupInfo::CustomTags(std::map<std::string, std::string> custom)
{ custom_ = std::move(custom); return *this; }

inline std::string ReadGroupInfo::Date() const
{ return date_; }

inline ReadGroupInfo& ReadGroupInfo::Date(std::string date)
{ date_ = std::move(date); return *this; }

inline std::string ReadGroupInfo::FlowOrder() const
{ return flowOrder_; }

inline ReadGroupInfo& ReadGroupInfo::FlowOrder(std::string order)
{ flowOrder_ = std::move(order); return *this; }

inline std::string ReadGroupInfo::FrameRateHz() const
{ return frameRateHz_; }

inline ReadGroupInfo& ReadGroupInfo::FrameRateHz(std::string frameRateHz)
{ frameRateHz_ = std::move(frameRateHz); return *this; }

inline bool ReadGroupInfo::HasBarcodeData() const
{ return hasBarcodeData_; }

inline bool ReadGroupInfo::HasBaseFeature(BaseFeature feature) const
{ return features_.find(feature) != features_.end(); }

inline std::string ReadGroupInfo::Id() const
{ return id_; }

inline ReadGroupInfo& ReadGroupInfo::Id(std::string id)
{ id_ = std::move(id); return *this; }

inline ReadGroupInfo& ReadGroupInfo::Id(const std::string& movieName,
                                        const std::string& readType)
{ id_ = MakeReadGroupId(movieName, readType); return *this; }

inline int32_t ReadGroupInfo::IdToInt(const std::string& rgId)
{
    const uint32_t rawid = std::stoul(rgId, nullptr, 16);
    return static_cast<int32_t>(rawid);
}

inline FrameCodec ReadGroupInfo::IpdCodec() const
{ return ipdCodec_; }

inline bool ReadGroupInfo::IsValid() const
{ return !id_.empty(); }

inline std::string ReadGroupInfo::KeySequence() const
{ return keySequence_; }

inline ReadGroupInfo& ReadGroupInfo::KeySequence(std::string sequence)
{ keySequence_ = std::move(sequence); return *this; }

inline std::string ReadGroupInfo::Library() const
{ return library_; }

inline ReadGroupInfo& ReadGroupInfo::Library(std::string library)
{ library_ = std::move(library); return *this; }

inline std::string ReadGroupInfo::MovieName() const
{ return movieName_; }

inline ReadGroupInfo& ReadGroupInfo::MovieName(std::string movieName)
{ movieName_ = std::move(movieName); return *this; }

inline std::string ReadGroupInfo::Platform() const
{ return std::string("GENEUS"); }

inline PlatformModelType ReadGroupInfo::PlatformModel() const
{ return platformModel_; }

inline ReadGroupInfo& ReadGroupInfo::PlatformModel(PlatformModelType platform)
{ platformModel_ = platform; return *this; }

inline std::string ReadGroupInfo::PredictedInsertSize() const
{ return predictedInsertSize_; }

inline ReadGroupInfo& ReadGroupInfo::PredictedInsertSize(std::string size)
{ predictedInsertSize_ = std::move(size); return *this; }

inline std::string ReadGroupInfo::Programs() const
{ return programs_; }

inline ReadGroupInfo& ReadGroupInfo::Programs(std::string programs)
{ programs_ = std::move(programs); return *this; }

inline FrameCodec ReadGroupInfo::PulseWidthCodec() const
{ return pulseWidthCodec_; }

inline std::string ReadGroupInfo::ReadType() const
{ return readType_; }

inline ReadGroupInfo& ReadGroupInfo::ReadType(std::string type)
{ readType_ = std::move(type); return *this; }

inline ReadGroupInfo& ReadGroupInfo::RemoveBaseFeature(BaseFeature feature)
{
    const auto iter = features_.find(feature);
    if (iter != features_.end())
        features_.erase(iter);
    return *this;
}

inline std::string ReadGroupInfo::Sample() const
{ return sample_; }

inline ReadGroupInfo& ReadGroupInfo::Sample(std::string sample)
{ sample_ = std::move(sample); return *this; }

inline std::string ReadGroupInfo::SequencingCenter() const
{ return sequencingCenter_; }

inline ReadGroupInfo& ReadGroupInfo::SequencingCenter(std::string center)
{ sequencingCenter_ = std::move(center); return *this; }

inline std::string ReadGroupInfo::SequencingChemistry() const
{   
    if (!sequencingChemistry_.empty()) return sequencingChemistry_;
    if (!model_.empty()) return sequencingChemistry_ = model_;
    return sequencingChemistry_ = SequencingChemistryFromTriple(BindingKit(),
                                                                SequencingKit(),
                                                                BasecallerVersion());
}

inline std::string ReadGroupInfo::SequencingKit() const
{ return sequencingKit_; }

inline ReadGroupInfo& ReadGroupInfo::SequencingKit(std::string kitNumber)
{ 
    if (sequencingKit_ != kitNumber) {
        sequencingKit_ = std::move(kitNumber);
        sequencingChemistry_.clear(); // reset cached chemistry name
    }
    return *this; }

inline std::string ReadGroupInfo::ToSam(const ReadGroupInfo& rg)
{ return rg.ToSam(); }

} // namespace BAM
} // namespace PacBio
