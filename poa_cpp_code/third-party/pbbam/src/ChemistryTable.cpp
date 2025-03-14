// Author: Lance Hepler

#include "PbbamInternalConfig.h"

#include "ChemistryTable.h"

#include <cstdlib>
#include <fstream>
#include <map>

#include "FileUtils.h"
#include "pbbam/exception/BundleChemistryMappingException.h"
#include "pugixml/pugixml.hpp"

namespace PacBio {
namespace BAM {
namespace internal {

extern const ChemistryTable BuiltInChemistryTable = {

    // BindingKit, SequencingKit, BasecallerVersion, Chemistry

    // {{"0.0.500", "0.0.1", "0.0.1", "S/P2-C2/6.0"}}
    {{"0.0.500", "0.0.1", "0.0.1", "Kit__500_Chem__1_BC__1_PW3_v2"}}
  
};

ChemistryTable ChemistryTableFromXml(const std::string& mappingXml)
{
    if (!FileUtils::Exists(mappingXml))
        throw BundleChemistryMappingException{
            mappingXml, "SMRT_CHEMISTRY_BUNDLE_DIR defined but file not found"};

    std::ifstream in(mappingXml);
    pugi::xml_document doc;
    const pugi::xml_parse_result loadResult = doc.load(in);
    if (loadResult.status != pugi::status_ok)
        throw BundleChemistryMappingException{
            mappingXml, "unparseable XML, error code:" + std::to_string(loadResult.status)};

    // parse top-level attributes
    pugi::xml_node rootNode = doc.document_element();
    if (rootNode == pugi::xml_node())
        throw BundleChemistryMappingException{mappingXml, "could not fetch XML root node"};

    if (std::string(rootNode.name()) != "MappingTable")
        throw BundleChemistryMappingException{mappingXml, "MappingTable not found"};

    ChemistryTable table;
    try {
        for (const auto& childNode : rootNode) {
            const std::string childName = childNode.name();
            if (childName != "Mapping") continue;
            table.emplace_back(
                std::array<std::string, 4>{{childNode.child("BindingKit").child_value(),
                                            childNode.child("SequencingKit").child_value(),
                                            childNode.child("SoftwareVersion").child_value(),
                                            childNode.child("SequencingChemistry").child_value()}});
        }
    } catch (std::exception& e) {
        const std::string msg = std::string{"Mapping entries unparseable - "} + e.what();
        throw BundleChemistryMappingException{mappingXml, msg};
    }
    return table;
}

const ChemistryTable& GetChemistryTableFromEnv()
{
    static const ChemistryTable empty{};
    static std::map<std::string, ChemistryTable> tableCache;

    std::string chemPath;
    const char* pth = getenv("SMRT_CHEMISTRY_BUNDLE_DIR");
    if (pth != nullptr && pth[0] != '\0')
        chemPath = pth;
    else
        return empty;

    auto it = tableCache.find(chemPath);
    if (it != tableCache.end()) return it->second;

    auto tbl = ChemistryTableFromXml(chemPath + "/chemistry.xml");
    it = tableCache.emplace(std::move(chemPath), std::move(tbl)).first;
    return it->second;
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio
