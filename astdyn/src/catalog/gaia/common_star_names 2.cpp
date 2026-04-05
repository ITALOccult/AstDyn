#include <astdyn/catalog/gaia/common_star_names.h>
#include <astdyn/catalog/gaia/iau_star_catalog_parser.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

namespace ioc::gaia {

bool CommonStarNames::loadDatabase(const std::string& database_path) {
    std::ifstream file(database_path);
    if (!file.is_open()) {
        std::cerr << "Could not open star names database: " << database_path << std::endl;
        return false;
    }
    
    std::string line;
    size_t line_number = 0;
    size_t loaded_count = 0;
    
    // Clear existing data
    source_id_to_info_.clear();
    name_to_source_id_.clear();
    sao_to_source_id_.clear();
    hd_to_source_id_.clear();
    hip_to_source_id_.clear();
    tycho2_to_source_id_.clear();
    constellation_to_source_ids_.clear();
    
    while (std::getline(file, line)) {
        line_number++;
        
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        CrossMatchInfo info;
        uint64_t source_id;
        
        if (parseLine(line, info, source_id)) {
            source_id_to_info_[source_id] = info;
            addToIndices(source_id, info);
            loaded_count++;
        } else {
            std::cerr << "Warning: Could not parse line " << line_number 
                      << " in " << database_path << std::endl;
        }
    }
    
    std::cout << "Loaded " << loaded_count << " star names from " << database_path << std::endl;
    return loaded_count > 0;
}

bool CommonStarNames::loadDefaultDatabase() {
    // Load embedded database
    loadEmbeddedDatabase();
    return !source_id_to_info_.empty();
}

std::optional<CrossMatchInfo> CommonStarNames::getCrossMatch(uint64_t source_id) const {
    auto it = source_id_to_info_.find(source_id);
    if (it != source_id_to_info_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<uint64_t> CommonStarNames::findByName(const std::string& name) const {
    std::string lower_name = toLower(name);
    auto it = name_to_source_id_.find(lower_name);
    if (it != name_to_source_id_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<uint64_t> CommonStarNames::findBySAO(const std::string& sao_number) const {
    auto it = sao_to_source_id_.find(sao_number);
    if (it != sao_to_source_id_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<uint64_t> CommonStarNames::findByHD(const std::string& hd_number) const {
    auto it = hd_to_source_id_.find(hd_number);
    if (it != hd_to_source_id_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<uint64_t> CommonStarNames::findByHipparcos(const std::string& hip_number) const {
    auto it = hip_to_source_id_.find(hip_number);
    if (it != hip_to_source_id_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<uint64_t> CommonStarNames::findByTycho2(const std::string& tycho2_id) const {
    auto it = tycho2_to_source_id_.find(tycho2_id);
    if (it != tycho2_to_source_id_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::vector<uint64_t> CommonStarNames::getConstellationStars(const std::string& constellation) const {
    std::string lower_constellation = toLower(constellation);
    auto it = constellation_to_source_ids_.find(lower_constellation);
    if (it != constellation_to_source_ids_.end()) {
        return it->second;
    }
    return {};
}

std::vector<std::string> CommonStarNames::getAllNames() const {
    std::vector<std::string> names;
    names.reserve(name_to_source_id_.size());
    
    for (const auto& [name, source_id] : name_to_source_id_) {
        // Find the original case version
        auto info_it = source_id_to_info_.find(source_id);
        if (info_it != source_id_to_info_.end()) {
            names.push_back(info_it->second.common_name);
        }
    }
    
    std::sort(names.begin(), names.end());
    return names;
}

bool CommonStarNames::parseLine(const std::string& line, CrossMatchInfo& info, uint64_t& source_id) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    
    // Split by comma
    while (std::getline(ss, field, ',')) {
        fields.push_back(field);
    }
    
    if (fields.size() < 7) {
        return false; // Need at least 7 fields
    }
    
    try {
        // Parse fields
        source_id = std::stoull(fields[0]);
        info.common_name = fields[1];
        info.sao_designation = fields[2];
        info.hd_designation = fields[3];
        info.hip_designation = fields[4];
        info.tycho2_designation = fields[5];
        info.constellation = fields[6];
        
        if (fields.size() > 7) {
            info.notes = fields[7];
        }
        
        return true;
        
    } catch (const std::exception&) {
        return false;
    }
}

std::string CommonStarNames::toLower(const std::string& str) const {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return result;
}

void CommonStarNames::addToIndices(uint64_t source_id, const CrossMatchInfo& info) {
    // Add to name index (case-insensitive)
    if (!info.common_name.empty()) {
        name_to_source_id_[toLower(info.common_name)] = source_id;
    }
    
    // Add to catalog indices
    if (!info.sao_designation.empty()) {
        sao_to_source_id_[info.sao_designation] = source_id;
    }
    
    if (!info.hd_designation.empty()) {
        hd_to_source_id_[info.hd_designation] = source_id;
    }
    
    if (!info.hip_designation.empty()) {
        hip_to_source_id_[info.hip_designation] = source_id;
    }
    
    if (!info.tycho2_designation.empty()) {
        tycho2_to_source_id_[info.tycho2_designation] = source_id;
    }
    
    // Add to constellation index
    if (!info.constellation.empty()) {
        std::string lower_constellation = toLower(info.constellation);
        constellation_to_source_ids_[lower_constellation].push_back(source_id);
    }
}

void CommonStarNames::loadEmbeddedDatabase() {
    // Embedded star database - brightest and most famous stars
    struct StarData {
        uint64_t source_id;
        const char* name;
        const char* sao;
        const char* hd;
        const char* hip;
        const char* tycho2;
        const char* constellation;
        const char* notes;
    };
    
    static const StarData stars[] = {
        // Brightest stars
        {2652764709496131200ULL, "Sirius", "151881", "48915", "32349", "5496-1026-1", "Canis Major", "Brightest star in night sky"},
        {4472832130942575872ULL, "Canopus", "234480", "45348", "30438", "8538-2822-1", "Carina", "Second brightest star"},
        {5340040567530378624ULL, "Arcturus", "100944", "124897", "69673", "1990-2527-1", "Bootes", "Fourth brightest star"},
        {2618709380455063424ULL, "Vega", "67174", "172167", "91262", "3105-2070-1", "Lyra", "Fifth brightest star"},
        {2483618161171672832ULL, "Capella", "40186", "34029", "24608", "1897-1116-1", "Auriga", "Sixth brightest star"},
        {6897537427869573760ULL, "Rigel", "131907", "34085", "24436", "5331-1752-1", "Orion", "Seventh brightest star"},
        {545827370388883840ULL, "Procyon", "115756", "61421", "37279", "187-2184-1", "Canis Minor", "Eighth brightest star"},
        {7001135808431135616ULL, "Betelgeuse", "113271", "39801", "27989", "4774-1136-1", "Orion", "Variable red supergiant"},
        {2090070885554436992ULL, "Achernar", "293818", "10144", "7588", "9159-684-1", "Eridanus", "Flattest known star"},
        {5453804681476640384ULL, "Altair", "125122", "187642", "97649", "1243-2021-1", "Aquila", "Twelfth brightest star"},
        {5132766955826987136ULL, "Aldebaran", "94027", "29139", "21421", "1266-1416-1", "Taurus", "Bull's eye star"},
        
        // Famous navigation stars
        {724150705952106752ULL, "Polaris", "308", "8890", "11767", "4628-237-1", "Ursa Minor", "North Star"},
        {5829642615214967936ULL, "Deneb", "19014", "197345", "102098", "3138-3060-1", "Cygnus", "Alpha Cygni"},
        {5533833017897066624ULL, "Spica", "157923", "116658", "65474", "6626-1182-1", "Virgo", "Blue giant binary"},
        {4427773526950133504ULL, "Antares", "184415", "148478", "80763", "6839-1741-1", "Scorpius", "Red supergiant"},
        {5853498713160606592ULL, "Pollux", "78716", "62509", "37826", "1985-460-1", "Gemini", "Orange giant"},
        {3555298692598179072ULL, "Regulus", "78416", "87901", "49669", "2467-1066-1", "Leo", "Heart of the Lion"},
        {5066002071299965696ULL, "Fomalhaut", "216956", "216956", "113368", "6697-0952-1", "Piscis Austrinus", "Autumn star"},
        
        // Orion constellation
        {2650825017058833152ULL, "Bellatrix", "112740", "35468", "25336", "4774-1848-1", "Orion", "Amazon Star"},
        {2542462002008550784ULL, "Mintaka", "132346", "36486", "26727", "4774-1588-1", "Orion", "Belt star"},
        {2444913574403008000ULL, "Alnilam", "132223", "37128", "26311", "4774-1160-1", "Orion", "Belt star"},
        {2539905936987806720ULL, "Alnitak", "132301", "37742", "26727", "4774-1659-1", "Orion", "Belt star"},
        
        // Multiple star systems
        {4293678397322932736ULL, "Rigil Kentaurus", "252838", "128620", "71683", "8470-819-1", "Centaurus", "Alpha Centauri A"},
        {5853498713160606592ULL, "Castor", "77984", "60179", "36850", "1985-456-1", "Gemini", "Six-star system"},
        
        // End marker
        {0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
    };
    
    // Clear existing data
    source_id_to_info_.clear();
    name_to_source_id_.clear();
    sao_to_source_id_.clear();
    hd_to_source_id_.clear();
    hip_to_source_id_.clear();
    tycho2_to_source_id_.clear();
    constellation_to_source_ids_.clear();
    
    // Load embedded data
    for (const StarData* star = stars; star->source_id != 0; ++star) {
        CrossMatchInfo info;
        info.common_name = star->name ? star->name : "";
        info.sao_designation = star->sao ? star->sao : "";
        info.hd_designation = star->hd ? star->hd : "";
        info.hip_designation = star->hip ? star->hip : "";
        info.tycho2_designation = star->tycho2 ? star->tycho2 : "";
        info.constellation = star->constellation ? star->constellation : "";
        info.notes = star->notes ? star->notes : "";
        
        source_id_to_info_[star->source_id] = info;
        addToIndices(star->source_id, info);
    }
    
    std::cout << "Loaded " << source_id_to_info_.size() << " stars from embedded database" << std::endl;
}

} // namespace ioc::gaia