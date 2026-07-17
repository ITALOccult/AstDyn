#include "astdyn/core/IOCConfig.hpp"
#include "astdyn/third_party/fkYAML.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>

namespace astdyn::core {

namespace {

/// fkYAML node -> nlohmann::json, recursively.
nlohmann::json yaml_to_json(const fkyaml::node& n) {
    using J = nlohmann::json;
    if (n.is_mapping()) {
        J j = J::object();
        for (const auto& kv : n.as_map()) j[kv.first.as_str()] = yaml_to_json(kv.second);
        return j;
    }
    if (n.is_sequence()) {
        J j = J::array();
        for (const auto& e : n.as_seq()) j.push_back(yaml_to_json(e));
        return j;
    }
    if (n.is_boolean())      return J(n.as_bool());
    if (n.is_integer())      return J(n.as_int());
    if (n.is_float_number()) return J(n.as_float());
    if (n.is_null())         return J(nullptr);
    return J(n.as_str());
}

std::string lowercase_extension(const std::string& filename) {
    const auto dot = filename.find_last_of('.');
    if (dot == std::string::npos) return {};
    std::string ext = filename.substr(dot + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

} // namespace

bool IOCConfig::load(const std::string& filename) {
    if (filename.empty()) return false;

    const std::string ext = lowercase_extension(filename);
    if (ext == "json") return load_json(filename);
    if (ext == "yaml" || ext == "yml") return load_yaml(filename);

    // Anything else is refused, loudly. The previous version fell back to the
    // OOP parser for every unknown extension, which meant a YAML file was
    // accepted and silently ignored: load() returned true, the configuration
    // came out empty, and every get() handed back its default. A configuration
    // that is not read is worse than one that fails to open.
    std::cerr << "IOCConfig: unsupported configuration format '" << ext << "' for "
              << filename << ". Use .yaml, .yml or .json.\n";
    return false;
}

bool IOCConfig::load_json(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open()) return false;
    try {
        f >> data_;
        return true;
    } catch (...) {
        return false;
    }
}

static std::string trim(const std::string& s) {
    size_t first = s.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = s.find_last_not_of(" \t\n\r");
    return s.substr(first, last - first + 1);
}

bool IOCConfig::load_yaml(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "IOCConfig: cannot open " << filename << "\n";
        return false;
    }
    std::stringstream ss;
    ss << f.rdbuf();
    return load_yaml_string(ss.str());
}

bool IOCConfig::load_yaml_string(const std::string& content) {
    try {
        std::istringstream in(content);
        data_ = yaml_to_json(fkyaml::node::deserialize(in));
    } catch (const std::exception& e) {
        std::cerr << "IOCConfig: malformed YAML: " << e.what() << "\n";
        return false;
    }
    if (!data_.is_object()) {
        std::cerr << "IOCConfig: the YAML root must be a mapping\n";
        return false;
    }
    return true;
}


nlohmann::json* IOCConfig::find_node(const std::string& path, bool create) {
    std::stringstream ss(path);
    std::string segment;
    nlohmann::json* current = &data_;

    while (std::getline(ss, segment, '.')) {
        if (segment.empty()) continue;
        
        if (!current->is_object()) {
            if (create) {
                *current = nlohmann::json::object();
            } else {
                return nullptr;
            }
        }

        if (!current->contains(segment)) {
            if (create) {
                (*current)[segment] = nlohmann::json::object();
            } else {
                return nullptr;
            }
        }
        current = &((*current)[segment]);
    }
    return current;
}

const nlohmann::json* IOCConfig::find_node(const std::string& path) const {
    std::stringstream ss(path);
    std::string segment;
    const nlohmann::json* current = &data_;

    while (std::getline(ss, segment, '.')) {
        if (segment.empty()) continue;
        if (!current->is_object() || !current->contains(segment)) {
            return nullptr;
        }
        current = &((*current)[segment]);
    }
    return current;
}

bool IOCConfig::has(const std::string& path) const {
    return find_node(path) != nullptr;
}

std::string IOCConfig::to_json_string(int indent) const {
    return data_.dump(indent);
}

// Explicit instantiations
template std::string IOCConfig::get<std::string>(const std::string&, const std::string&) const;
template double IOCConfig::get<double>(const std::string&, const double&) const;
template int IOCConfig::get<int>(const std::string&, const int&) const;
template bool IOCConfig::get<bool>(const std::string&, const bool&) const;

template void IOCConfig::set<std::string>(const std::string&, const std::string&);
template void IOCConfig::set<double>(const std::string&, const double&);
template void IOCConfig::set<int>(const std::string&, const int&);
template void IOCConfig::set<bool>(const std::string&, const bool&);

} // namespace astdyn::core
