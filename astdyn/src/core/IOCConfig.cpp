#include "astdyn/core/IOCConfig.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>

namespace astdyn::core {

bool IOCConfig::load(const std::string& filename) {
    if (filename.empty()) return false;

    std::string ext = "";
    size_t dot_pos = filename.find_last_of(".");
    if (dot_pos != std::string::npos) {
        ext = filename.substr(dot_pos + 1);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    }

    if (ext == "json") {
        return load_json(filename);
    } else {
        // Assume YAML/OOP for .yaml, .yml, .oop, .cfg, .txt or no extension
        return load_yaml(filename);
    }
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
    if (!f.is_open()) return false;
    
    std::stringstream ss;
    ss << f.rdbuf();
    return load_oop_string(ss.str());
}

bool IOCConfig::load_oop_string(const std::string& content) {
    std::stringstream ss(content);
    std::string line;
    std::vector<std::string> scope;

    while (std::getline(ss, line)) {
        // Remove comments
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);
        
        line = trim(line);
        if (line.empty()) continue;

        // Handle scope entry/exit
        if (line.back() == '{') {
            std::string name = trim(line.substr(0, line.length() - 1));
            if (!name.empty()) scope.push_back(name);
            continue;
        }
        if (line == "}") {
            if (!scope.empty()) scope.pop_back();
            continue;
        }

        // Handle key = value
        size_t eq_pos = line.find('=');
        if (eq_pos != std::string::npos) {
            std::string key = trim(line.substr(0, eq_pos));
            std::string val = trim(line.substr(eq_pos + 1));
            
            // Build full path
            std::string full_key = "";
            for (const auto& s : scope) full_key += s + ".";
            full_key += key;

            // Type deduction
            if (val == "true" || val == "True") {
                set(full_key, true);
            } else if (val == "false" || val == "False") {
                set(full_key, false);
            } else {
                // Try number
                try {
                    size_t pos;
                    double d = std::stod(val, &pos);
                    if (pos == val.length()) {
                        // It's a number. Check if it's integer for better precision
                        if (val.find('.') == std::string::npos && val.find('e') == std::string::npos && val.find('E') == std::string::npos) {
                            set(full_key, std::stoll(val));
                        } else {
                            set(full_key, d);
                        }
                    } else {
                        set(full_key, val);
                    }
                } catch (...) {
                    set(full_key, val);
                }
            }
        }
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
