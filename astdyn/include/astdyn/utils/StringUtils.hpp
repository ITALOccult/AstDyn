/**
 * @file StringUtils.hpp
 * @brief String manipulation utilities
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Replaces Fortran char_str.f90 module functionality.
 */

#ifndef ASTDYN_UTILS_STRINGUTILS_HPP
#define ASTDYN_UTILS_STRINGUTILS_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

namespace astdyn {
namespace utils {

/**
 * @brief Trim whitespace from both ends of string
 */
inline std::string trim(const std::string& str) {
    auto start = std::find_if_not(str.begin(), str.end(),
                                   [](unsigned char ch) { return std::isspace(ch); });
    auto end = std::find_if_not(str.rbegin(), str.rend(),
                                 [](unsigned char ch) { return std::isspace(ch); }).base();
    return (start < end) ? std::string(start, end) : std::string();
}

/**
 * @brief Convert string to lowercase
 */
inline std::string to_lower(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return str;
}

/**
 * @brief Convert string to uppercase
 */
inline std::string to_upper(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return str;
}

/**
 * @brief Split string by delimiter
 */
inline std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream stream(str);
    
    while (std::getline(stream, token, delimiter)) {
        tokens.push_back(token);
    }
    
    return tokens;
}

/**
 * @brief Check if string contains only digits
 */
inline bool is_numeric(const std::string& str) {
    return !str.empty() && std::all_of(str.begin(), str.end(),
                                        [](unsigned char c) { return std::isdigit(c); });
}

/**
 * @brief Check if string contains only letters
 */
inline bool is_alpha(const std::string& str) {
    return !str.empty() && std::all_of(str.begin(), str.end(),
                                        [](unsigned char c) { return std::isalpha(c); });
}

/**
 * @brief Replace all occurrences of substring
 */
inline std::string replace_all(std::string str, const std::string& from, const std::string& to) {
    size_t pos = 0;
    while ((pos = str.find(from, pos)) != std::string::npos) {
        str.replace(pos, from.length(), to);
        pos += to.length();
    }
    return str;
}

/**
 * @brief Check if string starts with prefix
 */
inline bool starts_with(const std::string& str, const std::string& prefix) {
    return str.size() >= prefix.size() &&
           str.compare(0, prefix.size(), prefix) == 0;
}

/**
 * @brief Check if string ends with suffix
 */
inline bool ends_with(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

/**
 * @brief Pad string to specified length with character
 */
inline std::string pad_left(const std::string& str, size_t width, char fill = ' ') {
    return (str.size() < width) ? std::string(width - str.size(), fill) + str : str;
}

inline std::string pad_right(const std::string& str, size_t width, char fill = ' ') {
    return (str.size() < width) ? str + std::string(width - str.size(), fill) : str;
}

} // namespace utils
} // namespace astdyn

#endif // ASTDYN_UTILS_STRINGUTILS_HPP
