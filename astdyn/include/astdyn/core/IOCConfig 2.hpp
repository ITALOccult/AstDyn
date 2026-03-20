#ifndef ASTDYN_CORE_IOCCONFIG_HPP
#define ASTDYN_CORE_IOCCONFIG_HPP

#include <string>
#include <vector>
#include <map>
#include <variant>
#include <optional>
#include <nlohmann/json.hpp>

namespace astdyn::core {

/**
 * @brief Advanced configuration manager supporting JSON, YAML, and OOP-style formats.
 * Designed after IOC_Config pattern.
 */
class IOCConfig {
public:
    using ConfigValue = std::variant<std::string, double, int, bool, nlohmann::json>;

    IOCConfig() = default;

    /**
     * @brief Load from file (auto-detects format)
     */
    bool load(const std::string& filename);

    /**
     * @brief Load from JSON format
     */
    bool load_json(const std::string& filename);

    /**
     * @brief Load from YAML/OOP format
     */
    bool load_yaml(const std::string& filename);

    /**
     * @brief Load from a string containing OOP-style config (key = value)
     */
    bool load_oop_string(const std::string& content);

    /**
     * @brief Get value with path-style key (e.g., "integrator.step_size")
     */
    template <typename T>
    T get(const std::string& path, const T& default_value) const {
        const nlohmann::json* node = find_node(path);
        if (!node || node->is_null()) {
            return default_value;
        }
        try {
            if constexpr (std::is_same_v<T, std::string>) {
                if (node->is_number()) {
                    return node->dump();
                }
            }
            return node->get<T>();
        } catch (...) {
            return default_value;
        }
    }

    /**
     * @brief Check if a key exists
     */
    bool has(const std::string& path) const;

    /**
     * @brief Set value for a specific path
     */
    template <typename T>
    void set(const std::string& path, const T& value) {
        nlohmann::json* node = find_node(path, true);
        if (node) {
            *node = value;
        }
    }

    /**
     * @brief Clear all data
     */
    void clear() { data_ = nlohmann::json::object(); }

    /**
     * @brief Export current config as JSON string
     */
    std::string to_json_string(int indent = 4) const;

    /**
     * @brief Access internal JSON data
     */
    const nlohmann::json& data() const { return data_; }

private:
    nlohmann::json data_ = nlohmann::json::object();

    // Helper for hierarchical access
    nlohmann::json* find_node(const std::string& path, bool create = false);
    const nlohmann::json* find_node(const std::string& path) const;
};

} // namespace astdyn::core

#endif // ASTDYN_CORE_IOCCONFIG_HPP
