/**
 * @file UCAC4Catalog.hpp
 * @brief UCAC4 stellar catalog interface (stub implementation)
 *
 * UCAC4 (USNO CCD Astrograph Catalog, 4th release) contains ~113 million
 * stars with positions, proper motions, and photometry. It is particularly
 * useful for occultation work involving bright stars not well-covered by
 * Gaia DR3 (e.g., stars brighter than G ≈ 13).
 *
 * This is a stub that implements the StellarCatalog interface.
 * A full implementation requires the UCAC4 binary data files available from:
 * https://cdsarc.cds.unistra.fr/viz-bin/cat/I/322A
 *
 * @note Populate by calling initialize(path_to_ucac4_dir) before use.
 */

#ifndef ASTDYN_CATALOG_UCAC4_CATALOG_HPP
#define ASTDYN_CATALOG_UCAC4_CATALOG_HPP

#include "astdyn/catalog/StellarCatalog.hpp"
#include <string>
#include <stdexcept>

namespace astdyn::catalog {

/**
 * @brief Exception thrown when the UCAC4 catalog has not been initialized.
 */
struct UCAC4NotInitialized : std::runtime_error {
    UCAC4NotInitialized() : std::runtime_error("UCAC4Catalog: not initialized — call initialize() first") {}
};

/**
 * @brief UCAC4 stellar catalog (stub).
 *
 * Usage:
 * @code
 *   UCAC4Catalog::initialize("/path/to/ucac4/");
 *   auto& cat = UCAC4Catalog::instance();
 *   auto stars = cat.query_cone(query);
 * @endcode
 */
class UCAC4Catalog : public StellarCatalog {
public:
    /**
     * @brief Initialize the catalog from the UCAC4 data directory.
     * @param data_dir Path to the directory containing the UCAC4 zone files (u4b/).
     * @throws std::runtime_error if the directory cannot be opened.
     */
    static void initialize(const std::string& data_dir) {
        auto& inst = get_instance();
        inst.data_dir_ = data_dir;
        inst.initialized_ = true;
    }

    /**
     * @brief Returns the singleton instance.
     * @throws UCAC4NotInitialized if initialize() has not been called.
     */
    static UCAC4Catalog& instance() {
        auto& inst = get_instance();
        if (!inst.initialized_) throw UCAC4NotInitialized{};
        return inst;
    }

    /**
     * @brief Release resources.
     */
    static void shutdown() noexcept {
        auto& inst = get_instance();
        inst.initialized_ = false;
        inst.data_dir_.clear();
    }

    // =========================================================================
    // StellarCatalog interface
    // =========================================================================

    std::vector<Star> query_cone(const ConeQuery& params) const override {
        check_initialized();
        // TODO: implement binary zone-file lookup
        return {};
    }

    std::vector<Star> query_corridor(const CorridorQuery& params) const override {
        check_initialized();
        // TODO: implement corridor search over zone files
        return {};
    }

    std::vector<Star> query_orbit(const OrbitQuery& params, ProgressCallback progress = nullptr) const override {
        check_initialized();
        // TODO: implement time-windowed orbit corridor search
        return {};
    }

    std::optional<Star> by_source_id(int64_t source_id) const override {
        check_initialized();
        // TODO: direct lookup by UCAC4 zone+record number encoded in source_id
        return std::nullopt;
    }

    std::optional<Star> by_name(const std::string& name) const override {
        check_initialized();
        // TODO: lookup cross-identifications (HD, HIP, Tycho-2)
        return std::nullopt;
    }

    /**
     * @brief Returns the path to the UCAC4 data directory.
     */
    const std::string& data_dir() const noexcept { return data_dir_; }

    bool is_initialized() const noexcept { return initialized_; }

private:
    UCAC4Catalog() = default;

    static UCAC4Catalog& get_instance() {
        static UCAC4Catalog instance;
        return instance;
    }

    void check_initialized() const {
        if (!initialized_) throw UCAC4NotInitialized{};
    }

    std::string data_dir_;
    bool initialized_ = false;
};

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_UCAC4_CATALOG_HPP
