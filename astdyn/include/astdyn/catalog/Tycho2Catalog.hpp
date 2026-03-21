/**
 * @file Tycho2Catalog.hpp
 * @brief Tycho-2 stellar catalog interface (stub implementation)
 *
 * Tycho-2 contains ~2.5 million stars with positions, proper motions,
 * and BT/VT photometry. Useful for occultation work involving bright stars
 * (V < 12) not well-covered by Gaia DR3 due to saturation.
 *
 * This is a stub that implements the StellarCatalog interface.
 * A full implementation requires the Tycho-2 catalog files available from:
 * https://cdsarc.cds.unistra.fr/viz-bin/cat/I/259
 *
 * @note Call initialize(path_to_tycho2_dir) before use.
 */

#ifndef ASTDYN_CATALOG_TYCHO2_CATALOG_HPP
#define ASTDYN_CATALOG_TYCHO2_CATALOG_HPP

#include "astdyn/catalog/StellarCatalog.hpp"
#include <string>
#include <stdexcept>

namespace astdyn::catalog {

struct Tycho2NotInitialized : std::runtime_error {
    Tycho2NotInitialized() : std::runtime_error("Tycho2Catalog: not initialized — call initialize() first") {}
};

/**
 * @brief Tycho-2 stellar catalog (stub).
 *
 * Usage:
 * @code
 *   Tycho2Catalog::initialize("/path/to/tycho2/");
 *   auto& cat = Tycho2Catalog::instance();
 *   auto stars = cat.query_cone(query);
 * @endcode
 */
class Tycho2Catalog : public StellarCatalog {
public:
    static void initialize(const std::string& data_dir) {
        auto& inst = get_instance();
        inst.data_dir_ = data_dir;
        inst.initialized_ = true;
    }

    static Tycho2Catalog& instance() {
        auto& inst = get_instance();
        if (!inst.initialized_) throw Tycho2NotInitialized{};
        return inst;
    }

    static void shutdown() noexcept {
        auto& inst = get_instance();
        inst.initialized_ = false;
        inst.data_dir_.clear();
    }

    std::vector<Star> query_cone(const ConeQuery& params) const override {
        check_initialized();
        // TODO: implement main table (tyc2.dat) lookup by HEALPix zone
        return {};
    }

    std::vector<Star> query_corridor(const CorridorQuery& params) const override {
        check_initialized();
        // TODO: implement corridor search
        return {};
    }

    std::vector<Star> query_orbit(const OrbitQuery& params, ProgressCallback progress = nullptr) const override {
        check_initialized();
        // TODO: implement time-windowed orbit corridor search
        return {};
    }

    std::optional<Star> by_source_id(int64_t source_id) const override {
        check_initialized();
        // TODO: decode Tycho-2 designation from source_id (TYC1-TYC2-TYC3)
        return std::nullopt;
    }

    std::optional<Star> by_name(const std::string& name) const override {
        check_initialized();
        // TODO: cross-match via Hipparcos number or HD designation
        return std::nullopt;
    }

    const std::string& data_dir() const noexcept { return data_dir_; }
    bool is_initialized() const noexcept { return initialized_; }

private:
    Tycho2Catalog() = default;

    static Tycho2Catalog& get_instance() {
        static Tycho2Catalog instance;
        return instance;
    }

    void check_initialized() const {
        if (!initialized_) throw Tycho2NotInitialized{};
    }

    std::string data_dir_;
    bool initialized_ = false;
};

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_TYCHO2_CATALOG_HPP
