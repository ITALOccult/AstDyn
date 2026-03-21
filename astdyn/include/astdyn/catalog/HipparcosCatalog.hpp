/**
 * @file HipparcosCatalog.hpp
 * @brief Hipparcos stellar catalog interface (stub implementation)
 *
 * The Hipparcos catalog contains ~118,000 stars with highly accurate
 * astrometry (parallax, proper motion) in the ICRS frame. Useful for
 * occultation work with bright stars (V < 9) and for astrometric calibration.
 *
 * This is a stub that implements the StellarCatalog interface.
 * A full implementation requires the Hipparcos catalog files available from:
 * https://cdsarc.cds.unistra.fr/viz-bin/cat/I/239
 *
 * @note Call initialize(path_to_hipparcos_dir) before use.
 */

#ifndef ASTDYN_CATALOG_HIPPARCOS_CATALOG_HPP
#define ASTDYN_CATALOG_HIPPARCOS_CATALOG_HPP

#include "astdyn/catalog/StellarCatalog.hpp"
#include <string>
#include <stdexcept>

namespace astdyn::catalog {

struct HipparcosNotInitialized : std::runtime_error {
    HipparcosNotInitialized() : std::runtime_error("HipparcosCatalog: not initialized — call initialize() first") {}
};

/**
 * @brief Hipparcos stellar catalog (stub).
 *
 * Usage:
 * @code
 *   HipparcosCatalog::initialize("/path/to/hipparcos/");
 *   auto& cat = HipparcosCatalog::instance();
 *   auto stars = cat.query_cone(query);
 * @endcode
 */
class HipparcosCatalog : public StellarCatalog {
public:
    static void initialize(const std::string& data_dir) {
        auto& inst = get_instance();
        inst.data_dir_ = data_dir;
        inst.initialized_ = true;
    }

    static HipparcosCatalog& instance() {
        auto& inst = get_instance();
        if (!inst.initialized_) throw HipparcosNotInitialized{};
        return inst;
    }

    static void shutdown() noexcept {
        auto& inst = get_instance();
        inst.initialized_ = false;
        inst.data_dir_.clear();
    }

    std::vector<Star> query_cone(const ConeQuery& params) const override {
        check_initialized();
        // TODO: implement hip_main.dat sequential or indexed lookup
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
        // source_id encodes the HIP number directly
        // TODO: binary search in sorted hip_main.dat
        return std::nullopt;
    }

    std::optional<Star> by_name(const std::string& name) const override {
        check_initialized();
        // TODO: lookup by HD/Bayer/Flamsteed designation via cross-index
        return std::nullopt;
    }

    const std::string& data_dir() const noexcept { return data_dir_; }
    bool is_initialized() const noexcept { return initialized_; }

private:
    HipparcosCatalog() = default;

    static HipparcosCatalog& get_instance() {
        static HipparcosCatalog instance;
        return instance;
    }

    void check_initialized() const {
        if (!initialized_) throw HipparcosNotInitialized{};
    }

    std::string data_dir_;
    bool initialized_ = false;
};

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_HIPPARCOS_CATALOG_HPP
