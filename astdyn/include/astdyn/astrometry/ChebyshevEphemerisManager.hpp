#ifndef ASTDYN_ASTROMETRY_CHEBYSHEV_EPHEMERIS_MANAGER_HPP
#define ASTDYN_ASTROMETRY_CHEBYSHEV_EPHEMERIS_MANAGER_HPP

#include "astdyn/astrometry/IChebyshevEphemeris.hpp"
#include "astdyn/astrometry/AsteroidChebyshevEphemeris.hpp"
#include "astdyn/astrometry/SPKChebyshevEphemeris.hpp"
#include <map>
#include <string>
#include <memory>
#include <vector>

namespace astdyn::astrometry {

/**
 * @brief High-level manager to handle Chebyshev polynomials for a collection of asteroids.
 *
 * This class coordinates the pre-calculation and retrieval of ephemerides for multiple
 * bodies (from Keplerian elements or SPK files), optimizing the transition from 
 * numerical propagation to analytical interpolation.
 */
class ChebyshevEphemerisManager {
public:
    ChebyshevEphemerisManager(const AstDynConfig& config) : config_(config) {}

    /**
     * @brief Add an asteroid to the managed list using Keplerian elements.
     */
    void add_asteroid(
        const std::string& id,
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
        time::EpochTDB start,
        time::EpochTDB end,
        int degree = 12
    );

    /**
     * @brief Add a body from an SPK file.
     * 
     * @param id Human-readable ID.
     * @param naif_id NAIF integer ID.
     * @param reader Initialized SPKReader.
     * @param start Start epoch.
     * @param end End epoch.
     * @param degree Polynomial degree.
     */
    void add_system_body(
        const std::string& id,
        int naif_id,
        io::SPKReader& reader,
        time::EpochTDB start,
        time::EpochTDB end,
        int degree = 12
    );

    /**
     * @brief Evaluate position and velocity for a specific body.
     */
    std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>>
    evaluate_full(const std::string& id, time::EpochTDB t) const;

    /**
     * @brief Check if a body is managed.
     */
    bool has_body(const std::string& id) const {
        return ephemerides_.find(id) != ephemerides_.end();
    }

    void clear() { ephemerides_.clear(); }
    size_t size() const { return ephemerides_.size(); }

    /**
     * @brief Get underlying ephemeris.
     */
    const IChebyshevEphemeris& get_ephemeris(const std::string& id) const;

    /**
     * @brief Set body diameter.
     */
    void set_diameter(const std::string& id, double diameter_km) {
        diameters_km_[id] = diameter_km;
    }

    /**
     * @brief Get body diameter (defaults to 100km if unknown).
     */
    double get_diameter(const std::string& id) const {
        auto it = diameters_km_.find(id);
        return (it != diameters_km_.end()) ? it->second : 100.0;
    }

private:
    const AstDynConfig& config_;
    std::map<std::string, std::unique_ptr<IChebyshevEphemeris>> ephemerides_;
    std::map<std::string, double> diameters_km_;
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_CHEBYSHEV_EPHEMERIS_MANAGER_HPP
