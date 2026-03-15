/**
 * @file ChebyshevEphemerisManager.cpp
 * @brief Implementation of ChebyshevEphemerisManager.
 */

#include "astdyn/astrometry/ChebyshevEphemerisManager.hpp"
#include <stdexcept>

namespace astdyn::astrometry {

void ChebyshevEphemerisManager::add_asteroid(
    const std::string& id,
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
    time::EpochTDB start,
    time::EpochTDB end,
    int degree)
{
    auto ephem = std::make_unique<AsteroidChebyshevEphemeris>(
        initial_elements, start, end, config_, degree
    );
    ephemerides_[id] = std::move(ephem);
}

std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>>
ChebyshevEphemerisManager::evaluate_full(const std::string& id, time::EpochTDB t) const
{
    auto it = ephemerides_.find(id);
    if (it == ephemerides_.end()) {
        throw std::out_of_range("Asteroid " + id + " not managed by ChebyshevEphemerisManager");
    }

    // Use evaluate_full from the underlying segment
    return it->second->evaluate_full(t);
}

} // namespace astdyn::astrometry
