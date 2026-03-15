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
    ephemerides_[id] = std::make_unique<AsteroidChebyshevEphemeris>(
        initial_elements, start, end, config_, degree
    );
}

void ChebyshevEphemerisManager::add_system_body(
    const std::string& id,
    int naif_id,
    io::SPKReader& reader,
    time::EpochTDB start,
    time::EpochTDB end,
    int degree)
{
    ephemerides_[id] = std::make_unique<SPKChebyshevEphemeris>(
        naif_id, reader, start, end, degree
    );
}

std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>>
ChebyshevEphemerisManager::evaluate_full(const std::string& id, time::EpochTDB t) const
{
    auto it = ephemerides_.find(id);
    if (it == ephemerides_.end()) {
        throw std::out_of_range("Body " + id + " not managed by ChebyshevEphemerisManager");
    }

    return it->second->evaluate_full(t);
}

const IChebyshevEphemeris& ChebyshevEphemerisManager::get_ephemeris(const std::string& id) const {
    auto it = ephemerides_.find(id);
    if (it == ephemerides_.end()) {
        throw std::out_of_range("Ephemeris for " + id + " not managed");
    }
    return *(it->second);
}

} // namespace astdyn::astrometry
