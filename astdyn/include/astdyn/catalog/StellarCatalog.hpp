/**
 * @file StellarCatalog.hpp
 * @brief Abstract interface for stellar catalogs.
 */

#ifndef ASTDYN_CATALOG_STELLAR_CATALOG_HPP
#define ASTDYN_CATALOG_STELLAR_CATALOG_HPP

#include "astdyn/catalog/CatalogTypes.hpp"
#include <string>
#include <vector>
#include <optional>

namespace astdyn::catalog {

/**
 * @brief Abstract base class for stellar catalogs.
 * 
 * Allows switching between Gaia DR3, UCAC4, and other catalogs seamlessly.
 */
class StellarCatalog {
public:
    virtual ~StellarCatalog() = default;

    virtual std::vector<Star> query_cone(const ConeQuery& params) const = 0;
    virtual std::vector<Star> query_corridor(const CorridorQuery& params) const = 0;
    virtual std::vector<Star> query_orbit(const OrbitQuery& params, ProgressCallback progress = nullptr) const = 0;

    virtual std::optional<Star> by_source_id(int64_t source_id) const = 0;
    virtual std::optional<Star> by_name(const std::string& name) const = 0;
};

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_STELLAR_CATALOG_HPP
