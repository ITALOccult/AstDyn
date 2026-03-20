/**
 * @file CatalogManager.hpp
 * @brief Factory and manager for stellar catalogs.
 */

#ifndef ASTDYN_CATALOG_MANAGER_HPP
#define ASTDYN_CATALOG_MANAGER_HPP

#include "astdyn/catalog/StellarCatalog.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <memory>
#include <string>
#include <map>

namespace astdyn::catalog {

/**
 * @brief Manages access to different stellar catalogs.
 */
class CatalogManager {
public:
    static CatalogManager& instance() {
        static CatalogManager manager;
        return manager;
    }

    /**
     * @brief Returns the preferred catalog based on a string identifier.
     */
    StellarCatalog* get_catalog(const std::string& type = "gaia_dr3") {
        if (type == "gaia_dr3") {
            try {
                return &GaiaDR3Catalog::instance();
            } catch (...) {
                return nullptr;
            }
        }
        // Placeholders for future catalogs
        return nullptr;
    }

private:
    CatalogManager() = default;
};

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_MANAGER_HPP
