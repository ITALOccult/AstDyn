/**
 * @file CatalogManager.hpp
 * @brief Factory and manager for stellar catalogs.
 */

#ifndef ASTDYN_CATALOG_MANAGER_HPP
#define ASTDYN_CATALOG_MANAGER_HPP

#include "astdyn/catalog/StellarCatalog.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/catalog/UCAC4Catalog.hpp"
#include "astdyn/catalog/Tycho2Catalog.hpp"
#include "astdyn/catalog/HipparcosCatalog.hpp"
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
        if (type == "ucac4") {
            try { return &UCAC4Catalog::instance(); } catch (...) { return nullptr; }
        }
        if (type == "tycho2") {
            try { return &Tycho2Catalog::instance(); } catch (...) { return nullptr; }
        }
        if (type == "hipparcos") {
            try { return &HipparcosCatalog::instance(); } catch (...) { return nullptr; }
        }
        return nullptr;
    }

private:
    CatalogManager() = default;
};

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_MANAGER_HPP
