#pragma once

#include "types.h"
#include <string>
#include <vector>
#include <optional>
#include <memory>

struct sqlite3;

namespace ioc::gaia {

/**
 * @brief Catalog implementation using an SQLite database with R*Tree indexing
 */
class GaiaSqliteCatalog {
public:
    /**
     * @brief Construct a new Gaia SQLite Catalog
     * @param db_path Path to the SQLite database file
     */
    explicit GaiaSqliteCatalog(const std::string& db_path);
    
    /**
     * @brief Destroy the Gaia SQLite Catalog
     */
    ~GaiaSqliteCatalog();

    // Disable copy
    GaiaSqliteCatalog(const GaiaSqliteCatalog&) = delete;
    GaiaSqliteCatalog& operator=(const GaiaSqliteCatalog&) = delete;

    /**
     * @brief Query stars in a cone
     */
    std::vector<GaiaStar> queryCone(double ra, double dec, double radius, double max_mag = 21.0);

    /**
     * @brief Query a star by its Gaia Source ID
     */
    std::optional<GaiaStar> queryBySourceId(uint64_t source_id);

    /**
     * @brief Query a star by various catalog designations
     * @param catalog Catalog type ("SAO", "HD", "HIP", "TYC", "NAME")
     * @param designation The identifier in that catalog
     */
    std::optional<GaiaStar> queryByDesignation(const std::string& catalog, const std::string& designation);

    /**
     * @brief Get total number of stars in the database
     */
    size_t getTotalStars() const;

    /**
     * @brief Check if the database is opened successfully
     */
    bool isOpen() const { return db_ != nullptr; }

private:
    sqlite3* db_ = nullptr;
    std::string db_path_;

    void openDatabase();
    void closeDatabase();
};

} // namespace ioc::gaia
