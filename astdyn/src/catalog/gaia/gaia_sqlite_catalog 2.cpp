#include <astdyn/catalog/gaia/gaia_sqlite_catalog.h>
#include <sqlite3.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <cmath>

namespace ioc::gaia {

GaiaSqliteCatalog::GaiaSqliteCatalog(const std::string& db_path) : db_path_(db_path) {
    openDatabase();
}

GaiaSqliteCatalog::~GaiaSqliteCatalog() {
    closeDatabase();
}

void GaiaSqliteCatalog::openDatabase() {
    int rc = sqlite3_open_v2(db_path_.c_str(), &db_, SQLITE_OPEN_READONLY, nullptr);
    if (rc != SQLITE_OK) {
        std::string error = sqlite3_errmsg(db_);
        closeDatabase();
        throw std::runtime_error("Failed to open Gaia SQLite database: " + error);
    }
}

void GaiaSqliteCatalog::closeDatabase() {
    if (db_) {
        sqlite3_close(db_);
        db_ = nullptr;
    }
}

static GaiaStar rowToStar(sqlite3_stmt* stmt) {
    GaiaStar star;
    star.source_id = sqlite3_column_int64(stmt, 0);
    star.ra = sqlite3_column_double(stmt, 1);
    star.dec = sqlite3_column_double(stmt, 2);
    star.pmra = sqlite3_column_double(stmt, 3);
    star.pmdec = sqlite3_column_double(stmt, 4);
    star.parallax = sqlite3_column_double(stmt, 5);
    star.phot_g_mean_mag = sqlite3_column_double(stmt, 6);
    star.ruwe = sqlite3_column_double(stmt, 7);
    
    const char* name = (const char*)sqlite3_column_text(stmt, 8);
    if (name) star.common_name = name;
    
    // Cross-references
    int hd = sqlite3_column_int(stmt, 12);
    if (hd > 0) star.hd_designation = "HD " + std::to_string(hd);
    
    int hip = sqlite3_column_int(stmt, 13);
    if (hip > 0) star.hip_designation = std::to_string(hip);
    
    int sao = sqlite3_column_int(stmt, 14);
    if (sao > 0) star.sao_designation = "SAO " + std::to_string(sao);
    
    return star;
}

std::vector<GaiaStar> GaiaSqliteCatalog::queryCone(double ra, double dec, double radius, double max_mag) {
    if (!db_) return {};

    std::vector<GaiaStar> results;
    
    // 1. Precise bounding box for R*Tree
    double min_ra = ra - radius / std::cos(dec * M_PI / 180.0);
    double max_ra = ra + radius / std::cos(dec * M_PI / 180.0);
    double min_dec = dec - radius;
    double max_dec = dec + radius;

    // Handle RA wrap-around if needed (simplified here for broad regions)
    // For narrow searches typical of occultations, wrap-around is rare 
    // but should be handled by two queries if cross 0/360 boundary.
    
    std::string sql;
    if (min_ra < 0 || max_ra > 360) {
        // Simple wrap-around handling
        sql = "SELECT s.* FROM stars s JOIN stars_spatial sp ON s.sid = sp.id "
              "WHERE ((sp.min_ra >= ? AND sp.max_ra <= 360) OR (sp.min_ra >= 0 AND sp.max_ra <= ?)) "
              "AND sp.min_dec >= ? AND sp.max_dec <= ? AND s.mag <= ?";
    } else {
        sql = "SELECT s.* FROM stars s JOIN stars_spatial sp ON s.sid = sp.id "
              "WHERE sp.min_ra >= ? AND sp.max_ra <= ? AND sp.min_dec >= ? AND sp.max_dec <= ? "
              "AND s.mag <= ?";
    }

    sqlite3_stmt* stmt;
    if (sqlite3_prepare_v2(db_, sql.c_str(), -1, &stmt, nullptr) != SQLITE_OK) {
        return {};
    }
    
    if (min_ra < 0 || max_ra > 360) {
        double wrapped_min = (min_ra < 0) ? min_ra + 360 : min_ra;
        double wrapped_max = (max_ra > 360) ? max_ra - 360 : max_ra;
        sqlite3_bind_double(stmt, 1, wrapped_min);
        sqlite3_bind_double(stmt, 2, wrapped_max);
        sqlite3_bind_double(stmt, 3, min_dec);
        sqlite3_bind_double(stmt, 4, max_dec);
        sqlite3_bind_double(stmt, 5, max_mag);
    } else {
        sqlite3_bind_double(stmt, 1, min_ra);
        sqlite3_bind_double(stmt, 2, max_ra);
        sqlite3_bind_double(stmt, 3, min_dec);
        sqlite3_bind_double(stmt, 4, max_dec);
        sqlite3_bind_double(stmt, 5, max_mag);
    }

    while (sqlite3_step(stmt) == SQLITE_ROW) {
        GaiaStar star = rowToStar(stmt);
        
        // Final angular distance filter
        double dist = angularDistance({ra, dec}, {star.ra, star.dec});
        if (dist <= radius) {
            results.push_back(star);
        }
    }

    sqlite3_finalize(stmt);
    return results;
}

std::optional<GaiaStar> GaiaSqliteCatalog::queryBySourceId(uint64_t source_id) {
    if (!db_) return std::nullopt;

    const char* sql = "SELECT * FROM stars WHERE sid = ?";
    sqlite3_stmt* stmt;
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return std::nullopt;
    }

    sqlite3_bind_int64(stmt, 1, source_id);

    std::optional<GaiaStar> result;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        result = rowToStar(stmt);
    }

    sqlite3_finalize(stmt);
    return result;
}

std::optional<GaiaStar> GaiaSqliteCatalog::queryByDesignation(const std::string& catalog, const std::string& designation) {
    if (!db_) return std::nullopt;

    std::string column;
    if (catalog == "SAO") column = "sao";
    else if (catalog == "HD") column = "hd";
    else if (catalog == "HIP") column = "hip";
    else if (catalog == "NAME") column = "name";
    else return std::nullopt;

    std::string sql = "SELECT * FROM stars WHERE " + column + " = ?";
    sqlite3_stmt* stmt;
    if (sqlite3_prepare_v2(db_, sql.c_str(), -1, &stmt, nullptr) != SQLITE_OK) {
        return std::nullopt;
    }

    if (catalog == "NAME") {
        sqlite3_bind_text(stmt, 1, designation.c_str(), -1, SQLITE_TRANSIENT);
    } else {
        try {
            // Remove prefix if present (e.g. "HD 1234" -> 1234)
            size_t last_space = designation.find_last_of(' ');
            std::string numeric = (last_space == std::string::npos) ? designation : designation.substr(last_space + 1);
            int id = std::stoi(numeric);
            sqlite3_bind_int(stmt, 1, id);
        } catch (...) {
            sqlite3_finalize(stmt);
            return std::nullopt;
        }
    }

    std::optional<GaiaStar> result;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        result = rowToStar(stmt);
    }

    sqlite3_finalize(stmt);
    return result;
}

size_t GaiaSqliteCatalog::getTotalStars() const {
    if (!db_) return 0;

    const char* sql = "SELECT COUNT(*) FROM stars";
    sqlite3_stmt* stmt;
    if (sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        return 0;
    }

    size_t count = 0;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        count = sqlite3_column_int64(stmt, 0);
    }

    sqlite3_finalize(stmt);
    return count;
}

} // namespace ioc::gaia
