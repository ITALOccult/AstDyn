#include <astdyn/AstDynEngine.hpp>
#include <astdyn/catalog/GaiaDR3Catalog.hpp>
#include <astdyn/catalog/CatalogIntegration.hpp>
#include <astdyn/core/Constants.hpp>
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    try {
        // Init catalog
        catalog::GaiaDR3Catalog::initialize("/Users/michelebigi/.ioccultcalc/catalog_config.json");
        auto& cat = catalog::GaiaDR3Catalog::instance();
        
        catalog::ConeQuery query;
        query.ra = Angle::from_deg(220.24375);
        query.dec = Angle::from_deg(14.673888);
        query.radius = Angle::from_arcsec(60.0);
        query.max_magnitude = 15.5;
        
        auto stars = cat.query_cone(query);
        std::cout << "Found " << stars.size() << " stars near Haumea event target.\n";
        for (const auto& s : stars) {
            std::cout << "Source ID: " << s.source_id 
                      << " RA: " << s.ra.to_deg() 
                      << " Dec: " << s.dec.to_deg() 
                      << " G_mag: " << s.g_mag << "\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
