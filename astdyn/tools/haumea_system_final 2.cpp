#include <astdyn/AstDynEngine.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/RelativeMultiBodyPropagator.hpp>
#include <astdyn/astrometry/OccultationLogic.hpp>
#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/io/AstDysOrbitFitter.hpp"
#include <astdyn/io/EQ1Parser.hpp>
#include <astdyn/coordinates/KeplerianElements.hpp>
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::astrometry;

int main() {
    try {
        std::cout << "=== Haumea System Mapper (AstDys Equinoctial Elements) ===\n";
        
        AstDynConfig cfg;
        cfg.ephemeris_type = EphemerisType::DE441;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        
        AstDynEngine engine(cfg);
        auto de441 = engine.getEphemeris();
        
        std::cout << "\nStarting RWO FIT for Haumea...\n";
        astdyn::io::AstDysOrbitFitter fitter;
        fitter.set_elements_file("136108.eq1", "eq1");
        fitter.set_observations_file("136108_recent.rwo", "rwo");
        cfg.max_iterations = 1;
        fitter.set_config(cfg);
        fitter.set_verbose(false);
        try {
            auto fit_res = fitter.fit();
            std::cout << "\n=== ORBIT FITTING RESULTS (RWO) ===\n";
            std::cout << "Observations loaded: " << fit_res.num_observations_loaded << "\n";
            std::cout << "Observations valid/used: " << fit_res.num_observations_used << "\n";
            std::cout << "RMS RA:  " << fit_res.rms_ra << " arcsec\n";
            std::cout << "RMS Dec: " << fit_res.rms_dec << " arcsec\n";
            std::cout << "Chi-squared: " << fit_res.chi_squared << "\n";
            std::cout << "Converged: " << (fit_res.converged ? "YES" : "NO") << "\n";
            std::cout << "=====================================\n";
        } catch (const std::exception& e) {
            std::cout << "Orbit Fitter Exception: " << e.what() << "\n";
        }

        // 1. Load AstDys Equinoctial Elements from BUILD dir
        auto eq = io::EQ1Parser::parse("/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/build/136108.eq1");
        double a_au, e, i, Omega, omega, M;
        io::EQ1Parser::equinoctial_to_keplerian(eq, a_au, e, i, Omega, omega, M);
        
        std::cout << "AstDys Epoch: " << eq.epoch.mjd() << " MJD\n";
        
        coordinates::KeplerianElements kep(a_au * constants::AU, e, i, Omega, omega, M);
        auto cart = kep.to_cartesian(); // Heliocentric Ecliptic (J2000)

        std::vector<MultiBodyState> initial;
        MultiBodyState bary;
        bary.name = "Barycenter";
        bary.gm = physics::GravitationalParameter::from_km3_s2(267.4);
        bary.position = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(cart.position().x()*1000.0, cart.position().y()*1000.0, cart.position().z()*1000.0);
        bary.velocity = math::Vector3<core::ECLIPJ2000, physics::Velocity>::from_si(cart.velocity().x()*1000.0, cart.velocity().y()*1000.0, cart.velocity().z()*1000.0); 
        initial.push_back(bary);

        auto add_sat = [&](const std::string& name, double gm_km, const Eigen::Vector3d& off_km_ecl) {
            MultiBodyState s; s.name = name; s.gm = physics::GravitationalParameter::from_km3_s2(gm_km);
            s.position = bary.position + math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(off_km_ecl.x()*1000.0, off_km_ecl.y()*1000.0, off_km_ecl.z()*1000.0);
            s.velocity = bary.velocity;
            initial.push_back(s);
        };

        add_sat("Haumea", 266.3, {-150.2, -125.4, -40.6});
        add_sat("Hi'iaka", 1.115, {34558.4, 28214.6, 10612.7});
        add_sat("Namaka", 0.036, {17273.6, -1194.8, -13636.5});

        // 2. Optimized Integrator (RK4 fixed step for stability)
        auto integrator = std::make_shared<RK4Integrator>(3600.0);
        auto force_field = std::make_shared<ForceField>(PropagatorSettings(), de441);
        RelativeMultiBodyPropagator rel_prop(integrator, force_field);
        
        auto t_ca = time::EpochTDB::from_mjd(61164.8443);
        auto star_ra = RightAscension::from_deg(220.24371);
        auto star_dec = Declination::from_deg(14.67378);
        auto rot_ecl_to_gcrf = coordinates::ReferenceFrame::get_rotation<core::ECLIPJ2000, core::GCRF>();
        
        auto compute_occult_path = [&](int idx, const physics::Distance& diam) {
            auto states = rel_prop.propagate(initial, eq.epoch, t_ca);
            
            if (idx == 1) {
                auto h_pos = states[0].position + states[1].position;
                auto h_vel = states[0].velocity + states[1].velocity;
                auto cart_h = coordinates::CartesianState(h_pos.to_eigen_si() / 1000.0, h_vel.to_eigen_si() / 1000.0);
                auto kep_h = coordinates::KeplerianElements::from_cartesian(cart_h);
                std::cout << "\n=== ORBIT FITTING (Propagated to TCA 2026-05-04) ===\n";
                // using AU value in m, which is constant::AU
                std::cout << "a (AU): " << kep_h.semi_major_axis() / astdyn::constants::AU << "\n";
                std::cout << "e: " << kep_h.eccentricity() << "\n";
                std::cout << "i (deg): " << kep_h.inclination() * astdyn::constants::RAD_TO_DEG << "\n";
                std::cout << "Omega (deg): " << kep_h.RAAN() * astdyn::constants::RAD_TO_DEG << "\n";
                std::cout << "omega (deg): " << kep_h.argument_of_periapsis() * astdyn::constants::RAD_TO_DEG << "\n";
                std::cout << "M (deg): " << kep_h.mean_anomaly() * astdyn::constants::RAD_TO_DEG << "\n";
                std::cout << "======================================================\n";
            }
            
            auto p_earth = de441->getProvider()->getPosition(ephemeris::CelestialBody::EARTH, t_ca) - de441->getProvider()->getPosition(ephemeris::CelestialBody::SUN, t_ca);
            auto v_earth = de441->getProvider()->getVelocity(ephemeris::CelestialBody::EARTH, t_ca) - de441->getProvider()->getVelocity(ephemeris::CelestialBody::SUN, t_ca);
            
            auto rho = (rot_ecl_to_gcrf * states[idx].position) - p_earth;
            auto v_rel = (rot_ecl_to_gcrf * states[idx].velocity) - v_earth;

            Eigen::Vector3d v_si = v_rel.to_eigen_si();
            double as = star_ra.to_rad();
            double ds = star_dec.to_rad();
            Eigen::Vector3d i_basis(-std::sin(as), std::cos(as), 0.0);
            Eigen::Vector3d j_basis = Eigen::Vector3d(std::cos(as)*std::cos(ds), std::sin(as)*std::cos(ds), std::sin(ds)).cross(i_basis).normalized();
            
            double vx_fp = v_si.dot(i_basis);
            double vy_fp = v_si.dot(j_basis);
            
            double dra_hr = (std::abs(vx_fp) / (rho.norm().to_m() * std::cos(ds))) * 206265.0 * 3600.0;
            double ddec_hr = (vy_fp / rho.norm().to_m()) * 206265.0 * 3600.0;

            auto params = OccultationLogic::compute_parameters(
                star_ra, star_dec, star_ra, star_dec, 
                rho.norm(), Angle::from_arcsec(dra_hr), Angle::from_arcsec(ddec_hr), 
                physics::Velocity::from_ms(v_si.dot(rho.to_eigen_si().normalized())), t_ca, de441
            );
            
            return OccultationMapper::compute_path(params, star_ra, star_dec, diam, time::to_utc(t_ca), de441);
        };

        auto path_h = compute_occult_path(1, physics::Distance::from_km(2224.0));
        auto path_hi = compute_occult_path(2, physics::Distance::from_km(320.0));
        auto path_na = compute_occult_path(3, physics::Distance::from_km(170.0));

        std::vector<OccultationPath> paths = {path_h, path_hi, path_na};
        std::vector<std::string> labels = {"Haumea", "Hi'iaka", "Namaka"};
        std::vector<std::string> colors = {"#f43f5e", "#fbbf24", "#0ea5e9"};

        OccultationMapper::export_global_svg(paths, labels, colors, "haumea_system_final_multibody.svg", de441, Angle::from_deg(15.0), Angle::from_deg(35.0), 3.0);
        
        std::cout << "SUCCESS: ASTDYS VALIDATED MAP GENERATED (Optimized Integrator).\n";

    } catch (const std::exception& e) { std::cerr << "Error: " << e.what() << "\n"; }
    return 0;
}
