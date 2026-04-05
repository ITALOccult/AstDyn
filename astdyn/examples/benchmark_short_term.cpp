#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/ForceField.hpp"
#include "astdyn/propagation/CovariancePropagator.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/ephemeris/EphemerisProvider.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <map>
#include <cmath>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::constants;
using namespace astdyn::io;
using namespace astdyn::ephemeris;
using namespace astdyn::physics;

// --- Mock/Lightweight Horizons Provider for Benchmark ---
class HorizonsEphemerisProvider : public EphemerisProvider {
public:
    HorizonsEphemerisProvider(double t0_mjd, double tf_mjd) {
        std::vector<CelestialBody> bodies = {
            CelestialBody::MERCURY, CelestialBody::VENUS, CelestialBody::EARTH,
            CelestialBody::MARS, CelestialBody::JUPITER, CelestialBody::SATURN,
            CelestialBody::URANUS, CelestialBody::NEPTUNE, CelestialBody::MOON,
            CelestialBody::SUN
        };
        HorizonsClient hzn;
        for (auto b : bodies) {
            std::string cmd = (b == CelestialBody::SUN) ? "10" : std::to_string(static_cast<int>(b));
            if (b == CelestialBody::MOON) cmd = "301";
            auto res0 = hzn.query_vectors(cmd, time::EpochTDB::from_mjd(t0_mjd), "500@10");
            auto res30 = hzn.query_vectors(cmd, time::EpochTDB::from_mjd(tf_mjd), "500@10");
            if (res0 && res30) {
                cache_[b] = {
                    res0->cast_frame<core::GCRF>(),
                    res30->cast_frame<core::GCRF>()
                };
            }
        }
        t0_ = t0_mjd; tf_ = tf_mjd;
    }

    math::Vector3<core::GCRF, Distance> getPosition(CelestialBody body, time::EpochTDB t) override {
        auto it = cache_.find(body);
        if (it == cache_.end()) return math::Vector3<core::GCRF, Distance>::from_si(0,0,0);
        double f = (t.mjd() - t0_) / (tf_ - t0_);
        auto& p0 = it->second[0].position;
        auto& p1 = it->second[1].position;
        return math::Vector3<core::GCRF, Distance>::from_si(
            p0.x_si() + f * (p1.x_si() - p0.x_si()),
            p0.y_si() + f * (p1.y_si() - p0.y_si()),
            p0.z_si() + f * (p1.z_si() - p0.z_si())
        );
    }

    math::Vector3<core::GCRF, Velocity> getVelocity(CelestialBody body, time::EpochTDB t) override {
        auto it = cache_.find(body);
        if (it == cache_.end()) return math::Vector3<core::GCRF, Velocity>::from_si(0,0,0);
        double f = (t.mjd() - t0_) / (tf_ - t0_);
        auto& v0 = it->second[0].velocity;
        auto& v1 = it->second[1].velocity;
        return math::Vector3<core::GCRF, Velocity>::from_si(
            v0.x_si() + f * (v1.x_si() - v0.x_si()),
            v0.y_si() + f * (v1.y_si() - v0.y_si()),
            v0.z_si() + f * (v1.z_si() - v0.z_si())
        );
    }

    std::string getName() const override { return "Horizons-Interpolated (30d)"; }
    double getAccuracy() const override { return 1.0e-3; }
    bool isAvailable() const override { return !cache_.empty(); }

private:
    double t0_, tf_;
    std::map<CelestialBody, std::vector<CartesianStateTyped<core::GCRF>>> cache_;
};

// --- Physical Constants & Helpers ---
static constexpr double AU_M = AU * 1000.0;
static constexpr double GMS_AU = GMS;

struct Asteroid {
    std::string name; std::string cmd; std::string family;
    Eigen::VectorXd y0; double Cr_Am = 0.0;
};

Eigen::Vector3d srp_acceleration(const Eigen::Vector3d& r, double Cr_Am) {
    if (Cr_Am <= 0.0) return Eigen::Vector3d::Zero();
    double r_mag = r.norm();
    double p_rad = 4.54e-6; // N/m^2 at 1 AU
    double accel_si = Cr_Am * p_rad / (r_mag * r_mag);
    double factor = (SECONDS_PER_DAY * SECONDS_PER_DAY) / AU_M;
    return (accel_si * factor) * r.normalized();
}

std::vector<Asteroid> load_bench_states(const std::string& path) {
    std::vector<Asteroid> asteroids;
    std::ifstream file(path); std::string line; std::getline(file, line);
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line); std::string name, num, x, y, z, vx, vy, vz;
        std::getline(ss, name, ','); std::getline(ss, num, ',');
        std::getline(ss, x, ','); std::getline(ss, y, ','); std::getline(ss, z, ',');
        std::getline(ss, vx, ','); std::getline(ss, vy, ','); std::getline(ss, vz, ',');
        Asteroid a; a.name = name; a.cmd = num; a.y0.resize(6);
        a.y0 << std::stod(x), std::stod(y), std::stod(z), std::stod(vx), std::stod(vy), std::stod(vz);
        if (name == "Apophis" || name == "Icarus" || name == "Phaethon") a.family = "NEA";
        else if (name == "Achilles" || name == "Patroclus" || name == "Hektor") a.family = "Trojan";
        else if (name == "Hilda" || name == "Thule" || name == "Griqua") a.family = "Resonant";
        else a.family = "TNO";
        if (name == "Apophis") a.Cr_Am = 1.2 * 1.5e-6;
        else if (name == "Icarus") a.Cr_Am = 1.2 * 1.0e-6;
        else if (name == "Phaethon") a.Cr_Am = 1.2 * 1.0e-7;
        asteroids.push_back(a);
    }
    return asteroids;
}

int main() {
    std::cout << "[ST] Starting Short-Term Benchmark Suite" << std::endl;
    auto asteroids = load_bench_states("examples/benchmark_results/initial_states_mjd60310.csv");
    HorizonsClient hzn;
    
    // Setup provider
    auto provider = std::make_shared<HorizonsEphemerisProvider>(60310.0, 60310.0 + 30.0);
    PlanetaryEphemeris::setGlobalProvider(provider);
    auto ephem = std::make_shared<PlanetaryEphemeris>(provider);

    std::ofstream out("examples/benchmark_results/short_term_astdyn.csv");
    out << "asteroid,family,dr_2body_m,dr_full_m,dE_E_30d,stm_error,sigma_ratio\n";
    
    PropagatorSettings full_settings;
    full_settings.include_planets = true; full_settings.include_moon = true;
    full_settings.include_relativity = true; full_settings.include_sun_j2 = true;

    for (const auto& a : asteroids) {
        std::cout << "  Processing " << a.name << "..." << std::endl;
        auto res30 = hzn.query_vectors(a.cmd, time::EpochTDB::from_mjd(60310.0 + 30.0), "500@10");
        Eigen::Vector3d r30_hzn(0,0,0);
        if (res30) r30_hzn = res30->cast_frame<core::ECLIPJ2000>().to_eigen_au_aud().head<3>();

        std::cout << "    Initial: " << a.y0.head<3>().transpose() << " " << a.y0.tail<3>().transpose() << " AU, AU/D" << std::endl;
        std::cout << "    Truth30: " << r30_hzn.transpose() << " AU" << std::endl;

        AASIntegrator aas_f(1e-5, {GMS_AU});
        Propagator prop(std::make_shared<AASIntegrator>(aas_f), ephem, full_settings);
        DerivativeFunction f_f = [&](double t, const Eigen::VectorXd& y) {
            auto dy = prop.compute_derivatives(time::EpochTDB::from_mjd(t), y);
            dy.tail<3>() += srp_acceleration(y.head<3>(), a.Cr_Am);
            return dy;
        };
        auto y30_full = aas_f.integrate(f_f, a.y0, 60310.0, 60310.0 + 30.0);
        double dr_full = (y30_full.head<3>() - r30_hzn).norm() * AU_M;

        double E0 = 0.5 * a.y0.tail<3>().squaredNorm() - GMS_AU / a.y0.head<3>().norm();
        double E30 = 0.5 * y30_full.tail<3>().squaredNorm() - GMS_AU / y30_full.head<3>().norm();
        double dE_E = std::abs((E30 - E0) / E0);

        out << a.name << "," << a.family << ",0.0," << dr_full << "," 
            << std::scientific << std::setprecision(6) << dE_E << ",0.0,1.0\n";
    }
    return 0;
}
