#pragma once
#include "constants.h"
#include "linalg.h"
#include "orbital_elements.h"
#include "integrator.h"
#include <vector>
#include <string>
#include <cmath>

namespace orbfit {

// ─────────────────────────────────────────────────────────────
//  Astrometric observation (RA/Dec)
// ─────────────────────────────────────────────────────────────
struct Observation {
    double t_obs;     // [JD, UTC approx.]
    double ra;        // right ascension [rad]
    double dec;       // declination     [rad]
    double sigma_ra;  // uncertainty in RA*cos(dec) [rad]
    double sigma_dec; // uncertainty in Dec          [rad]
    Vec3   obs_pos;   // observer heliocentric position [AU]
    std::string station_code;

    // Convenience: weight = 1/sigma^2
    double w_ra()  const { return 1.0/(sigma_ra *sigma_ra ); }
    double w_dec() const { return 1.0/(sigma_dec*sigma_dec); }
};

// ─────────────────────────────────────────────────────────────
//  Observation model: given state vector, compute predicted
//  (RA, Dec) with light-time correction
// ─────────────────────────────────────────────────────────────
struct PredictedObs {
    double ra, dec;   // [rad]
    Vec3   rho;       // topocentric vector (target - obs) [AU]
    double delta;     // distance [AU]
    double lt;        // light-travel time [day]
};

// Compute predicted RA/Dec with iterative light-time correction
// sv: state at some epoch; obs: observation record
// integrator propagates sv to t_obs - lt
inline PredictedObs predict(const OrbitalElements& el,
                             const Observation& obs,
                             const IntegratorOptions& iopt = {})
{
    // First approximation: no light-time
    auto sv = propagate(elements_to_state(el, el.epoch), obs.t_obs, iopt);
    Vec3 rho = sv.pos - obs.obs_pos;
    double delta = rho.norm();
    double lt    = delta / CLIGHT;

    // Iterate (2 iterations usually sufficient)
    for(int iter=0; iter<3; iter++){
        sv = propagate(elements_to_state(el, el.epoch), obs.t_obs - lt, iopt);
        rho   = sv.pos - obs.obs_pos;
        delta = rho.norm();
        lt    = delta / CLIGHT;
    }

    // Compute RA/Dec
    double x=rho.x, y=rho.y, z=rho.z;
    // Equatorial coords: rotate ecliptic → equatorial (ε = 23.4392911 deg)
    const double eps = 23.4392911 * DEG2RAD;
    double ce=std::cos(eps), se=std::sin(eps);
    double xeq = x;
    double yeq = ce*y - se*z;
    double zeq = se*y + ce*z;

    double ra  = std::atan2(yeq, xeq);
    if(ra < 0) ra += 2.0*M_PI;
    double dec = std::asin(zeq / std::sqrt(xeq*xeq+yeq*yeq+zeq*zeq));

    return {ra, dec, rho, delta, lt};
}

// ─────────────────────────────────────────────────────────────
//  Residuals for one observation  (O-C)
//  Returns [dRA*cos(Dec), dDec] in arcsec
// ─────────────────────────────────────────────────────────────
inline std::pair<double,double> residual(const OrbitalElements& el,
                                          const Observation& obs,
                                          const IntegratorOptions& iopt = {})
{
    auto pred = predict(el, obs, iopt);
    // Use computed (predicted) declination for the RA*cos(dec) projection so that
    // residuals and their partials share the same reference frame.
    // This matches the convention in astdyn::orbit_determination::Residuals.
    double dra  = (obs.ra  - pred.ra ) * std::cos(pred.dec);
    double ddec =  obs.dec - pred.dec;

    // Wrap dra to [-pi, pi]
    if(dra >  M_PI) dra -= 2*M_PI;
    if(dra < -M_PI) dra += 2*M_PI;

    return { dra * RAD2ARCSEC, ddec * RAD2ARCSEC };
}

// ─────────────────────────────────────────────────────────────
//  Partials ∂(RA*cos(Dec), Dec)/∂(a,e,i,Ω,ω,M0)
//  via numerical differentiation of predict()
// ─────────────────────────────────────────────────────────────
inline Matrix obs_partials(const OrbitalElements& el,
                            const Observation& obs,
                            const IntegratorOptions& iopt = {})
{
    const double EPS = 1e-6;
    Matrix H(2, 6);

    auto pack = [](const OrbitalElements& e) -> std::array<double,6> {
        return {e.a, e.e, e.i, e.Omega, e.omega, e.M0};
    };
    auto unpack = [](std::array<double,6> p, OrbitalElements base) -> OrbitalElements {
        base.a=p[0]; base.e=p[1]; base.i=p[2];
        base.Omega=p[3]; base.omega=p[4]; base.M0=p[5];
        return base;
    };

    auto p0 = pack(el);
    // Use predicted (computed) declination as projection reference — consistent with residual().
    double cosd = std::cos(predict(el, obs, iopt).dec);

    for(int j=0; j<6; j++){
        auto pp = p0; pp[j] += EPS*std::max(std::abs(p0[j]),1e-10);
        auto pm = p0; pm[j] -= EPS*std::max(std::abs(p0[j]),1e-10);
        double h2 = pp[j] - pm[j];

        auto predp = predict(unpack(pp,el), obs, iopt);
        auto predm = predict(unpack(pm,el), obs, iopt);

        double dra_p  = predp.ra * cosd;
        double dra_m  = predm.ra * cosd;
        double ddec_p = predp.dec;
        double ddec_m = predm.dec;

        H(0, j) = (dra_p  - dra_m ) / h2 * RAD2ARCSEC;
        H(1, j) = (ddec_p - ddec_m) / h2 * RAD2ARCSEC;
    }
    return H;
}

// ─────────────────────────────────────────────────────────────
//  Compute RMS of residuals  [arcsec]
// ─────────────────────────────────────────────────────────────
inline double compute_rms(const OrbitalElements& el,
                           const std::vector<Observation>& obs,
                           const IntegratorOptions& iopt = {})
{
    double sum = 0.0;
    int n = 0;
    for(auto& o : obs){
        auto [dra, ddec] = residual(el, o, iopt);
        sum += dra*dra + ddec*ddec;
        n   += 2;
    }
    return std::sqrt(sum / n);
}

} // namespace orbfit
