#pragma once
#include "constants.h"
#include "linalg.h"
#include "orbital_elements.h"
#include "integrator.h"
#include "differential_corrections.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace orbfit {

// ─────────────────────────────────────────────────────────────
//  Confidence ellipsoid in element space
//  Defined by covariance matrix C:  Δx^T C^{-1} Δx ≤ χ²_{6,p}
// ─────────────────────────────────────────────────────────────
struct ConfidenceEllipsoid {
    Matrix covariance;   // 6×6
    double sigma_level;  // e.g. 3.0 for 3σ

    // Sample k virtual objects along the Line Of Variations (LOV)
    // LOV = first eigenvector of C (direction of max uncertainty)
    std::vector<OrbitalElements> lov_sample(const OrbitalElements& nominal,
                                             int k = 101) const
    {
        // Eigenvalue decomposition via power iteration for dominant eigenvector
        // For small 6×6 we use Jacobi method
        int n = 6;
        Matrix V = Matrix::identity(n);
        Matrix A = covariance;  // working copy

        // Jacobi iterations
        for(int sweep=0; sweep<100; sweep++){
            double off = 0;
            for(int i=0;i<n;i++)
                for(int j=i+1;j<n;j++) off+=A(i,j)*A(i,j);
            if(std::sqrt(off)<1e-15) break;

            for(int p=0;p<n;p++){
                for(int q=p+1;q<n;q++){
                    double Apq=A(p,q);
                    if(std::abs(Apq)<1e-14) continue;
                    double App=A(p,p), Aqq=A(q,q);
                    double theta=(Aqq-App)/(2*Apq);
                    double t=std::copysign(1.0,theta)
                             /(std::abs(theta)+std::sqrt(theta*theta+1));
                    double c=1/std::sqrt(t*t+1), s=t*c;
                    // Rotate A
                    A(p,p)=App-t*Apq; A(q,q)=Aqq+t*Apq; A(p,q)=A(q,p)=0;
                    for(int r=0;r<n;r++){
                        if(r==p||r==q) continue;
                        double Arp=A(r,p), Arq=A(r,q);
                        A(r,p)=A(p,r)=c*Arp-s*Arq;
                        A(r,q)=A(q,r)=s*Arp+c*Arq;
                    }
                    // Update eigenvectors
                    for(int r=0;r<n;r++){
                        double Vrp=V(r,p), Vrq=V(r,q);
                        V(r,p)=c*Vrp-s*Vrq;
                        V(r,q)=s*Vrp+c*Vrq;
                    }
                }
            }
        }

        // Find max eigenvalue index
        int imax=0;
        double lmax=A(0,0);
        for(int i=1;i<n;i++) if(A(i,i)>lmax){lmax=A(i,i);imax=i;}

        // LOV direction = eigenvector corresponding to lmax
        std::vector<double> lov(n);
        for(int i=0;i<n;i++) lov[i]=V(i,imax);

        double lov_scale=sigma_level*std::sqrt(lmax);

        // Generate k points from -1 to +1 along LOV
        std::vector<OrbitalElements> pts;
        pts.reserve(k);
        for(int j=0;j<k;j++){
            double s = lov_scale*(-1.0 + 2.0*j/(k-1));
            OrbitalElements el = nominal;
            el.a     += s*lov[0];
            el.e     += s*lov[1];
            el.i     += s*lov[2];
            el.Omega += s*lov[3];
            el.omega += s*lov[4];
            el.M0    += s*lov[5];
            if(el.a>0 && el.e>=0 && el.e<1.0) pts.push_back(el);
        }
        return pts;
    }
};

// ─────────────────────────────────────────────────────────────
//  Minimum orbital intersection distance (MOID)
//  Between target orbit and Earth orbit
//  Uses grid search + Brent refinement
// ─────────────────────────────────────────────────────────────
inline double moid_estimate(const OrbitalElements& el,
                              double t_ref = JD_J2000,
                              int grid = 360)
{
    // Earth orbital elements (approximate)
    OrbitalElements earth;
    earth.a=1.0; earth.e=0.0167; earth.i=0.0; earth.Omega=0.0;
    earth.omega=1.7966; earth.M0=0.0; earth.epoch=t_ref;

    double moid = 1e10;
    // Grid search over mean anomalies
    for(int j=0;j<grid;j++){
        OrbitalElements el2=el;
        el2.M0 = 2.0*M_PI*j/grid;
        auto sv1=elements_to_state(el2, t_ref);

        for(int k=0;k<grid;k++){
            OrbitalElements e2=earth;
            e2.M0=2.0*M_PI*k/grid;
            auto sv2=elements_to_state(e2, t_ref);
            double d=(sv1.pos-sv2.pos).norm();
            if(d<moid) moid=d;
        }
    }
    return moid;
}

// ─────────────────────────────────────────────────────────────
//  Close approach / impact detection
// ─────────────────────────────────────────────────────────────
struct CloseApproach {
    double t_ca;          // [JD] time of closest approach
    double distance;      // [AU]
    double v_inf;         // [km/s] hyperbolic excess speed
    bool   possible_impact; // distance < Earth radius
};

inline std::vector<CloseApproach>
find_close_approaches(const OrbitalElements& el,
                       double t_start, double t_end,
                       double planet_gm   = GM_EARTH + GM_MOON,
                       int    planet_id   = 2,      // Earth
                       double alert_dist  = 0.05,   // AU
                       const IntegratorOptions& iopt = {})
{
    std::vector<CloseApproach> results;
    double dt_coarse = 1.0;   // [day] coarse scan step

    StateVector sv = elements_to_state(el, el.epoch);

    double d_prev = 1e10;
    double t = t_start;

    while(t <= t_end){
        auto sv_t  = propagate(sv, t, iopt);
        auto pl    = planet_state(planet_id, t);
        double d   = (sv_t.pos - pl.pos).norm();

        if(d < alert_dist){
            if(d < d_prev){
                d_prev = d; t += dt_coarse; continue;
            }
            // Local minimum found near t-dt_coarse
            double t_ca = t - dt_coarse;
            // Refine with bisection
            double ta=t-2*dt_coarse, tb=t;
            for(int ref=0;ref<50;ref++){
                double tm=(ta+tb)/2;
                auto sa=propagate(sv,ta,iopt), sb=propagate(sv,tb,iopt);
                auto sm=propagate(sv,tm,iopt);
                auto pla=planet_state(planet_id,ta);
                auto plb=planet_state(planet_id,tb);
                auto plm=planet_state(planet_id,tm);
                double da=(sa.pos-pla.pos).norm();
                double db=(sb.pos-plb.pos).norm();
                double dm=(sm.pos-plm.pos).norm();
                if(dm<da && dm<db){ta=ta;tb=tb;}
                else if(da<db){tb=tm;}
                else {ta=tm;}
                if(tb-ta<1e-5) break;
            }
            t_ca=(ta+tb)/2;
            auto sv_ca = propagate(sv, t_ca, iopt);
            auto pl_ca = planet_state(planet_id, t_ca);
            double d_ca= (sv_ca.pos - pl_ca.pos).norm();
            Vec3 dv    = sv_ca.vel - pl_ca.vel;
            double vinf= dv.norm() * AU_KM / 86400.0;  // km/s

            CloseApproach ca;
            ca.t_ca   = t_ca;
            ca.distance = d_ca;
            ca.v_inf  = vinf;
            ca.possible_impact = (d_ca < EARTH_RADIUS * 2.0);
            results.push_back(ca);
        }
        d_prev = d;
        t += dt_coarse;
    }
    return results;
}

// ─────────────────────────────────────────────────────────────
//  Palermo Scale (impact risk metric)
//  PS = log10( p_impact / (f_background * T_warning) )
//  where f_background = 0.03 * E^{-0.8} [impacts/yr] for energy E [MT]
// ─────────────────────────────────────────────────────────────
inline double palermo_scale(double p_impact, double t_warning_yr,
                              double diameter_km = 0.1)
{
    // Kinetic energy ~ (d/0.1km)^3 * 10 MT
    double energy_MT = std::pow(diameter_km/0.1, 3.0) * 10.0;
    double f_bg      = 0.03 * std::pow(energy_MT, -0.8);
    double p_bg      = f_bg * t_warning_yr;
    if(p_bg <= 0 || p_impact <= 0) return -1e10;
    return std::log10(p_impact / p_bg);
}

// ─────────────────────────────────────────────────────────────
//  Covariance propagation via STM
//  C_out = Φ C_in Φ^T
// ─────────────────────────────────────────────────────────────
inline Matrix propagate_covariance(const Matrix& C_in,
                                    const OrbitalElements& el,
                                    double t_prop,
                                    const IntegratorOptions& iopt = {})
{
    auto sv0 = elements_to_state(el, el.epoch);
    auto stm = propagate_stm(sv0, t_prop, iopt);

    // Build Phi matrix
    Matrix Phi(6,6);
    for(int i=0;i<6;i++)
        for(int j=0;j<6;j++)
            Phi(i,j)=stm.Phi[i*6+j];

    // Jacobian ∂(r,v)/∂(a,e,i,Ω,ω,M0) at t0
    Matrix J = elements_jacobian(el, el.epoch);
    Matrix Jinv = invert(J);

    // Transform C from element space to Cartesian
    Matrix C_cart = J * C_in * J.transpose();

    // Propagate
    Matrix C_cart_out = Phi * C_cart * Phi.transpose();

    // Transform back
    // J_out at t_prop
    StateVector sv_out; sv_out.epoch=t_prop;
    sv_out.pos=Vec3(stm.y[0],stm.y[1],stm.y[2]);
    sv_out.vel=Vec3(stm.y[3],stm.y[4],stm.y[5]);
    OrbitalElements el_out = state_to_elements(sv_out);
    Matrix Jout = elements_jacobian(el_out, t_prop);
    Matrix Jout_inv = invert(Jout);

    return Jout_inv * C_cart_out * Jout_inv.transpose();
}

} // namespace orbfit
