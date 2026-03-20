#pragma once
#include "constants.h"
#include "linalg.h"
#include "orbital_elements.h"
#include "observations.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace orbfit {

// ─────────────────────────────────────────────────────────────
//  Fit result
// ─────────────────────────────────────────────────────────────
struct FitResult {
    OrbitalElements elements;
    Matrix          covariance;   // 6×6 covariance in element space
    Matrix          normal;       // 6×6 normal matrix (C^-1)
    double          rms;          // weighted RMS [arcsec]
    double          chi2;         // χ²
    int             n_obs;        // number of observations used
    int             iterations;   // DC iterations used
    bool            converged;

    // Sigma (1-sigma uncertainty) for each element
    std::array<double,6> sigma() const {
        std::array<double,6> s;
        for(int i=0;i<6;i++) s[i]=std::sqrt(covariance(i,i));
        return s;
    }

    void print() const {
        auto s = sigma();
        std::cout << std::fixed;
        std::cout << "=== OrbFit Result ===\n";
        std::cout << std::setw(8) << "a"     << " = " << std::setprecision(8) << elements.a
                  << " ± " << s[0] << " AU\n";
        std::cout << std::setw(8) << "e"     << " = " << std::setprecision(8) << elements.e
                  << " ± " << s[1] << "\n";
        std::cout << std::setw(8) << "i"     << " = " << std::setprecision(6) << elements.i*RAD2DEG
                  << " ± " << s[2]*RAD2DEG << " deg\n";
        std::cout << std::setw(8) << "Omega" << " = " << std::setprecision(6) << elements.Omega*RAD2DEG
                  << " ± " << s[3]*RAD2DEG << " deg\n";
        std::cout << std::setw(8) << "omega" << " = " << std::setprecision(6) << elements.omega*RAD2DEG
                  << " ± " << s[4]*RAD2DEG << " deg\n";
        std::cout << std::setw(8) << "M0"    << " = " << std::setprecision(6) << elements.M0*RAD2DEG
                  << " ± " << s[5]*RAD2DEG << " deg\n";
        std::cout << std::setw(8) << "epoch" << " = JD " << std::setprecision(3) << elements.epoch << "\n";
        std::cout << "RMS = " << std::setprecision(4) << rms << " arcsec,  "
                  << "chi2 = " << chi2 << ",  nobs = " << n_obs
                  << ",  iter = " << iterations << "\n";
        std::cout << "Converged: " << (converged?"YES":"NO") << "\n";
    }
};

// ─────────────────────────────────────────────────────────────
//  Differential Corrections options
// ─────────────────────────────────────────────────────────────
struct DCOptions {
    int    max_iter     = 50;
    double conv_rms     = 1e-5;   // convergence: RMS change [arcsec]
    double conv_delta   = 1e-10;  // convergence: element change norm
    double lambda_init  = 1e-3;   // Levenberg-Marquardt initial damping
    bool   use_lm       = true;   // Levenberg-Marquardt (robust) vs pure GN
    bool   verbose      = false;
    IntegratorOptions iopt;
};

// ─────────────────────────────────────────────────────────────
//  Differential Corrections (Gauss-Newton + Levenberg-Marquardt)
//  Minimises  χ² = Σ [(dRA/σ_RA)² + (dDec/σ_Dec)²]
//
//  Elements parameterised as (a, e, i, Ω, ω, M0)
//  Normal equations:  (H^T W H) Δx = H^T W ξ
// ─────────────────────────────────────────────────────────────
inline FitResult differential_corrections(OrbitalElements el0,
                                           const std::vector<Observation>& obs,
                                           const DCOptions& dco = {})
{
    int nobs = obs.size();
    if(nobs < 3)
        throw std::runtime_error("Need at least 3 observations");

    int m = 2*nobs;  // number of residuals
    int n = 6;       // number of parameters

    OrbitalElements el = el0;
    double rms_prev = 1e30;
    double lambda   = dco.lambda_init;
    bool converged  = false;
    int  iter;

    FitResult res;
    res.n_obs = nobs;
    res.covariance = Matrix(6,6,0.0);
    res.normal     = Matrix(6,6,0.0);

    for(iter = 0; iter < dco.max_iter; iter++){
        // Build design matrix H (m×6) and residual vector xi (m×1)
        Matrix H(m, n, 0.0);
        std::vector<double> xi(m, 0.0);
        std::vector<double> w(m, 0.0);   // weights

        double chi2 = 0.0;

        for(int k=0; k<nobs; k++){
            const auto& o = obs[k];
            auto [dra, ddec] = residual(el, o, dco.iopt);
            auto Hk = obs_partials(el, o, dco.iopt);  // 2×6

            xi[2*k  ] = dra;
            xi[2*k+1] = ddec;
            w[2*k  ]  = o.w_ra();
            w[2*k+1]  = o.w_dec();

            chi2 += dra*dra*o.w_ra() + ddec*ddec*o.w_dec();

            for(int j=0;j<6;j++){
                H(2*k,  j) = Hk(0,j);
                H(2*k+1,j) = Hk(1,j);
            }
        }

        double rms = std::sqrt(chi2 / m);

        if(dco.verbose){
            std::cout << "  iter " << iter
                      << "  rms=" << std::setprecision(6) << rms << " arcsec\n";
        }

        // Normal matrix  C = H^T W H
        Matrix C(n, n, 0.0);
        std::vector<double> b(n, 0.0);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                double s=0;
                for(int k=0;k<m;k++) s+=H(k,i)*w[k]*H(k,j);
                C(i,j)=s;
            }
            double s=0;
            for(int k=0;k<m;k++) s+=H(k,i)*w[k]*xi[k];
            b[i]=s;
        }

        // Levenberg-Marquardt damping
        Matrix Clm = C;
        if(dco.use_lm){
            for(int i=0;i<n;i++) Clm(i,i) += lambda*C(i,i);
        }

        // Solve for correction Δx
        std::vector<double> dx;
        try {
            dx = solve_linear(Clm, b);
        } catch(...) {
            if(dco.use_lm) { lambda *= 10; continue; }
            break;
        }

        // Compute norm of step
        double dnorm = 0;
        for(int i=0;i<6;i++) dnorm += dx[i]*dx[i];
        dnorm = std::sqrt(dnorm);

        // Apply correction
        OrbitalElements el_new = el;
        el_new.a     += dx[0];
        el_new.e     += dx[1];
        el_new.i     += dx[2];
        el_new.Omega += dx[3];
        el_new.omega += dx[4];
        el_new.M0    += dx[5];

        // Sanity checks
        if(el_new.a < 0.01 || el_new.a > 1e4) { lambda*=10; continue; }
        if(el_new.e < 0.0 || el_new.e >= 1.0) { lambda*=10; continue; }

        double rms_new = compute_rms(el_new, obs, dco.iopt);

        if(dco.use_lm){
            if(rms_new < rms){
                lambda = std::max(lambda/3.0, 1e-10);
                el = el_new;
            } else {
                lambda = std::min(lambda*3.0, 1e6);
                continue;
            }
        } else {
            el = el_new;
        }

        // Convergence check
        double drms = std::abs(rms_prev - rms_new);
        if(drms < dco.conv_rms && dnorm < dco.conv_delta){
            converged = true;
            if(dco.verbose) std::cout << "  Converged at iter " << iter+1 << "\n";
            break;
        }
        rms_prev = rms_new;
    }

    // Recompute final normal matrix and covariance
    {
        Matrix C(n,n,0.0);
        double chi2=0;
        for(int k=0;k<nobs;k++){
            const auto& o=obs[k];
            auto [dra,ddec]=residual(el,o,dco.iopt);
            auto Hk=obs_partials(el,o,dco.iopt);
            chi2+=dra*dra*o.w_ra()+ddec*ddec*o.w_dec();
            for(int i=0;i<n;i++)
                for(int j=0;j<n;j++){
                    C(i,j)+=Hk(0,i)*o.w_ra()*Hk(0,j)
                            +Hk(1,i)*o.w_dec()*Hk(1,j);
                }
        }
        res.normal = C;
        try { res.covariance = invert(C); }
        catch(...) { res.covariance = Matrix(6,6,0.0); }
        res.chi2 = chi2;
        res.rms  = std::sqrt(chi2/(2*nobs));
    }

    res.elements   = el;
    res.iterations = iter;
    res.converged  = converged;
    return res;
}

} // namespace orbfit
