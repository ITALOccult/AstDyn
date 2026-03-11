/**
 * OrbFit C++ — Self-contained verification test (pure two-body)
 * No planetary perturbations → fast execution
 *
 * Build: g++ -O3 -std=c++17 -I../include -o test_2body test_2body.cpp
 */

#include "constants.h"
#include "linalg.h"
#include "orbital_elements.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <vector>

using namespace orbfit;

// ── Two-body propagation (no planet calls) ────────────────────
using State6 = std::array<double,6>;

State6 eom_2body(const State6& s, double /*t*/) {
    Vec3 r{s[0],s[1],s[2]};
    double rv3 = std::pow(r.norm(),3.0);
    Vec3 acc = r*(-GM_SUN/rv3);
    Vec3 v{s[3],s[4],s[5]};
    return {v.x,v.y,v.z,acc.x,acc.y,acc.z};
}

State6 rk4_2body(const State6& y, double t, double h) {
    auto add=[](const State6& a, const State6& b, double s)->State6{
        State6 r; for(int i=0;i<6;i++) r[i]=a[i]+b[i]*s; return r;
    };
    auto k1=eom_2body(y,t);
    auto k2=eom_2body(add(y,k1,h/2),t+h/2);
    auto k3=eom_2body(add(y,k2,h/2),t+h/2);
    auto k4=eom_2body(add(y,k3,h),  t+h  );
    State6 yn; for(int i=0;i<6;i++) yn[i]=y[i]+h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    return yn;
}

StateVector propagate_2body(const StateVector& sv0, double t_end, double h=0.1) {
    State6 y={sv0.pos.x,sv0.pos.y,sv0.pos.z,sv0.vel.x,sv0.vel.y,sv0.vel.z};
    double t=sv0.epoch;
    double sgn=(t_end>t)?1.0:-1.0;
    h=std::abs(h)*sgn;
    while(sgn*(t_end-t)>1e-10){
        if(std::abs(h)>std::abs(t_end-t)) h=t_end-t;
        y=rk4_2body(y,t,h);
        t+=h;
    }
    StateVector sv; sv.epoch=t_end;
    sv.pos=Vec3(y[0],y[1],y[2]); sv.vel=Vec3(y[3],y[4],y[5]);
    return sv;
}

// ── Observation model (2-body, heliocentric observer approx) ──
struct Obs2 {
    double t, ra, dec;
    double sra, sdec;
    Vec3 obs_pos;
};

std::pair<double,double> predict_2body(const OrbitalElements& el, const Obs2& o){
    auto sv=propagate_2body(elements_to_state(el,el.epoch), o.t, 0.05);
    Vec3 rho=sv.pos-o.obs_pos;
    double rn=rho.norm();
    // ecliptic→equatorial
    const double eps=23.4392911*DEG2RAD;
    double ce=std::cos(eps),se=std::sin(eps);
    double xeq=rho.x, yeq=ce*rho.y-se*rho.z, zeq=se*rho.y+ce*rho.z;
    double ra=std::atan2(yeq,xeq); if(ra<0) ra+=2*M_PI;
    double dec=std::asin(zeq/rn);
    return {ra,dec};
}

// Residual [arcsec]
std::pair<double,double> resid(const OrbitalElements& el, const Obs2& o){
    auto [pra,pdec]=predict_2body(el,o);
    double dra=(o.ra-pra)*std::cos(o.dec);
    double ddec=o.dec-pdec;
    if(dra> M_PI) dra-=2*M_PI;
    if(dra<-M_PI) dra+=2*M_PI;
    return {dra*RAD2ARCSEC, ddec*RAD2ARCSEC};
}

// 2×6 partials
Matrix partials_2body(const OrbitalElements& el, const Obs2& o){
    const double EPS=1e-5;
    Matrix H(2,6);
    auto pack=[](const OrbitalElements& e)->std::array<double,6>{
        return {e.a,e.e,e.i,e.Omega,e.omega,e.M0};};
    auto unpack=[](std::array<double,6> p, OrbitalElements b)->OrbitalElements{
        b.a=p[0];b.e=p[1];b.i=p[2];b.Omega=p[3];b.omega=p[4];b.M0=p[5];return b;};
    auto p0=pack(el);
    auto [pra0,pdec0]=predict_2body(el,o);
    double cd=std::cos(o.dec);
    for(int j=0;j<6;j++){
        auto pp=p0,pm=p0;
        double dp=EPS*std::max(std::abs(p0[j]),1e-9);
        pp[j]+=dp; pm[j]-=dp;
        auto [rp,dp2]=predict_2body(unpack(pp,el),o);
        auto [rm,dm] =predict_2body(unpack(pm,el),o);
        H(0,j)=(rp*cd-rm*cd)/(2*dp)*RAD2ARCSEC;
        H(1,j)=(dp2-dm)/(2*dp)*RAD2ARCSEC;
    }
    return H;
}

double rms_2body(const OrbitalElements& el, const std::vector<Obs2>& obs){
    double s=0; int n=0;
    for(auto& o:obs){auto[a,b]=resid(el,o);s+=a*a+b*b;n+=2;}
    return std::sqrt(s/n);
}

// ── DC (Gauss-Newton) ──────────────────────────────────────────
OrbitalElements dc_2body(OrbitalElements el, const std::vector<Obs2>& obs,
                          int maxiter=20, bool verbose=true)
{
    int m=2*obs.size(), n=6;
    for(int it=0;it<maxiter;it++){
        Matrix C(n,n,0.0); std::vector<double> b(n,0.0);
        double chi2=0;
        for(int k=0;k<(int)obs.size();k++){
            auto[dra,ddec]=resid(el,obs[k]);
            auto H=partials_2body(el,obs[k]);
            double wr=1/(obs[k].sra*obs[k].sra*RAD2ARCSEC*RAD2ARCSEC);
            double wd=1/(obs[k].sdec*obs[k].sdec*RAD2ARCSEC*RAD2ARCSEC);
            chi2+=dra*dra*wr+ddec*ddec*wd;
            for(int i=0;i<n;i++) for(int j=0;j<n;j++)
                C(i,j)+=H(0,i)*wr*H(0,j)+H(1,i)*wd*H(1,j);
            for(int i=0;i<n;i++) b[i]+=H(0,i)*wr*dra+H(1,i)*wd*ddec;
        }
        double rms=std::sqrt(chi2/m);
        if(verbose) std::cout<<"  iter "<<it<<"  rms="<<std::setprecision(5)<<rms<<" arcsec\n";
        auto dx=solve_linear(C,b);
        el.a+=dx[0]; el.e+=dx[1]; el.i+=dx[2];
        el.Omega+=dx[3]; el.omega+=dx[4]; el.M0+=dx[5];
        double dnorm=0; for(auto x:dx) dnorm+=x*x;
        if(std::sqrt(dnorm)<1e-12 && rms<1e-3) break;
    }
    return el;
}

// ── TESTS ─────────────────────────────────────────────────────
void test_kepler(){
    std::cout<<"[1] Kepler equation... ";
    for(double e:{0.0,0.2,0.7,0.99})
        for(double M:{0.1,1.0,2.5,5.5}){
            double E=solve_kepler(M,e);
            double M2=E-e*std::sin(E);
            double Mn=std::fmod(M,2*M_PI); if(Mn<0)Mn+=2*M_PI;
            assert(std::abs(M2-Mn)<1e-11);
        }
    std::cout<<"PASS\n";
}

void test_roundtrip(){
    std::cout<<"[2] Elements↔State roundtrip... ";
    OrbitalElements el; el.a=2.5;el.e=0.15;el.i=10*DEG2RAD;
    el.Omega=120*DEG2RAD;el.omega=250*DEG2RAD;el.M0=30*DEG2RAD;el.epoch=JD_J2000;
    auto sv=elements_to_state(el,JD_J2000);
    auto el2=state_to_elements(sv);
    assert(std::abs(el.a-el2.a)<1e-10);
    assert(std::abs(el.e-el2.e)<1e-10);
    assert(std::abs(el.i-el2.i)<1e-8);
    std::cout<<"PASS  (Δa="<<std::abs(el.a-el2.a)<<")\n";
}

void test_energy(){
    std::cout<<"[3] Energy conservation 1 yr... ";
    OrbitalElements el; el.a=1.5;el.e=0.2;el.i=5*DEG2RAD;
    el.Omega=0;el.omega=0;el.M0=0;el.epoch=JD_J2000;
    auto sv0=elements_to_state(el,JD_J2000);
    auto sv1=propagate_2body(sv0,JD_J2000+365.25,0.05);
    double E0=0.5*sv0.vel.norm2()-GM_SUN/sv0.pos.norm();
    double E1=0.5*sv1.vel.norm2()-GM_SUN/sv1.pos.norm();
    double dE=std::abs(E1-E0)/std::abs(E0);
    assert(dE<1e-7);
    std::cout<<"PASS  (ΔE/E="<<std::scientific<<dE<<")\n";
}

void test_dc(){
    std::cout<<"[4] Differential Corrections (pure 2-body)...\n";
    OrbitalElements true_el;
    true_el.a=2.0;true_el.e=0.1;true_el.i=15*DEG2RAD;
    true_el.Omega=80*DEG2RAD;true_el.omega=130*DEG2RAD;true_el.M0=60*DEG2RAD;
    true_el.epoch=JD_J2000;

    // Observer at 1 AU (fixed, simplified)
    Vec3 obs_pos{1.0,0.0,0.0};
    std::vector<Obs2> obs;
    for(int k=0;k<10;k++){
        double t=JD_J2000+200+k*3.0;
        auto sv=propagate_2body(elements_to_state(true_el,true_el.epoch),t,0.05);
        Vec3 rho=sv.pos-obs_pos;
        double rn=rho.norm();
        const double eps=23.4392911*DEG2RAD;
        double ce=std::cos(eps),se=std::sin(eps);
        double xeq=rho.x,yeq=ce*rho.y-se*rho.z,zeq=se*rho.y+ce*rho.z;
        double ra=std::atan2(yeq,xeq); if(ra<0)ra+=2*M_PI;
        double dec=std::asin(zeq/rn);
        obs.push_back({t,ra,dec,1.0,1.0,obs_pos});
    }

    // Perturb initial elements slightly
    OrbitalElements el0=true_el;
    el0.a+=0.01; el0.e+=0.005; el0.i+=0.01;

    auto el_fit=dc_2body(el0,obs,15,true);

    double da=std::abs(el_fit.a-true_el.a);
    double de=std::abs(el_fit.e-true_el.e);
    std::cout<<"  Δa="<<std::scientific<<da<<"  Δe="<<de<<"\n";
    assert(da<1e-5 && "a not recovered");
    assert(de<1e-5 && "e not recovered");
    std::cout<<"  PASS\n";
}

int main(){
    std::cout<<"╔══════════════════════════════════════╗\n";
    std::cout<<"║  OrbFit C++ — Two-body unit tests    ║\n";
    std::cout<<"╚══════════════════════════════════════╝\n\n";
    test_kepler();
    test_roundtrip();
    test_energy();
    test_dc();
    std::cout<<"\n✓  All tests passed.\n";
    return 0;
}
