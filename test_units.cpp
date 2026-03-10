
#include <iostream>
#include <iomanip>
#include "astdyn/core/physics_state_au.hpp"
#include "astdyn/core/physics_types.hpp"

using namespace astdyn;
using namespace astdyn::physics;

int main() {
    CartesianStateAU<> state;
    state.x = 1.0; state.y = 0.0; state.z = 0.0;
    state.vx = 0.01; state.vy = 0.0; state.vz = 0.0;
    
    DerivativeAU d;
    d.dx = VelocityAUD::from_au_d(0.1); // AU/d
    d.dvx = AccelerationAUD2::from_au_d2(0.01); // AU/d^2
    
    double h = 0.5;
    DerivativeAU dh = d * h;
    
    std::cout << "d.dx: " << d.dx.to_au_d() << " AU/d" << std::endl;
    std::cout << "dh.dx: " << dh.dx.to_au_d() << " AU" << std::endl;
    
    CartesianStateAU<> next = state + dh;
    
    std::cout << "state.x: " << state.x << std::endl;
    std::cout << "next.x: " << next.x << " (Expected 1.05)" << std::endl;
    std::cout << "next.vx: " << next.vx << " (Expected 0.005)" << std::endl;

    return 0;
}
