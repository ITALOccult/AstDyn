#ifndef ASTDYN_CONFIG_HPP
#define ASTDYN_CONFIG_HPP

// Build configuration
#define ASTDYN_USE_SPICE
/* #undef ASTDYN_ENABLE_PROFILING */

// Build type
#define ASTDYN_BUILD_TYPE "Release"

namespace astdyn {

struct Config {
    static constexpr const char* build_type = ASTDYN_BUILD_TYPE;
    
#ifdef ASTDYN_USE_SPICE
    static constexpr bool use_spice = true;
#else
    static constexpr bool use_spice = false;
#endif

#ifdef ASTDYN_ENABLE_PROFILING
    static constexpr bool profiling_enabled = true;
#else
    static constexpr bool profiling_enabled = false;
#endif
};

} // namespace astdyn

#endif // ASTDYN_CONFIG_HPP
