#ifndef ASTDYN_VERSION_HPP
#define ASTDYN_VERSION_HPP

#define ASTDYN_VERSION_MAJOR 1
#define ASTDYN_VERSION_MINOR 0
#define ASTDYN_VERSION_PATCH 0
#define ASTDYN_VERSION_STRING "1.0.0"

namespace astdyn {

struct Version {
    static constexpr int major = ASTDYN_VERSION_MAJOR;
    static constexpr int minor = ASTDYN_VERSION_MINOR;
    static constexpr int patch = ASTDYN_VERSION_PATCH;
    static constexpr const char* string = ASTDYN_VERSION_STRING;
};

} // namespace astdyn

#endif // ASTDYN_VERSION_HPP
