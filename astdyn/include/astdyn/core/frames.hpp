#ifndef ASTDYN_CORE_FRAMES_HPP
#define ASTDYN_CORE_FRAMES_HPP

namespace astdyn::core {

/**
 * @brief Enumeration of supported reference frames.
 */
enum class ReferenceFrame {
    GCRF,
    ITRF,
    TEME,
    ECLIPJ2000
};

/**
 * @brief Contextual information for a reference frame at a specific epoch.
 */
struct FrameContext {
    const ReferenceFrame frame;
    const double epoch_mjd;

    explicit constexpr FrameContext(const ReferenceFrame f, const double mjd) 
        : frame(f), epoch_mjd(mjd) {}
};

} // namespace astdyn::core

#endif // ASTDYN_CORE_FRAMES_HPP
