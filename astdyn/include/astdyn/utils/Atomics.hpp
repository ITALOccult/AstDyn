#ifndef ASTDYN_UTILS_ATOMICS_HPP
#define ASTDYN_UTILS_ATOMICS_HPP

#include <atomic>
#include <algorithm>

namespace astdyn::utils {

/**
 * @brief Thread-safe minimum update for an atomic variable.
 */
template<typename T>
inline void atomic_min(std::atomic<T>& atomic_val, T val) {
    T old_val = atomic_val.load();
    while (val < old_val && !atomic_val.compare_exchange_weak(old_val, val));
}

/**
 * @brief Thread-safe maximum update for an atomic variable.
 */
template<typename T>
inline void atomic_max(std::atomic<T>& atomic_val, T val) {
    T old_val = atomic_val.load();
    while (val > old_val && !atomic_val.compare_exchange_weak(old_val, val));
}

} // namespace astdyn::utils

#endif
