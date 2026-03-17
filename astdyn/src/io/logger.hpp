#ifndef ASTDYN_IO_LOGGER_HPP
#define ASTDYN_IO_LOGGER_HPP

#include "src/core/messages.hpp"
#include <iostream>

namespace astdyn::io {

using core::MessageKey;
using core::Language;

/**
 * @brief Centralized logger for the AstDyn engine.
 * CTFYH adherence: No free text strings allowed. Only localized MessageKeys.
 */
class Logger {
public:
    /** @brief Logs a message to stdout based on the key and preferred language. */
    static void info(const MessageKey key, const Language lang = Language::English) noexcept {
        const auto message = core::get_message(key, lang);
        std::cout << "[AstDyn INFO] " << message << "\n";
    }

    /** @brief Logs an error message to stderr. */
    static void error(const MessageKey key, const Language lang = Language::English) noexcept {
        const auto message = core::get_message(key, lang);
        std::cerr << "[AstDyn ERROR] " << message << "\n";
    }
};

} // namespace astdyn::io

#endif // ASTDYN_IO_LOGGER_HPP
