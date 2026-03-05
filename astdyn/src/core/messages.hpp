#ifndef ASTDYN_CORE_MESSAGES_HPP
#define ASTDYN_CORE_MESSAGES_HPP

#include <string_view>

namespace astdyn::core {

enum class Language {
    English,
    Italian
};

enum class MessageKey {
    ConvergenceError,
    StepSizeMinimumReached,
    MaxStepsExceeded,
    EnergyDriftDetected,
    IntegrationStarted,
    BackPropagationNotSupported,
    FileOpenError,
    InvalidFormat,
    ValueOutOfRange
};

/**
 * @brief Retrieves a localized message string based on the provided key and language.
 * Follows the "Fit in Your Head" protocol: no hardcoded strings in logic.
 */
[[nodiscard]] constexpr std::string_view get_message(const MessageKey key, const Language lang) noexcept {
    if (lang == Language::Italian) {
        switch (key) {
            case MessageKey::ConvergenceError: return "Errore di convergenza nel sistema implicito.";
            case MessageKey::StepSizeMinimumReached: return "Passo di integrazione inferiore al minimo consentito.";
            case MessageKey::MaxStepsExceeded: return "Raggiunto il numero massimo di passi di integrazione.";
            case MessageKey::EnergyDriftDetected: return "Rilevata deriva energetica significativa.";
            case MessageKey::IntegrationStarted: return "Inizio della propagazione orbitale.";
            case MessageKey::BackPropagationNotSupported: return "Propagazione all'indietro non supportata da questo motore.";
            case MessageKey::FileOpenError: return "Impossibile aprire il file richiesto.";
            case MessageKey::InvalidFormat: return "Formato del file non valido o corrotto.";
            case MessageKey::ValueOutOfRange: return "Valore fisico fuori dai limiti consentiti.";
        }
    }

    switch (key) {
        case MessageKey::ConvergenceError: return "Convergence error in the implicit system.";
        case MessageKey::StepSizeMinimumReached: return "Integration step size fell below minimum.";
        case MessageKey::MaxStepsExceeded: return "Maximum number of integration steps reached.";
        case MessageKey::EnergyDriftDetected: return "Significant energy drift detected.";
        case MessageKey::IntegrationStarted: return "Starting orbital propagation.";
        case MessageKey::BackPropagationNotSupported: return "Back-propagation not supported by this engine.";
        case MessageKey::FileOpenError: return "Failed to open the requested file.";
        case MessageKey::InvalidFormat: return "Invalid or corrupted file format.";
        case MessageKey::ValueOutOfRange: return "Physical value out of acceptable range.";
    }

    return "Unknown Message";
}

} // namespace astdyn::core

#endif // ASTDYN_CORE_MESSAGES_HPP
