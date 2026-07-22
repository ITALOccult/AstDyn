/**
 * @file OccultationJSONIO.hpp
 * @brief Serializzazione degli eventi di occultazione in JSON nativo.
 *
 * Formato interno del nostro orchestratore: contiene TUTTI i campi di
 * OccultationEvent (piu' ricco dell'XML Occult4, che ne espone solo alcuni in
 * forma posizionale). L'orchestratore Python lo usa per lo screening/raffinamento
 * (campo object.number esplicito -> estrazione robusta dei positivi) e come
 * formato di export alternativo all'XML.
 */
#ifndef ASTDYN_IO_OCCULTATION_JSON_IO_HPP
#define ASTDYN_IO_OCCULTATION_JSON_IO_HPP

#include <string>
#include <vector>
#include "astdyn/astrometry/OccultationEvent.hpp"

namespace astdyn::io {

using namespace astdyn::astrometry;

class OccultationJSONIO {
public:
    /// Serializza gli eventi in una stringa JSON (oggetto con chiave "events").
    static std::string write_string(const std::vector<OccultationEvent>& events);

    /// Scrive gli eventi in un file JSON. Ritorna true se ok.
    static bool write_file(const std::vector<OccultationEvent>& events,
                           const std::string& filename);
};

} // namespace astdyn::io

#endif // ASTDYN_IO_OCCULTATION_JSON_IO_HPP
