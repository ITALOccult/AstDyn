/**
 * @file OccultationXMLIO.hpp
 * @brief I/O handler for the Occult4-style XML format for occultation events.
 */

#ifndef ASTDYN_IO_OCCULTATIONXMLIO_HPP
#define ASTDYN_IO_OCCULTATIONXMLIO_HPP

#include "astdyn/astrometry/OccultationEvent.hpp"
#include <string>
#include <vector>

namespace astdyn::io {

using namespace astdyn::astrometry;

/**
 * @brief Class for reading and writing OccultationEvent data in XML format.
 */
class OccultationXMLIO {
public:
    /**
     * @brief Reads a list of occultation events from an XML string.
     * @param xml_content The XML content as a string.
     * @return A vector of OccultationEvent structs.
     */
    static std::vector<OccultationEvent> read_string(const std::string& xml_content);

    /**
     * @brief Reads a list of occultation events from an XML file.
     * @param filename Path to the XML file.
     * @return A vector of OccultationEvent structs.
     */
    static std::vector<OccultationEvent> read_file(const std::string& filename);

    /**
     * @brief Writes a list of occultation events to an XML string.
     * @param events The list of events to write.
     * @return A string containing the formatted XML.
     */
    static std::string write_string(const std::vector<OccultationEvent>& events);

    /**
     * @brief Writes a list of occultation events to an XML file.
     * @param events The list of events to write.
     * @param filename Path to the output XML file.
     * @return True if successful, false otherwise.
     */
    static bool write_file(const std::vector<OccultationEvent>& events, const std::string& filename);

private:
    static OccultationEvent parse_event_node(const std::string& event_xml);
    static std::string format_event_node(const OccultationEvent& event);
    static std::vector<std::string> split_csv(const std::string& csv);
};

} // namespace astdyn::io

#endif // ASTDYN_IO_OCCULTATIONXMLIO_HPP
