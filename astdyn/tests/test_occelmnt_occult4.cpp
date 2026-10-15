/**
 * @file test_occelmnt_occult4.cpp
 * @brief The occelmnt XML we write must survive Occult4's own reader.
 *
 * OccultWatcher and Occult4 exchange this file, and Occult4's reader
 * (OccultationElements.Read_XMLElements and the Parse_* methods) works line by
 * line and calls double.Parse / int.Parse on most fields: an empty field is an
 * exception, not a missing value. The rules below are transcribed from those
 * methods; the writer is checked against them on a complete event and on the
 * worst case (no B/R magnitudes, no covariance).
 */

#include <gtest/gtest.h>
#include "astdyn/io/OccultationXMLIO.hpp"
#include <cstdlib>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using astdyn::astrometry::OccultationEvent;
using astdyn::io::OccultationXMLIO;

namespace {

enum class Rule { Double, Int, Any };

/// Which fields of each tag Occult4 parses strictly.
const std::map<std::string, std::vector<Rule>> kOccult4Rules = {
    {"Elements", {Rule::Any, Rule::Double, Rule::Int, Rule::Int, Rule::Int, Rule::Double, Rule::Double,
                  Rule::Double, Rule::Double, Rule::Double, Rule::Double, Rule::Double, Rule::Double, Rule::Double}},
    {"Earth", {Rule::Double, Rule::Double, Rule::Double, Rule::Double, Rule::Any}},
    {"Star", {Rule::Any, Rule::Double, Rule::Double, Rule::Any, Rule::Any, Rule::Any, Rule::Double, Rule::Int,
              Rule::Any, Rule::Double, Rule::Double, Rule::Double, Rule::Double, Rule::Any, Rule::Int, Rule::Int}},
    {"Object", {Rule::Any, Rule::Any, Rule::Double, Rule::Double, Rule::Double, Rule::Int, Rule::Int,
                Rule::Double, Rule::Double, Rule::Any, Rule::Double, Rule::Any, Rule::Double, Rule::Double}},
    {"Orbit", {Rule::Double, Rule::Double, Rule::Int, Rule::Int, Rule::Double, Rule::Double, Rule::Double,
               Rule::Double, Rule::Double, Rule::Double, Rule::Double, Rule::Double, Rule::Double, Rule::Double}},
    {"Errors", {Rule::Double, Rule::Double, Rule::Double, Rule::Double, Rule::Double, Rule::Any, Rule::Double,
                Rule::Int, Rule::Int, Rule::Int}},
    {"ID", {Rule::Any, Rule::Double}},
};

bool parses_as_double(const std::string& field) {
    if (field.empty()) return false;
    char* end = nullptr;
    std::strtod(field.c_str(), &end);
    return *end == '\0';
}

bool parses_as_int(const std::string& field) {
    if (field.empty()) return false;
    char* end = nullptr;
    std::strtol(field.c_str(), &end, 10);
    return *end == '\0';
}

std::vector<std::string> split_fields(const std::string& csv) {
    std::vector<std::string> out;
    std::string item;
    std::stringstream ss(csv);
    while (std::getline(ss, item, ',')) out.push_back(item);
    if (!csv.empty() && csv.back() == ',') out.push_back("");
    return out;
}

/// Every violation of Occult4's reader rules found in @p xml, one per line.
std::vector<std::string> occult4_violations(const std::string& xml) {
    std::vector<std::string> problems;
    for (const auto& [tag, rules] : kOccult4Rules) {
        const std::string open = "<" + tag + ">", close = "</" + tag + ">";
        const auto a = xml.find(open), b = xml.find(close);
        if (a == std::string::npos || b == std::string::npos) { problems.push_back(tag + ": missing"); continue; }
        const auto fields = split_fields(xml.substr(a + open.size(), b - a - open.size()));
        for (size_t i = 0; i < rules.size() && i < fields.size(); ++i) {
            const bool ok = rules[i] == Rule::Any || (rules[i] == Rule::Double ? parses_as_double(fields[i])
                                                                                : parses_as_int(fields[i]));
            if (!ok) problems.push_back(tag + "[" + std::to_string(i) + "] = '" + fields[i] + "'");
        }
    }
    return problems;
}

OccultationEvent complete_event() {
    OccultationEvent ev;
    ev.elements_source = "AstDyn";
    ev.duration_s = 0.4271;
    ev.year = 2027; ev.month = 10; ev.day = 21; ev.ut_closest_h = 16.9513;
    ev.x = 0.47; ev.y = 1.12; ev.dx = 23.7; ev.dy = -9.9;
    ev.star_id = "J161525.81-211018.4";
    ev.star_ra_h = 16.257; ev.star_dec_deg = -21.17; ev.star_app_ra_h = 16.26; ev.star_app_dec_deg = -21.2;
    ev.mag_b = 10.1; ev.mag_v = 9.59; ev.mag_r = 9.0; ev.mag_drop_v = 7.91; ev.mag_drop_r = 7.5;
    ev.object_number = "3830"; ev.object_name = "Trelleborg"; ev.object_mag = 17.5;
    ev.mag_v_asteroid = 17.5; ev.mag_r_asteroid = 17.05;
    ev.diameter_km = 17.789; ev.distance_au = 4.17;
    ev.epoch_year = 2026; ev.epoch_month = 6; ev.epoch_day = 9;
    ev.err_path_widths = 1.53; ev.err_major_as = 0.022; ev.err_minor_as = 0.010; ev.err_pa_deg = 80;
    ev.error_basis = "Star+PEU"; ev.reliability = 1.03;
    ev.event_id = "20271021_1018.4"; ev.prediction_mjd = 61324.79;
    return ev;
}

OccultationEvent worst_case_event() {
    OccultationEvent ev = complete_event();
    ev.mag_b = 99.0; ev.mag_r = 99.0;            // catalogue with G only
    ev.mag_drop_r = 0.0;
    ev.mag_v_asteroid = 0.0; ev.mag_r_asteroid = 0.0;
    ev.err_path_widths = 0.0; ev.err_major_as = 0.0; ev.err_minor_as = 0.0;
    ev.error_basis = "Star+Assumed"; ev.reliability = -1.0;
    return ev;
}

}  // namespace

TEST(Occelmnt, CompleteEventSurvivesOccult4Reader) {
    const auto problems = occult4_violations(OccultationXMLIO::write_string({complete_event()}));
    EXPECT_TRUE(problems.empty()) << problems.front();
}

TEST(Occelmnt, WorstCaseEventSurvivesOccult4Reader) {
    const auto problems = occult4_violations(OccultationXMLIO::write_string({worst_case_event()}));
    EXPECT_TRUE(problems.empty()) << problems.front();
}

TEST(Occelmnt, ObjectLineHasOccult4TrailingCommaAndAsteroidMagnitudes) {
    // Occult4 writes "{0:f2}," after MagR too; with more than 13 fields its
    // reader parses fields 12 and 13 as doubles.
    const std::string xml = OccultationXMLIO::write_string({complete_event()});
    EXPECT_NE(xml.find(",17.50,17.05,</Object>"), std::string::npos) << xml;
}

TEST(Occelmnt, LinesEndWithCrlfAndNoXmlDeclaration) {
    const std::string xml = OccultationXMLIO::write_string({complete_event()});
    EXPECT_EQ(xml.rfind("<Occultations>\r\n", 0), 0u);
    EXPECT_EQ(xml.find("<?xml"), std::string::npos);
    size_t lf = 0, crlf = 0;
    for (size_t i = 0; i < xml.size(); ++i) {
        if (xml[i] != '\n') continue;
        ++lf;
        if (i > 0 && xml[i - 1] == '\r') ++crlf;
    }
    EXPECT_EQ(lf, crlf);
}

TEST(Occelmnt, DecimalsFollowOccult4) {
    const std::string xml = OccultationXMLIO::write_string({complete_event()});
    EXPECT_NE(xml.find("<Elements>AstDyn,0.43,2027,10,21,"), std::string::npos);
    EXPECT_NE(xml.find(",16.25700000,-21.1700000,10.10,9.59,9.00,0.0,"), std::string::npos) << xml;
}

TEST(Occelmnt, MissingBlueMagnitudeDerivedWithOccult4Relation) {
    auto ev = complete_event();
    ev.mag_b = 99.0;                               // B from V and R: (1.54 V - R) / 0.54
    const std::string xml = OccultationXMLIO::write_string({ev});
    EXPECT_NE(xml.find(",10.68,9.59,9.00,"), std::string::npos) << xml;
}

TEST(Occelmnt, PathErrorRecordCarriesOccult4OffsetOfOne) {
    auto ev = complete_event();
    ev.err_path_widths = 0.020;                    // OWC showed it as -0.980 without the offset
    const std::string xml = OccultationXMLIO::write_string({ev});
    EXPECT_NE(xml.find("<Errors>1.020,"), std::string::npos) << xml;
    EXPECT_NEAR(OccultationXMLIO::read_string(xml)[0].err_path_widths, 0.020, 1e-9);
}

TEST(Occelmnt, UncomputedPathErrorStaysZero) {
    const std::string xml = OccultationXMLIO::write_string({worst_case_event()});
    EXPECT_NE(xml.find("<Errors>0.000,"), std::string::npos) << xml;
}

TEST(Occelmnt, RoundTripKeepsTheGeometry) {
    const auto back = OccultationXMLIO::read_string(OccultationXMLIO::write_string({complete_event()}));
    ASSERT_EQ(back.size(), 1u);
    EXPECT_NEAR(back[0].ut_closest_h, 16.9513, 1e-7);
    EXPECT_NEAR(back[0].dy, -9.9, 1e-7);
    EXPECT_EQ(back[0].object_name, "Trelleborg");
    EXPECT_NEAR(back[0].mag_r_asteroid, 17.05, 1e-9);
}
