/**
 * @file SPKReader.hpp
 * @brief Native reader for JPL SPK (Binary SPICE Kernel) files
 * @author AstDyn Team
 * 
 * Implements a lightweight, dependency-free reader for DAF/SPK files (Type 2 and 3).
 * Compatible with JPL DE4xx ephemerides.
 */

#ifndef ASTDYN_IO_SPK_READER_HPP
#define ASTDYN_IO_SPK_READER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <memory>
#include <Eigen/Dense>

namespace astdyn::io {

/**
 * @brief Represents a segment within an SPK file
 */
struct SPKSegment {
    int body_id;            ///< Target body NAIF ID
    int center_id;          ///< Center body NAIF ID (e.g. 0 for SSB)
    int frame_id;           ///< Reference frame ID (1=J2000)
    int type;               ///< Data type (2=Chebyshev Pos, 3=Chebyshev Pos+Vel)
    double start_et;        ///< Start epoch (seconds past J2000)
    double end_et;          ///< End epoch (seconds past J2000)
    int start_addr;         ///< Start address of data in file (1-based record index -> needs conversion)
    int end_addr;           ///< End address
    
    // Type 2/3 specific parameters
    double init_sec;        ///< Initial epoch of first record
    double intlen;          ///< Length of interval in seconds
    int rsize;              ///< Record size (number of coeffs)
    int order;              ///< Polynomial order (N)
    int n_comp;             ///< Number of components (3)
};

/**
 * @brief Reader for binary SPK files (DAF format)
 */
class SPKReader {
public:
    explicit SPKReader(const std::string& filename);
    ~SPKReader();

    /**
     * @brief Get position and velocity of a target body relative to a center
     * @param target_id NAIF ID of target body
     * @param et Ephemeris time (seconds past J2000 TDB)
     * @return State vector [x, y, z, vx, vy, vz] in km and km/s
     */
    Eigen::VectorXd getState(int target_id, double et);

    bool isOpen() const { return file_.is_open(); }

private:
    std::ifstream file_;
    std::string filename_;
    bool big_endian_; // File endianness

    // Index of segments mapped by target ID
    // A body might have multiple segments covering different time ranges
    std::multimap<int, SPKSegment> segments_;

    void loadIndex();
    void readHeader();
    
    // Low-level read functions
    void readDouble(double& out);
    void readInt(int& out);
    void seekToRecord(int record_idx); // 1-based record index (1 record = 1024 bytes)
    
    // DAF internal parsing
    void parseDAFRecord(const std::vector<char>& buffer);

    // Chebyshev evaluation
    Eigen::VectorXd evaluateType2(const SPKSegment& seg, double et);
};

} // namespace astdyn::io

#endif // ASTDYN_IO_SPK_READER_HPP
