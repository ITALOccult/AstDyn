#ifndef ASTDYN_IO_SPKREADER_HPP
#define ASTDYN_IO_SPKREADER_HPP

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <Eigen/Dense>

namespace astdyn::io {

/**
 * @brief Native SPK (BSP) file reader
 * 
 * Supports Type 2 (Chebyshev) and Type 13 (Hermite) segments.
 * Independent of CSPICE library.
 */
class SPKReader {
public:
    explicit SPKReader(const std::string& filename);
    ~SPKReader();

    /**
     * @brief Get state vector (Pos, Vel) for target body at time ET
     * @param target_id NAIF ID of the body
     * @param et Ephemeris time (seconds from J2000 epoch)
     * @return 6D vector [km, km/s] in Ecliptic J2000 frame
     */
    Eigen::Matrix<double, 6, 1> getState(int target_id, double et);

private:
    struct SPKSegment {
        int body_id;
        int center_id;
        int frame_id;
        int type;
        int start_addr;
        int end_addr;
        double start_et;
        double end_et;
        
        // Type specific params
        double init_sec = 0;
        double intlen = 0;
        int rsize = 0;
        int order = 0;
        int n_comp = 3;
    };

    struct Type2Cache {
        const SPKSegment* seg = nullptr;
        int rec_idx = -1;
        std::vector<double> coeffs;
    };

    std::string filename_;
    std::ifstream file_;
    bool big_endian_ = false;
    std::multimap<int, SPKSegment> segments_;
    std::map<int, Type2Cache> type2_cache_;

    void readHeader();
    void loadIndex();
    void seekToRecord(int record_idx);
    
    
    Eigen::Matrix<double, 6, 1> evaluateType2(const SPKSegment& seg, double et);
    Eigen::Matrix<double, 6, 1> evaluateType13(const SPKSegment& seg, double et);

    // Buffers to avoid allocations
    std::vector<double> T_buf_;
    std::vector<double> U_buf_;
};

} // namespace astdyn::io

#endif // ASTDYN_IO_SPKREADER_HPP
