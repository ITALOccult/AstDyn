/**
 * @file SPKReader.cpp
 * @brief Implementation of Native SPK Reader
 */

#include "astdyn/io/SPKReader.hpp"
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cmath>

namespace astdyn::io {

// DAF standard record size in bytes
// DAF standard record size in bytes
constexpr int RECORD_SIZE = 1024;
// constexpr int DBL_SIZE = 8; // Unused

SPKReader::SPKReader(const std::string& filename) : filename_(filename) {
    file_.open(filename, std::ios::binary);
    if (!file_) {
        throw std::runtime_error("Could not open SPK file: " + filename);
    }
    readHeader();
    loadIndex();
}

SPKReader::~SPKReader() {
    if (file_.is_open()) file_.close();
}

// Helper to determine system endianness
static bool isSystemBigEndian() {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};
    return bint.c[0] == 1;
}

// Swap bytes for endian conversion
template <typename T>
T swapEndian(T u) {
    union {
        T u;
        unsigned char u8[sizeof(T)];
    } source, dest;
    source.u = u;
    for (size_t k = 0; k < sizeof(T); k++)
        dest.u8[k] = source.u8[sizeof(T) - k - 1];
    return dest.u;
}

void SPKReader::readHeader() {
    // Read file record (first 1024 bytes)
    std::vector<char> buffer(RECORD_SIZE);
    file_.read(buffer.data(), RECORD_SIZE);
    
    // Check FTP string at byte 700 to detect endianness
    // Or check standard DAF architecture string
    std::string arch(buffer.data(), 8);
    // Usually "DAF/SPK "
    
    // Heuristic: Read nd (num doubles) and ni (num ints) at specific offsets
    // DAF file record structure:
    // Bytes 0-7: ID word
    // Bytes 8-11: ND (int)
    // Bytes 12-15: NI (int)
    
    int nd, ni;
    std::memcpy(&nd, buffer.data() + 8, 4);
    std::memcpy(&ni, buffer.data() + 12, 4);
    
    // Sanity check: ND usually 2, NI usually 6 for SPK
    if (nd == 2 && ni == 6) {
        big_endian_ = isSystemBigEndian(); // File matches system
    } else {
        // Try swapping
        int nd_swapped = swapEndian(nd);
        int ni_swapped = swapEndian(ni);
        if (nd_swapped == 2 && ni_swapped == 6) {
            big_endian_ = !isSystemBigEndian();
        } else {
            // Try legacy LocID check or assume LTL-IEEE (Little Endian)
             // Defaulting to Little Endian for PC compatibility
             big_endian_ = false;
    }
    }
}

void SPKReader::loadIndex() {
    file_.seekg(0, std::ios::beg);
    std::vector<char> record(RECORD_SIZE);
    file_.read(record.data(), RECORD_SIZE);
    
    // DAF Record 1 layout checks
    int fwd;
    std::memcpy(&fwd, record.data() + 76, 4);
    if (big_endian_ != isSystemBigEndian()) fwd = swapEndian(fwd);
    
    if (fwd <= 0) fwd = 2; // Fallback
    
    int current_rec = fwd;
    // int total_segments = 0; // Unused
    
    while (current_rec > 0) {
        seekToRecord(current_rec);
        file_.read(record.data(), RECORD_SIZE);
        
        double next_d_ptr, ns_d;
        std::memcpy(&next_d_ptr, record.data(), 8);
        std::memcpy(&ns_d, record.data() + 16, 8);
        
        if (big_endian_ != isSystemBigEndian()) {
            next_d_ptr = swapEndian(next_d_ptr);
            ns_d = swapEndian(ns_d);
        }
        
        int next = (int)next_d_ptr;
        int ns = (int)ns_d;
        
        // Safety break
        if (ns > 250 || ns < 0) break; 
        
        int sum_size = 2 * 8 + 6 * 4; // 40 bytes
        
        for (int i = 0; i < ns; ++i) {
            int offset = 24 + i * sum_size;
            
            double start_et, end_et;
            std::memcpy(&start_et, record.data() + offset, 8);
            std::memcpy(&end_et, record.data() + offset + 8, 8);
            
            if (big_endian_ != isSystemBigEndian()) {
                start_et = swapEndian(start_et);
                end_et = swapEndian(end_et);
            }
            
            int ints[6];
            std::memcpy(ints, record.data() + offset + 16, 6 * 4);
            
            if (big_endian_ != isSystemBigEndian()) {
                for (int& v : ints) v = swapEndian(v);
            }
            
            SPKSegment seg;
            seg.start_et = start_et;
            seg.end_et = end_et;
            seg.body_id = ints[0];
            seg.center_id = ints[1];
            seg.frame_id = ints[2];
            seg.type = ints[3];
            seg.start_addr = ints[4];
            seg.end_addr = ints[5];
            
            // Read segment params (lazy or eager)
            // Eager reading for Type 2
            if (seg.type == 2) {
                // int end_byte = (seg.end_addr - 1) * 8; // Unused
                // Read last 4 doubles
                // Check bounds?
                 
                double params[4];
                file_.seekg((seg.end_addr - 4) * 8, std::ios::beg);
                char pbuf[32];
                file_.read(pbuf, 32);
                std::memcpy(params, pbuf, 32);
                
                if (big_endian_ != isSystemBigEndian()) {
                     for (double& v : params) v = swapEndian(v);
                }
                
                seg.init_sec = params[0];
                seg.intlen = params[1];
                seg.rsize = (int)params[2];
                int n_coeffs = (seg.rsize - 2) / 3;
                seg.order = n_coeffs - 1;
                seg.n_comp = 3;
            }
            // Use else if for Type 13 if specific initialization needed, but currently Type 13 is lazy.
            else if (seg.type == 13) {
                seg.n_comp = 3;
            }
            
            segments_.insert({seg.body_id, seg});
        }
        
        current_rec = next;
    }
}

void SPKReader::seekToRecord(int record_idx) {
    file_.seekg((record_idx - 1) * RECORD_SIZE, std::ios::beg);
}

Eigen::VectorXd SPKReader::getState(int target_id, double et) {
    auto range = segments_.equal_range(target_id);
    bool found_body = false;
    for (auto it = range.first; it != range.second; ++it) {
        found_body = true;
        if (et >= it->second.start_et && et <= it->second.end_et) {
            if (it->second.type == 2) {
                return evaluateType2(it->second, et);
            } else if (it->second.type == 13) {
                return evaluateType13(it->second, et);
            } else {
                throw std::runtime_error("Unsupported SPK segment type: " + std::to_string(it->second.type));
            }
        }
    }
    
    std::string msg = "No SPK data found for body " + std::to_string(target_id) + " at time " + std::to_string(et);
    if (found_body) {
        msg += ". Available ranges: ";
        for (auto it = range.first; it != range.second; ++it) {
             msg += "[" + std::to_string(it->second.start_et) + " : " + std::to_string(it->second.end_et) + "] ";
        }
    } else {
        msg += ". Body ID not present in file.";
    }
    throw std::runtime_error(msg);
}

Eigen::VectorXd SPKReader::evaluateType2(const SPKSegment& seg, double et) {
    // Type 2: Chebyshev for Position. Velocity via differentiation.
    
    // 1. Locate record
    // Index = floor((et - init) / intlen)
    int rec_idx = (int)std::floor((et - seg.init_sec) / seg.intlen);
    
    // Address of this record (1-based double index)
    // start_addr is start of DATA.
    int rec_start_addr = seg.start_addr + rec_idx * seg.rsize;
    
    // Read record data (Midpoint, Radius, Coeffs...)
    std::vector<double> buf(seg.rsize);
    // Seek: (addr - 1) * 8
    file_.seekg((rec_start_addr - 1) * 8, std::ios::beg);
    
    std::vector<char> byte_buf(seg.rsize * 8);
    file_.read(byte_buf.data(), seg.rsize * 8);
    std::memcpy(buf.data(), byte_buf.data(), seg.rsize * 8);
    
    if (big_endian_ != isSystemBigEndian()) {
        for (double& v : buf) v = swapEndian(v);
    }
    
    double mid = buf[0];
    double rad = buf[1];
    
    // Normalized time x in [-1, 1]
    double x = (et - mid) / rad;
    
    // Evaluate Chebyshev polynomials
    int n_coeffs = (seg.rsize - 2) / 3; // Number of coefficients per component
    
    std::vector<double> T(n_coeffs); // Tk(x)
    std::vector<double> U(n_coeffs); // T'k(x)
    
    T[0] = 1.0;
    T[1] = x;
    U[0] = 0.0;
    U[1] = 1.0;
    
    for (int k = 2; k < n_coeffs; ++k) {
        T[k] = 2.0 * x * T[k-1] - T[k-2];
        U[k] = 2.0 * x * U[k-1] - U[k-2] + 2.0 * T[k-1];
    }
    
    // Sum coefficients
    // Layout: [Mid, Rad, X_coeffs..., Y_coeffs..., Z_coeffs...]
    double pos[3] = {0, 0, 0};
    double vel[3] = {0, 0, 0};
    
    for (int comp = 0; comp < 3; ++comp) {
        int offset = 2 + comp * n_coeffs;
        for (int k = 0; k < n_coeffs; ++k) {
            double c = buf[offset + k];
            pos[comp] += c * T[k];
            vel[comp] += c * U[k];
        }
    }
    
    // Velocity scaling: dx/dt = 1/rad
    for (int i=0; i<3; ++i) vel[i] /= rad;
    
    Eigen::VectorXd state(6);
    state << pos[0], pos[1], pos[2], vel[0], vel[1], vel[2];
    return state;
}

Eigen::VectorXd SPKReader::evaluateType13(const SPKSegment& seg, double et) {
    // Read Metadata from End
    // Layout: ... [Directory N doubles] [Degree] [N] (Last)
    double meta[2];
    file_.seekg((seg.end_addr - 2) * 8, std::ios::beg);
    char mbuf[16];
    file_.read(mbuf, 16);
    std::memcpy(meta, mbuf, 16);
    if (big_endian_ != isSystemBigEndian()) {
        for (double& v : meta) v = swapEndian(v);
    }
    
    // N (number of records) is last
    // Degree (window size?) is penultimate
    int N = (int)meta[1];
    // int degree = (int)meta[0]; // Not used for Hermite cubic
    
    if (N <= 0) throw std::runtime_error("Invalid Type 13 N: " + std::to_string(N));

    // Calculate RSIZE
    // Total Interval Doubles = (end_addr - start_addr + 1)
    // Size = N * RSIZE + N (Directory) + 2 (Meta)
    // N * RSIZE = Total - N - 2
    int total_doubles = seg.end_addr - seg.start_addr + 1;
    int data_doubles = total_doubles - N - 2;
    int rsize = data_doubles / N;
    
    if (rsize != 6) {
        // Fallback or Error?
        // Type 13 usually is 6 (Pos+Vel). If not 6, we might crash.
        // But proceed assuming first 6 are state.
    }

    // Read Directory (N epochs)
    // Located at [End - 2 - N] to [End - 3]
    std::vector<double> epochs(N);
    // int directory_start_addr = seg.end_addr - 2 - N + 1; // Unused 
    // seg.end_addr is 1-based index of last double.
    // Address of meta[1] (N) is end_addr. (index N-1 in buffer sort of)
    // Seek to (directory_start - 1) * 8
    
    // Optimization: Just seek to End - 2 - N
    // Absolute offset in file doubles: (seg.end_addr - 2 - N)
    file_.seekg((seg.end_addr - 2 - N) * 8, std::ios::beg);
    std::vector<char> dbuf(N * 8);
    file_.read(dbuf.data(), N * 8);
    std::memcpy(epochs.data(), dbuf.data(), N * 8);
    
    if (big_endian_ != isSystemBigEndian()) {
        for (double& v : epochs) v = swapEndian(v);
    }
    
    // Find Interval t1 <= et <= t2
    // std::upper_bound finds first element > et
    auto it = std::upper_bound(epochs.begin(), epochs.end(), et);
    
    int idx = std::distance(epochs.begin(), it) - 1;
    // idx points to t1 such that t1 <= et. t2 is at idx+1.
    
    // Boundary checks
    if (idx < 0) idx = 0;
    if (idx >= N - 1) idx = N - 2;
    
    double t1 = epochs[idx];
    double t2 = epochs[idx+1];
    
    // Read State 1 and State 2
    // Data starts at seg.start_addr
    // Record i is at start + i * rsize
    int rec1_addr = seg.start_addr + idx * rsize;
    
    file_.seekg((rec1_addr - 1) * 8, std::ios::beg);
    char sbuf[12 * 8]; // Read 2 records ? RSIZE might be 6.
    // if records are adjacent? Yes.
    file_.read(sbuf, 2 * rsize * 8);
    
    std::vector<double> states(2 * rsize);
    std::memcpy(states.data(), sbuf, 2 * rsize * 8);
    if (big_endian_ != isSystemBigEndian()) {
        for (double& v : states) v = swapEndian(v);
    }
    
    // P1, V1 from rec1
    // P2, V2 from rec2 (which is states[rsize...])
    
    Eigen::Vector3d p1(states[0], states[1], states[2]);
    Eigen::Vector3d v1(states[3], states[4], states[5]);
    
    Eigen::Vector3d p2(states[rsize], states[rsize+1], states[rsize+2]);
    Eigen::Vector3d v2(states[rsize+3], states[rsize+4], states[rsize+5]);
    
    // Cubic Hermite Interpolation
    double dt = t2 - t1;
    double tau = (et - t1) / dt;
    double tau2 = tau * tau;
    double tau3 = tau2 * tau;
    
    // Hermite Basis Functions
    // H00 = 2t^3 - 3t^2 + 1
    // H10 = t^3 - 2t^2 + t
    // H01 = -2t^3 + 3t^2
    // H11 = t^3 - t^2
    
    double h00 = 2*tau3 - 3*tau2 + 1;
    double h10 = tau3 - 2*tau2 + tau;
    double h01 = -2*tau3 + 3*tau2;
    double h11 = tau3 - tau2;
    
    // Position
    Eigen::Vector3d pos = h00 * p1 + h10 * dt * v1 + h01 * p2 + h11 * dt * v2;
    
    // Velocity (Derivatives of Basis * 1/dt)
    // dH00 = 6t^2 - 6t
    // dH10 = 3t^2 - 4t + 1
    // dH01 = -6t^2 + 6t
    // dH11 = 3t^2 - 2t
    
    double dh00 = 6*tau2 - 6*tau;
    double dh10 = 3*tau2 - 4*tau + 1;
    double dh01 = -6*tau2 + 6*tau;
    double dh11 = 3*tau2 - 2*tau;
    
    Eigen::Vector3d vel = (dh00 * p1 + dh10 * dt * v1 + dh01 * p2 + dh11 * dt * v2) / dt;
    
    Eigen::VectorXd res(6);
    res << pos, vel;
    return res;
}

} // namespace astdyn::io
