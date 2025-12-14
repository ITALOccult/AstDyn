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
constexpr int RECORD_SIZE = 1024;
constexpr int DBL_SIZE = 8;

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
    
    // std::cout << "[SPKReader] FWD Ptr from Header: " << fwd << "\n";
    
    // Fallback: IF fwd is 0 or weird, start at record 2 (standard for most .bsp)
    if (fwd <= 0) fwd = 2;
    
    int current_rec = fwd;
    int total_segments = 0;
    
    while (current_rec > 0) {
        seekToRecord(current_rec);
        file_.read(record.data(), RECORD_SIZE);
        
        double next_d_ptr, prev_d_ptr, ns_d;
        std::memcpy(&next_d_ptr, record.data(), 8);
        std::memcpy(&ns_d, record.data() + 16, 8);
        
        if (big_endian_ != isSystemBigEndian()) {
            next_d_ptr = swapEndian(next_d_ptr);
            ns_d = swapEndian(ns_d);
        }
        
        int next = (int)next_d_ptr;
        int ns = (int)ns_d;
        
        // Safety break
        if (ns > 250 || ns < 0) break; // Garbage read
        
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
                int end_byte = (seg.end_addr - 1) * 8;
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
            
            segments_.insert({seg.body_id, seg});
            total_segments++;
        }
        
        current_rec = next;
    }
    std::cout << "[SPKReader] Loaded " << total_segments << " segments.\n";
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

} // namespace astdyn::io
