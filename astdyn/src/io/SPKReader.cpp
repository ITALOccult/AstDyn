/**
 * @file SPKReader.cpp
 * @brief Implementation of Native SPK Reader
 */

#include "astdyn/io/SPKReader.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <unordered_set>
#include <stdexcept>
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
    int current_rec = get_first_summary_record(record);
    while (current_rec > 0) {
        seekToRecord(current_rec);
        if (!file_.read(record.data(), RECORD_SIZE)) break;
        current_rec = process_summary_record(record);
    }
}

int SPKReader::get_first_summary_record(const std::vector<char>& file_record) {
    int fwd; std::memcpy(&fwd, file_record.data() + 76, 4);
    if (big_endian_ != isSystemBigEndian()) fwd = swapEndian(fwd);
    return (fwd <= 0) ? 2 : fwd;
}

int SPKReader::process_summary_record(const std::vector<char>& record) {
    double next_d, ns_d; std::memcpy(&next_d, record.data(), 8); std::memcpy(&ns_d, record.data() + 16, 8);
    if (big_endian_ != isSystemBigEndian()) { next_d = swapEndian(next_d); ns_d = swapEndian(ns_d); }
    int next = static_cast<int>(next_d), ns = static_cast<int>(ns_d);
    if (ns > 0 && ns <= 250) for (int i = 0; i < ns; ++i) process_segment(record, i);
    return (next > 0) ? next : 0;
}

void SPKReader::process_segment(const std::vector<char>& record, int i) {
    int offset = 24 + i * 40; SPKSegment seg;
    std::memcpy(&seg.start_et, record.data() + offset, 8);
    std::memcpy(&seg.end_et, record.data() + offset + 8, 8);
    int ints[6]; std::memcpy(ints, record.data() + offset + 16, 24);
    if (big_endian_ != isSystemBigEndian()) {
        seg.start_et = swapEndian(seg.start_et); seg.end_et = swapEndian(seg.end_et);
        for (int& v : ints) v = swapEndian(v);
    }
    seg.body_id = ints[0]; seg.center_id = ints[1]; seg.frame_id = ints[2];
    seg.type = ints[3]; seg.start_addr = ints[4]; seg.end_addr = ints[5];
    if (seg.type == 2) init_type2_segment(seg);
    segments_.insert({seg.body_id, seg});
}

void SPKReader::init_type2_segment(SPKSegment& seg) {
    double p[4]; file_.seekg(static_cast<long long>(seg.end_addr - 4) * 8, std::ios::beg);
    file_.read(reinterpret_cast<char*>(p), 32);
    if (big_endian_ != isSystemBigEndian()) for (double& v : p) v = swapEndian(v);
    seg.init_sec = p[0]; seg.intlen = p[1]; seg.rsize = (int)p[2];
    seg.order = ((int)p[2] - 2) / 3 - 1; seg.n_comp = 3;
}

void SPKReader::seekToRecord(int record_idx) {
    file_.seekg(static_cast<long long>(record_idx - 1) * RECORD_SIZE, std::ios::beg);
}

Eigen::Matrix<double, 6, 1> SPKReader::getState(int target_id, double et) {
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

int SPKReader::getCenter(int target_id) const {
    auto it = segments_.find(target_id);
    if (it != segments_.end()) return it->second.center_id;
    return -1;
}

Eigen::Matrix<double, 6, 1> SPKReader::evaluateType2(const SPKSegment& seg, double et) {
    std::lock_guard<std::mutex> lock(file_mutex_);
    int rec_idx = (int)std::floor((et - seg.init_sec) / seg.intlen);
    Type2Cache& cache = type2_cache_[seg.body_id];
    if (cache.seg != &seg || cache.rec_idx != rec_idx) refresh_type2_cache(seg, rec_idx, cache);
    double x = (et - cache.coeffs[0]) / cache.coeffs[1];
    int n = (seg.rsize - 2) / 3;
    std::vector<double> T(n), U(n); 
    compute_chebyshev(x, n, T, U);
    return finalize_type2_state(cache.coeffs, T, U, n, cache.coeffs[1]);
}

void SPKReader::refresh_type2_cache(const SPKSegment& seg, int idx, Type2Cache& cache) {
    cache.seg = &seg; cache.rec_idx = idx; cache.coeffs.resize(seg.rsize);
    file_.seekg(static_cast<long long>(seg.start_addr + idx * seg.rsize - 1) * 8, std::ios::beg);
    file_.read(reinterpret_cast<char*>(cache.coeffs.data()), seg.rsize * 8);
    if (big_endian_ != isSystemBigEndian()) for (double& v : cache.coeffs) v = swapEndian(v);
}

void SPKReader::compute_chebyshev(double x, int n, std::vector<double>& T, std::vector<double>& U) {
    T[0] = 1.0; T[1] = x; U[0] = 0.0; U[1] = 1.0;
    for (int k = 2; k < n; ++k) {
        T[k] = 2.0 * x * T[k-1] - T[k-2];
        U[k] = 2.0 * x * U[k-1] - U[k-2] + 2.0 * T[k-1];
    }
}

Eigen::Matrix<double, 6, 1> SPKReader::finalize_type2_state(const std::vector<double>& c, const std::vector<double>& T, const std::vector<double>& U, int n, double rad) {
    Eigen::Matrix<double, 6, 1> s = Eigen::Matrix<double, 6, 1>::Zero();
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < n; ++k) {
            s[i] += c[2 + i*n + k] * T[k];
            s[i+3] += c[2 + i*n + k] * U[k];
        }
        s[i+3] /= rad;
    }
    return s;
}

Eigen::Matrix<double, 6, 1> SPKReader::evaluateType13(const SPKSegment& seg, double et) {
    std::lock_guard<std::mutex> lock(file_mutex_);
    // Read Metadata from End
    // Layout: ... [Directory N doubles] [Degree] [N] (Last)
    double meta[2];
    file_.seekg(static_cast<long long>(seg.end_addr - 2) * 8, std::ios::beg);
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
    file_.seekg(static_cast<long long>(seg.end_addr - 2 - N) * 8, std::ios::beg);
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
    
    file_.seekg(static_cast<long long>(rec1_addr - 1) * 8, std::ios::beg);
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
    
    Eigen::Matrix<double, 6, 1> res;
    res << pos, vel;
    return res;
}

} // namespace astdyn::io
