/**
 * @file Logger.hpp
 * @brief Simple logging system for OrbFit
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 */

#ifndef ASTDYN_UTILS_LOGGER_HPP
#define ASTDYN_UTILS_LOGGER_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>

namespace astdyn {
namespace utils {

enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR,
    FATAL
};

class Logger {
public:
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }
    
    void setLogLevel(LogLevel level) { min_level_ = level; }
    void setLogFile(const std::string& filename) {
        if (file_.is_open()) {
            file_.close();
        }
        file_.open(filename, std::ios::app);
        use_file_ = file_.is_open();
    }
    
    void log(LogLevel level, const std::string& message) {
        if (level < min_level_) return;
        
        std::string level_str;
        switch (level) {
            case LogLevel::DEBUG:   level_str = "DEBUG"; break;
            case LogLevel::INFO:    level_str = "INFO"; break;
            case LogLevel::WARNING: level_str = "WARNING"; break;
            case LogLevel::ERROR:   level_str = "ERROR"; break;
            case LogLevel::FATAL:   level_str = "FATAL"; break;
        }
        
        std::ostringstream oss;
        oss << "[" << getTimestamp() << "] " 
            << level_str << ": " << message;
        
        std::string formatted = oss.str();
        
        std::cout << formatted << std::endl;
        
        if (use_file_ && file_.is_open()) {
            file_ << formatted << std::endl;
            file_.flush();
        }
    }
    
    void debug(const std::string& msg) { log(LogLevel::DEBUG, msg); }
    void info(const std::string& msg) { log(LogLevel::INFO, msg); }
    void warning(const std::string& msg) { log(LogLevel::WARNING, msg); }
    void error(const std::string& msg) { log(LogLevel::ERROR, msg); }
    void fatal(const std::string& msg) { log(LogLevel::FATAL, msg); }
    
private:
    Logger() : min_level_(LogLevel::INFO), use_file_(false) {}
    ~Logger() { if (file_.is_open()) file_.close(); }
    
    std::string getTimestamp() {
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        std::ostringstream oss;
        oss << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S");
        return oss.str();
    }
    
    LogLevel min_level_;
    std::ofstream file_;
    bool use_file_;
};

// Global logging macros
#define LOG_DEBUG(msg) astdyn::utils::Logger::getInstance().debug(msg)
#define LOG_INFO(msg) astdyn::utils::Logger::getInstance().info(msg)
#define LOG_WARNING(msg) astdyn::utils::Logger::getInstance().warning(msg)
#define LOG_ERROR(msg) astdyn::utils::Logger::getInstance().error(msg)
#define LOG_FATAL(msg) astdyn::utils::Logger::getInstance().fatal(msg)

} // namespace utils
} // namespace astdyn

#endif // ASTDYN_UTILS_LOGGER_HPP
