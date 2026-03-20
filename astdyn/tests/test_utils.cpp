/**
 * @file test_utils.cpp
 * @brief Unit tests for StringUtils and Logger
 */

#include <gtest/gtest.h>
#include "astdyn/utils/StringUtils.hpp"
#include "astdyn/utils/Logger.hpp"

using namespace astdyn::utils;

// ========== StringUtils Tests ==========

TEST(StringUtilsTest, Trim) {
    EXPECT_EQ(trim("  hello  "), "hello");
    EXPECT_EQ(trim("\t\ntest\n\t"), "test");
    EXPECT_EQ(trim("no-spaces"), "no-spaces");
    EXPECT_EQ(trim(""), "");
    EXPECT_EQ(trim("   "), "");
}

TEST(StringUtilsTest, CaseConversion) {
    EXPECT_EQ(to_lower("HELLO"), "hello");
    EXPECT_EQ(to_lower("MiXeD"), "mixed");
    EXPECT_EQ(to_upper("hello"), "HELLO");
    EXPECT_EQ(to_upper("MiXeD"), "MIXED");
}

TEST(StringUtilsTest, Split) {
    auto tokens = split("a,b,c", ',');
    ASSERT_EQ(tokens.size(), 3);
    EXPECT_EQ(tokens[0], "a");
    EXPECT_EQ(tokens[1], "b");
    EXPECT_EQ(tokens[2], "c");
    
    auto tokens2 = split("one:two:three", ':');
    ASSERT_EQ(tokens2.size(), 3);
    EXPECT_EQ(tokens2[0], "one");
    EXPECT_EQ(tokens2[1], "two");
    EXPECT_EQ(tokens2[2], "three");
}

TEST(StringUtilsTest, IsNumeric) {
    EXPECT_TRUE(is_numeric("12345"));
    EXPECT_TRUE(is_numeric("0"));
    EXPECT_FALSE(is_numeric("123.45"));
    EXPECT_FALSE(is_numeric("abc"));
    EXPECT_FALSE(is_numeric(""));
}

TEST(StringUtilsTest, IsAlpha) {
    EXPECT_TRUE(is_alpha("hello"));
    EXPECT_TRUE(is_alpha("ABC"));
    EXPECT_FALSE(is_alpha("hello123"));
    EXPECT_FALSE(is_alpha("123"));
    EXPECT_FALSE(is_alpha(""));
}

TEST(StringUtilsTest, ReplaceAll) {
    EXPECT_EQ(replace_all("hello world", "world", "there"), "hello there");
    EXPECT_EQ(replace_all("aaa", "a", "b"), "bbb");
    EXPECT_EQ(replace_all("test", "x", "y"), "test");
}

TEST(StringUtilsTest, StartsWithEndsWith) {
    EXPECT_TRUE(starts_with("hello world", "hello"));
    EXPECT_FALSE(starts_with("hello world", "world"));
    EXPECT_TRUE(ends_with("hello world", "world"));
    EXPECT_FALSE(ends_with("hello world", "hello"));
}

TEST(StringUtilsTest, Padding) {
    EXPECT_EQ(pad_left("5", 3, '0'), "005");
    EXPECT_EQ(pad_right("5", 3, '0'), "500");
    EXPECT_EQ(pad_left("test", 2), "test"); // Already longer
}

// ========== Logger Tests ==========

TEST(LoggerTest, Singleton) {
    auto& logger1 = Logger::getInstance();
    auto& logger2 = Logger::getInstance();
    
    EXPECT_EQ(&logger1, &logger2);
}

TEST(LoggerTest, LogLevels) {
    auto& logger = Logger::getInstance();
    
    // Should not throw
    EXPECT_NO_THROW(logger.debug("Debug message"));
    EXPECT_NO_THROW(logger.info("Info message"));
    EXPECT_NO_THROW(logger.warning("Warning message"));
    EXPECT_NO_THROW(logger.error("Error message"));
}

TEST(LoggerTest, SetLogLevel) {
    auto& logger = Logger::getInstance();
    
    // Should not throw
    EXPECT_NO_THROW(logger.setLogLevel(LogLevel::WARNING));
    EXPECT_NO_THROW(logger.setLogLevel(LogLevel::DEBUG));
}

TEST(LoggerTest, Macros) {
    // Should not throw
    EXPECT_NO_THROW(LOG_DEBUG("Debug via macro"));
    EXPECT_NO_THROW(LOG_INFO("Info via macro"));
    EXPECT_NO_THROW(LOG_WARNING("Warning via macro"));
    EXPECT_NO_THROW(LOG_ERROR("Error via macro"));
}

