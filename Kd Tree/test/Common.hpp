#ifndef COMMON_HPP
#define COMMON_HPP
#include "Logging.hpp"
#include <gtest/gtest.h>

// GTEST Utilities
inline char const* TestName() { return ::testing::UnitTest::GetInstance()->current_test_info()->name(); }
inline char const* TestSuiteName() { return ::testing::UnitTest::GetInstance()->current_test_info()->test_suite_name(); }

#endif // COMMON_HPP
