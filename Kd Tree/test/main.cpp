#include "Utils.hpp"
#include <gtest/gtest.h>

int main(int argc, char* argv[]) {
    CS350::ChangeWorkdir();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
