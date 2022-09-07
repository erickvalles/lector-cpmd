#include <gtest/gtest.h>
#include "utils/utils.h"

TEST(utils, evaluaCajaRecursivo){
    GTEST_ASSERT_EQ(evaluaCajaRecursivo(3.5,2.3,1.15),1.2);
}

int main(int argc, char* argv[]){
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}