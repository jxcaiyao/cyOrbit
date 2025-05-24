#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <string>
#include "EOPDataHandler.cpp"

// 测试用例：加载数据
TEST(EOPDataHandlerTest, LoadEOPData)
{
    std::string testFilename = "data/EOP-Last5Years.txt"; // 直接使用现有文件

    EOPDataHandler handler;
    handler.loadEOPData(); // 不需要输入参数

    std::string expectedData = "2020 01 01 58849  0.076614  0.282309 -0.1771665  0.0004417 -0.108563 -0.006596  0.000358 -0.000007  37";
    std::string dateKey = "2020-1-1";
    EXPECT_EQ(handler.eopDataMap[dateKey], expectedData);
}

// 测试用例：查找数据
TEST(EOPDataHandlerTest, FindEOPDataByDate)
{
    std::string testFilename = "data/EOP-Last5Years.txt"; // 直接使用现有文件

    EOPDataHandler handler;
    handler.loadEOPData(); // 不需要输入参数

    std::vector<std::string> result = handler.findEOPDataByDate(2020, 1, 1);
    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(result[0], "2020 01 01 58849  0.076614  0.282309 -0.1771665  0.0004417 -0.108563 -0.006596  0.000358 -0.000007  37");
}

// 测试用例：插值数据
TEST(EOPDataHandlerTest, InterpolateEOPDataByDateTime)
{
    std::string testFilename = "data/EOP-Last5Years.txt"; // 直接使用现有文件

    EOPDataHandler handler;
    handler.loadEOPData(); // 不需要输入参数

    std::vector<double> interpolatedValues = handler.interpolateEOPDataByDateTime(2020, 1, 1, 12, 0, 0.0);
    EXPECT_EQ(interpolatedValues.size(), 13);
    EXPECT_NEAR(interpolatedValues[2], 0.075650, 1e-6);
    EXPECT_NEAR(interpolatedValues[3], 0.2825015, 1e-6);
    EXPECT_NEAR(interpolatedValues[4], -0.17739965, 1e-6);
    EXPECT_NEAR(interpolatedValues[5], 0.00046225, 1e-6);
    EXPECT_NEAR(interpolatedValues[6], -0.1085655, 1e-6);
    EXPECT_NEAR(interpolatedValues[7], -0.0066245, 1e-6);
    EXPECT_NEAR(interpolatedValues[8], 0.0003775, 1e-6);
    EXPECT_NEAR(interpolatedValues[9], 0.0000085, 1e-6);
    EXPECT_NEAR(interpolatedValues[10], 0.0003275, 1e-6);
    EXPECT_NEAR(interpolatedValues[11], 0.0000085, 1e-6);
    EXPECT_NEAR(interpolatedValues[12], 37.0, 1e-6);
}
