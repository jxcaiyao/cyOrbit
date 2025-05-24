#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cmath>

class EOPDataHandler
{
public:
    std::unordered_map<std::string, std::string> eopDataMap;

public:
    void loadEOPData()
    {
        std::ifstream file("./data/EOP-Last5Years.txt");
        std::string line;

        if (!file.is_open())
        {
            std::cerr << "Failed to open EOP-Last5Years.txt" << std::endl;
            return;
        }

        while (std::getline(file, line))
        {
            // 删除结尾的 \r
            if (!line.empty() && line.back() == '\r')
            {
                line.pop_back();
            }

            std::istringstream iss(line);
            int fileYear, fileMonth, fileDay;

            // 跳过注释行和非数据行
            if (line.empty() || line[0] == '#' || line.find("VERSION") != std::string::npos || line.find("UPDATED") != std::string::npos ||
                line.find("BEGIN OBSERVED") != std::string::npos || line.find("END OBSERVED") != std::string::npos ||
                line.find("BEGIN PREDICTED") != std::string::npos || line.find("END PREDICTED") != std::string::npos ||
                line.find("NUM_OBSERVED_POINTS") != std::string::npos || line.find("NUM_PREDICTED_POINTS") != std::string::npos)
            {
                continue;
            }
            // 提取日期信息
            if (iss >> fileYear >> fileMonth >> fileDay)
            {
                std::string dateKey = std::to_string(fileYear) + "-" + std::to_string(fileMonth) + "-" + std::to_string(fileDay);
                eopDataMap[dateKey] = line;
            }
        }

        file.close();
    }

    std::vector<std::string> findEOPDataByDate(int year, int month, int day)
    {
        std::vector<std::string> result;
        std::string dateKey = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);

        if (eopDataMap.find(dateKey) != eopDataMap.end())
        {
            result.push_back(eopDataMap[dateKey]);
        }

        return result;
    }

    std::vector<double> interpolateEOPDataByDateTime(int year, int month, int day, int hour, int minute, double second)
    {
        // 获取当天和下一天的数据
        std::vector<std::string> todayData = findEOPDataByDate(year, month, day);
        std::vector<std::string> tomorrowData = findEOPDataByDate(year, month, day + 1);

        if (todayData.empty() || tomorrowData.empty())
        {
            std::cerr << "Data not found for the given date or next day." << std::endl;
            return {};
        }

        // 解析当天和下一天的数据
        std::istringstream issToday(todayData[0]);
        std::istringstream issTomorrow(tomorrowData[0]);

        std::vector<double> todayValues;
        std::vector<double> tomorrowValues;
        double value;

        while (issToday >> value)
        {
            todayValues.push_back(value);
        }

        while (issTomorrow >> value)
        {
            tomorrowValues.push_back(value);
        }

        // 计算插值因子
        double fractionOfDay = (hour * 3600.0 + minute * 60.0 + second) / (24.0 * 3600.0);

        // 进行线性插值
        std::vector<double> interpolatedValues;
        for (size_t i = 0; i < todayValues.size(); ++i)
        {
            double interpolatedValue = todayValues[i] + fractionOfDay * (tomorrowValues[i] - todayValues[i]);
            interpolatedValues.push_back(interpolatedValue);
        }

        return interpolatedValues;
    }
};
