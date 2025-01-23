#pragma once

#include <Eigen/Dense>
#include <fstream>

class Constants
{
public:
    static constexpr double mu_E = 3.986004418e14; // m^3/s^2 Earth's gravitational parameter

public:
    static constexpr double GM_EGM2008 = 3.986004418e14; // m^3/s^2 Earth's gravitational parameter
    static constexpr double Re_EGM2008 = 6378137;
    static int importEGM2008(std::string filepath, const uint64_t degree, const uint64_t order, Eigen::MatrixXd &C, Eigen::MatrixXd &S);
    static int importEGM2008(const uint64_t degree, const uint64_t order, Eigen::MatrixXd &C, Eigen::MatrixXd &S);
};
