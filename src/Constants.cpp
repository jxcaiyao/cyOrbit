#include <iostream>
#include <sstream>
#include <unistd.h>
#include "Constants.h"
int Constants::importEGM2008(std::string filepath, const uint64_t degree, const uint64_t order, Eigen::MatrixXd &C, Eigen::MatrixXd &S)
{
    std::ifstream file(filepath);
    if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        return -1;
    }

    C = Eigen::MatrixXd::Zero(degree + 1, order + 1);
    S = Eigen::MatrixXd::Zero(degree + 1, order + 1);

    std::string line;
    while (std::getline(file, line))
    {
        std::replace(line.begin(), line.end(), 'D', 'E');

        std::istringstream iss(line);
        uint64_t n, m;
        double Cnm, Snm;
        double Cnm_err, Snm_err;

        if (!(iss >> n >> m >> Cnm >> Snm >> Cnm_err >> Snm_err))
        {
            std::cerr << "Failed to parse line: " << line << std::endl;
            return -1;
        }

        if (n <= degree && m <= order)
        {
            C(n, m) = Cnm;
            S(n, m) = Snm;
        }
        else if (n > degree)
        {
            break;
        }
    }

    file.close();

    return 0;
}

int Constants::importEGM2008(const uint64_t degree, const uint64_t order, Eigen::MatrixXd &C, Eigen::MatrixXd &S)
{
    assert(degree <= 120 && order <= 120);
    std::string relativePath = "./data/EGM2008_to120_TideFree";

    char cwd[256];
    if (getcwd(cwd, sizeof(cwd)) != nullptr)
    {
        std::cout << "Current directory: " << cwd << std::endl;
    }
    else
    {
        perror("getcwd");
    }

    return Constants::importEGM2008(relativePath, degree, order, C, S);
}
