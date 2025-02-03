#include <iostream>
#include <iomanip>
#include <cyOrbit.h>
#include <Constants.h>

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    double mu = Constants::mu_E;
    auto twoBodyModel = genTwoBodyModel(mu);

    Vector<double, 6> El0(42166000, 0.0001, 0.5, 1, 2, 1.5);
    // Vector<double, 6> El0(6378000, 0.0001, 0.5, 1, 2, 1.5);
    VectorXd RV0 = El2RV(El0, mu);

    MatrixXd RVt = ode4(twoBodyModel, RV0, VectorXd::LinSpaced(86400 / 60 + 1, 0, 86400));
    // MatrixXd RVt = ode4(twoBodyModel, RV0, VectorXd::LinSpaced(100 / 10 + 1, 0, 100));

    MatrixXd Elt = RV2El(RVt.transpose()).transpose();

    auto func = [twoBodyModel, El0]()
    {
        VectorXd RV0 = El2RV(El0);
        MatrixXd RVt = ode4(twoBodyModel, RV0, VectorXd::LinSpaced(86400 / 60 + 1, 0, 86400));
    };

    cout << "timeit: " << timeit(func) << endl;

    auto gaussPtbModel = genGaussPtbModel();
    MatrixXd Elt2 = ode4(gaussPtbModel, El0, VectorXd::LinSpaced(86400 / 60 + 1, 0, 86400));
    MatrixXd RVt2 = El2RV(Elt2.transpose()).transpose();

    auto func1 = [gaussPtbModel, El0]()
    {
        MatrixXd Elt2 = ode4(gaussPtbModel, El0, VectorXd::LinSpaced(86400 / 60 + 1, 0, 86400));
        // MatrixXd RVt2 = El2RV(Elt2.transpose()).transpose();
        // VectorXd El1 = RV2El(RV0);
    };

    cout << "timeit: " << timeit(func1) << endl;

    cout << "RVf_err = " << RVt.row(RVt.rows() - 1) - RVt2.row(RVt2.rows() - 1) << endl;
    cout << "Elf_err = " << Elt.row(Elt.rows() - 1) - Elt2.row(Elt2.rows() - 1) << endl;

    auto centralBodyModel = genCentralBodyAccelModel(21, 21, DJ00);

    RV0 << -1697165.15298332, -6564028.937726, -1157311.63506205, 6707.50456351712, -1132.47915606244, -3417.69654930092;
    VectorXd af = centralBodyModel(0, RV0);
    cout << "RV0 = " << RV0.format(IOFormat(FullPrecision, DontAlignCols, ", ", ", ", "", "", "[", "]")) << endl;
    cout << "af = " << af.format(IOFormat(FullPrecision, DontAlignCols, ", ", ", ", "", "", "[", "]")) << endl;

    auto func2 = [centralBodyModel, RV0]()
    {
        VectorXd af = centralBodyModel(0, RV0);
    };
    cout << "timeit: " << timeit(func2) << endl;

    cout << setprecision(15) << "JD0 = " << mjulianSecond(2025, 1, 23, 8, 52, 34.678) << endl;

    double JD0 = 0.0;
    double mJD = 0.0;
    julianDate(2024, 6, 13, 4, 0, 0, JD0, mJD);

    int year, month, day, hour, minute;
    double second;
    JD2Cal(JD0, mJD, year, month, day, hour, minute, second);

    cout << year << " " << month << " " << day << " " << hour << " " << minute << " " << second << endl;

    auto Afi = eci2ecef(JD0, mJD);

    auto func3 = [JD0, mJD]()
    {
        eci2ecef(JD0, mJD);
    };

    cout << "timeit: " << timeit(func3) << endl;

    cout << "Afi = " << Afi << endl;

    Vector3d R0(33068056.167, -26162779.601, -77426.421);
    Vector3d Rf = Afi * R0.head(3);

    cout << "Rf = " << Rf.format(IOFormat(FullPrecision, DontAlignCols, ", ", ", ", "", "", "[", "]")) << endl;

    // ode4 nonsphere test

    // VectorXd x0(6);
    // x0 << 33068056.167, -26162779.601, -77426.421, 1907.672875, 2411.188377, -4.602025;
    El0 << 7000000, 0.001, 50.0 / 180 * M_PI, 2.0 / 180 * M_PI, 3.0 / 180 * M_PI, 6.0 / 180 * M_PI;
    VectorXd x0 = El2RV(El0);
    cout << "x0 = " << x0 << endl;
    centralBodyModel = genCentralBodyAccelModel(21, 21, JD0 + mJD);

    VectorXd tspan = VectorXd::LinSpaced(86400.0 / 60 + 1, 0, 86400);
    RVt2 = ode4(centralBodyModel, x0, tspan);

    // 将tspan和RVt2拼接成7列矩阵
    MatrixXd tRV(RVt2.rows(), 7);
    tRV.col(0) = tspan;
    tRV.rightCols(6) = RVt2;

    // 输出到CSV文件
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
    ofstream file("orbit_data.csv");
    file << tRV.format(CSVFormat);
    file.close();

    cout << "RVt2 = " << tRV << endl;

    auto func4 = [centralBodyModel, x0, tspan]()
    {
        MatrixXd RVt2 = ode4(centralBodyModel, x0, tspan);
    };
    cout << "timeit: " << timeit(func4) << endl;

    return 0;
}
