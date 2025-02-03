#pragma once
#include <iostream>
#include <functional>
#include <Eigen/Dense>
#include <sofa.h>
#include <sofam.h>
#include <Constants.h>
#include <chrono>

Eigen::Matrix3d rotx(double angle);
Eigen::Matrix3d roty(double angle);
Eigen::Matrix3d rotz(double angle);

Eigen::MatrixXd ode4(std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> func, const Eigen::VectorXd &x0, const Eigen::VectorXd &tspan);
Eigen::VectorXd ode4Iter(std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> func, const Eigen::VectorXd &x0, double dt);

double M2f(double M, double e);
double f2M(double f, double e);

Eigen::MatrixXd El2RV(const Eigen::MatrixXd &El, const double mu);
Eigen::MatrixXd El2RV(const Eigen::MatrixXd &El);
Eigen::MatrixXd RV2El(const Eigen::MatrixXd &RV, const double mu);
Eigen::MatrixXd RV2El(const Eigen::MatrixXd &RV);

Eigen::Matrix3d El2Aoi(const Eigen::VectorXd &El);

int julianDate(const int year, const int month, const int day, const int hour, const int minute, const double second, double &d1, double &d2);
double julianDate(const int year, const int month, const int day, const int hour, const int minute, const double second);
int julianDate(const Eigen::VectorXd &cal, double &d1, double &d2);
double julianDate(const Eigen::VectorXd &cal);

double mjulianSecond(const int year, const int month, const int day, const int hour, const int minute, const double second);

int JD2Cal(const double dj1, const double dj2, int &year, int &month, int &day, int &hour, int &minute, double &second);
int JD2Cal(const double JD, int &year, int &month, int &day, int &hour, int &minute, double &second);
Eigen::VectorXd JD2Cal(const double dj1, const double dj2);
Eigen::VectorXd JD2Cal(const double JD);

Eigen::Matrix3d eci2ecef(const double JD0, const double JD1);
Eigen::Matrix3d eci2ecef(const double JD);

Eigen::Vector3d twoBodyAccel(const Eigen::Vector3d &r, const double mu);
std::function<Eigen::VectorXd(double JD0, const Eigen::VectorXd &x0)> genTwoBodyModel(const double mu);
std::function<Eigen::VectorXd(double JD0, const Eigen::VectorXd &x0)> genTwoBodyModel();

Eigen::Vector3d CentralBodyAccel(const Eigen::Vector3d &Rf, const uint64_t degree, const uint64_t order, const Eigen::MatrixXd &C, const Eigen::MatrixXd &S, const double Re, const double GM);
int Associted_Legendre(const uint64_t N, const Eigen::Vector3d &R, Eigen::MatrixXd &P, Eigen::MatrixXd &scaleFactor);
std::function<Eigen::VectorXd(double JD0, const Eigen::VectorXd &x0)> genCentralBodyAccelModel(const uint64_t degree, const uint64_t order, const Eigen::MatrixXd &C, const Eigen::MatrixXd &S, const double Re, const double mu, const double JD0);
std::function<Eigen::VectorXd(double JD0, const Eigen::VectorXd &x0)> genCentralBodyAccelModel(const uint64_t degree, const uint64_t order, const double JD0);

Eigen::Vector3d AtmosphereAccel(const Eigen::Vector3d &R, const Eigen::Vector3d &V, const double CD, const double AMR, const double JD0);
Eigen::Vector3d SolarRadiationAccel(const Eigen::Vector3d &R, const double Cr, const double AMR, const double JD0);
Eigen::Vector3d ThirdBodyGravityAccel(const Eigen::Vector3d &R, const double JD0);

Eigen::Vector<double, 6> GaussPtb(const Eigen::Vector<double, 6> &El, const Eigen::Vector3d &fi, const double mu);
Eigen::Vector<double, 6> GaussPtb(const Eigen::Vector<double, 6> &El, const Eigen::Vector3d &fi);
std::function<Eigen::VectorXd(double JD0, const Eigen::VectorXd &x0)> genGaussPtbModel();

double timeit(std::function<void()> func);