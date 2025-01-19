#pragma once
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

Eigen::Vector3d twoBodyAccel(const Eigen::Vector3d &r, const double mu);
std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> genTwoBodyModel(const double mu);

Eigen::Vector<double, 6> GaussPtb(const Eigen::Vector<double, 6> &El, const Eigen::Vector3d &fi, const double mu);
Eigen::Vector<double, 6> GaussPtb(const Eigen::Vector<double, 6> &El, const Eigen::Vector3d &fi);
std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> genGaussPtbModel();

double timeit(std::function<void()> func);