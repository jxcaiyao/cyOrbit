#pragma once
#include <functional>
#include <Eigen/Dense>
#include <sofa.h>
#include <sofam.h>

Eigen::MatrixXd ode4(std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> func, Eigen::VectorXd x0, Eigen::VectorXd tspan);
Eigen::VectorXd ode4Iter(std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> func, Eigen::VectorXd x0, double dt);

Eigen::MatrixXd El2RV(const double mu, Eigen::MatrixXd El);
Eigen::MatrixXd RV2El(const double mu, Eigen::MatrixXd RV);

Eigen::Vector3d twoBodyAccel(const Eigen::Vector3d &r, const double mu);

std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> genTwoBodyModel(const double mu);
