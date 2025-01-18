#include <iostream>
#include <cyOrbit.h>
#include <Constants.h>

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    double mu = Constants::mu_E;
    auto twoBodyModel = genTwoBodyModel(mu);

    Vector<double, 6> El0(42166000, 0.0001, -0.5, -1, -2, -1.5);
    VectorXd RV0 = El2RV(mu, El0);

    // MatrixXd RVt = ode4(twoBodyModel, RV0, VectorXd::LinSpaced(86400 / 60 + 1, 0, 86400));
    MatrixXd RVt = ode4(twoBodyModel, RV0, VectorXd::LinSpaced(100 / 10 + 1, 0, 100));

    cout << "RVt: " << endl;
    cout << RVt << endl;

    MatrixXd Elt = RV2El(mu, RVt.transpose()).transpose();

    cout << "Elt: " << endl;
    cout << Elt << endl;

    return 0;
}