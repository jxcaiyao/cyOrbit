#include <cyOrbit.h>

Eigen::MatrixXd ode4(std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> func, Eigen::VectorXd x0, Eigen::VectorXd tspan)
{
    Eigen::Index n = tspan.size();
    Eigen::Index m = x0.size();
    Eigen::MatrixXd x(n, m);

    x.row(0) = x0;

    for (Eigen::Index i = 1; i < n; i++)
    {
        x.row(i) = ode4Iter(func, x.row(i - 1), tspan(i) - tspan(i - 1));
    }

    return x;
}

Eigen::VectorXd ode4Iter(std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> func, Eigen::VectorXd x0, double dt)
{
    Eigen::VectorXd k1 = func(0, x0);
    Eigen::VectorXd k2 = func(dt / 2, x0 + dt / 2 * k1);
    Eigen::VectorXd k3 = func(dt / 2, x0 + dt / 2 * k2);
    Eigen::VectorXd k4 = func(dt, x0 + dt * k3);

    x0 = x0 + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

    return x0;
}

Eigen::MatrixXd El2RV(const double mu, Eigen::MatrixXd El)
{
    // transform orbit elements [a,e,i,W,w,M] to position and velocity [rx, ry, rz, vx, vy, vz]
    assert(El.rows() == 6);

    Eigen::MatrixXd RV(6, El.cols());

    for (Eigen::Index n = 0; n < El.cols(); n++)
    {
        double a = El(0, n);
        double e = El(1, n);
        double i = El(2, n);
        double W = El(3, n);
        double w = El(4, n);
        double M = El(5, n);

        // double E = M;
        // for (int j = 0; j < 5; j++)
        // {
        //     E = M + e * sin(E);
        // }

        // Use Newton's method to solve for E
        double E = M;
        double tol = 1e-15; // tolerance
        for (int j = 0; j < 10; j++)
        {
            double f_E = E - e * sin(E) - M;
            double f_E_prime = 1 - e * cos(E);
            double dE = f_E / f_E_prime;
            E -= dE;
            if (fabs(dE) < tol)
                break;
        }

        double f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));

        double cf = cos(f);
        double sf = sin(f);
        double ci = cos(i);
        double si = sin(i);
        double cw = cos(w);
        double sw = sin(w);
        double cW = cos(W);
        double sW = sin(W);

        double c11 = cW * cw - sW * sw * ci;
        double c12 = -cW * sw - sW * cw * ci;
        double c21 = sW * cw + cW * sw * ci;
        double c22 = -sW * sw + cW * cw * ci;
        double c31 = sw * si;
        double c32 = cw * si;

        double p = a * (1 - e * e);
        double r = p / (1 + e * cos(f));
        double rx = r * cf;
        double ry = r * sf;
        double v = sqrt(mu / p);
        double vx = -v * sf;
        double vy = v * (e + cf);

        RV(0, n) = c11 * rx + c12 * ry;
        RV(1, n) = c21 * rx + c22 * ry;
        RV(2, n) = c31 * rx + c32 * ry;
        RV(3, n) = c11 * vx + c12 * vy;
        RV(4, n) = c21 * vx + c22 * vy;
        RV(5, n) = c31 * vx + c32 * vy;
    }

    return RV;
}

Eigen::MatrixXd RV2El(const double mu, Eigen::MatrixXd RV)
{
    // transform position and velocity [rx, ry, rz, vx, vy, vz] to orbit elements [a,e,i,W,w,M]

    assert(RV.rows() == 6);

    Eigen::MatrixXd El(6, RV.cols());

    for (Eigen::Index n = 0; n < RV.cols(); n++)
    {
        double rx = RV(0, n);
        double ry = RV(1, n);
        double rz = RV(2, n);
        double vx = RV(3, n);
        double vy = RV(4, n);
        double vz = RV(5, n);

        Eigen::Vector3d r_vec(rx, ry, rz);
        Eigen::Vector3d v_vec(vx, vy, vz);
        Eigen::Vector3d h_vec = r_vec.cross(v_vec);

        double r = r_vec.norm();
        double v = v_vec.norm();
        double h = h_vec.norm();

        double a = 1 / (2 / r - v * v / mu);
        double e2 = 1 - (h * h) / (a * mu);
        if (e2 < 0 && e2 >= -1e-15)
        {
            e2 = 0;
        }
        double e = sqrt(e2);
        double i = acos(h_vec(2) / h);
        double W = atan2(h_vec(0), -h_vec(1));

        Eigen::Vector3d n_vec = Eigen::Vector3d(-h_vec(1), h_vec(0), 0).normalized();
        Eigen::Vector3d e_vec = -r_vec / r + v_vec.cross(h_vec) / mu;
        Eigen::Vector3d c_vec = h_vec.cross(n_vec) / h;
        double p = h * h / mu;

        double w = atan2(c_vec.dot(e_vec), n_vec.dot(e_vec));
        double f = atan2(r_vec.dot(v_vec) / sqrt(mu * p), p / r - 1);
        double E = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(f / 2));
        double M = E - e * sin(E);
        M = fmod(M, 2 * M_PI);
        if (M < 0)
        {
            M += 2 * M_PI;
        }

        El(0, n) = a;
        El(1, n) = e;
        El(2, n) = i;
        El(3, n) = W;
        El(4, n) = w;
        El(5, n) = M;
    }

    return El;
}

Eigen::Vector3d twoBodyAccel(const Eigen::Vector3d &r, const double mu)
{
    Eigen::Vector3d acc;
    acc = -mu * r / pow(r.norm(), 3);
    return acc;
}

std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> genTwoBodyModel(const double mu)
{
    return [mu](double t, const Eigen::VectorXd &x) -> Eigen::VectorXd
    {
        Eigen::VectorXd dxdt(x.size());
        dxdt.segment(0, 3) = x.segment(3, 3);
        dxdt.segment(3, 3) = twoBodyAccel(x.segment(0, 3), mu);
        return dxdt;
    };
}