#include <cyOrbit.h>

Eigen::Matrix3d rotx(double angle)
{
    // vector rotation matrix about x-axis
    Eigen::Matrix3d R;
    R << 1, 0, 0,
        0, cos(angle), -sin(angle),
        0, sin(angle), cos(angle);
    return R;
}

Eigen::Matrix3d roty(double angle)
{
    // vector rotation matrix about y-axis
    Eigen::Matrix3d R;
    R << cos(angle), 0, sin(angle),
        0, 1, 0,
        -sin(angle), 0, cos(angle);
    return Eigen::Matrix3d();
}

Eigen::Matrix3d rotz(double angle)
{
    // vector rotation matrix about z-axis
    Eigen::Matrix3d R;
    R << cos(angle), -sin(angle), 0,
        sin(angle), cos(angle), 0,
        0, 0, 1;
    return R;
}

Eigen::MatrixXd ode4(std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> func, const Eigen::VectorXd &x0, const Eigen::VectorXd &tspan)
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

Eigen::VectorXd ode4Iter(std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> func, const Eigen::VectorXd &x0, double dt)
{
    Eigen::VectorXd k1 = func(0, x0);
    Eigen::VectorXd k2 = func(dt / 2, x0 + dt / 2 * k1);
    Eigen::VectorXd k3 = func(dt / 2, x0 + dt / 2 * k2);
    Eigen::VectorXd k4 = func(dt, x0 + dt * k3);

    Eigen::VectorXd x1 = x0 + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

    return x1;
}

double M2f(double M, double e)
{
    assert(e >= 0 && e < 1);

    // M2E
    // double E = M;
    // for (uint64_t j = 0; j < 5; j++)
    // {
    //     E = M + e * sin(E);
    // }

    // Use Newton's method to solve for E
    double E = M;
    double tol = 1e-15; // tolerance
    for (uint64_t j = 0; j < 10; j++)
    {
        double f_E = E - e * sin(E) - M;
        double f_E_prime = 1 - e * cos(E);
        double dE = f_E / f_E_prime;
        E -= dE;
        if (fabs(dE) < tol)
            break;
    }

    // E2f
    double f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));

    return f;
}

double f2M(double f, double e)
{
    assert(e >= 0 && e < 1);

    double E = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(f / 2));
    double M = E - e * sin(E);
    return M;
}

Eigen::MatrixXd El2RV(const Eigen::MatrixXd &El, const double mu)
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
        double f = M2f(M, e);

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
        double r = p / (1 + e * cf);
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

Eigen::MatrixXd El2RV(const Eigen::MatrixXd &El)
{
    return El2RV(El, Constants::mu_E);
}

Eigen::MatrixXd RV2El(const Eigen::MatrixXd &RV, const double mu)
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
        double M = f2M(f, e);
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

Eigen::MatrixXd RV2El(const Eigen::MatrixXd &RV)
{
    return RV2El(RV, Constants::mu_E);
}

Eigen::Matrix3d El2Aoi(const Eigen::VectorXd &El)
{
    // transform orbit elements [a,e,i,W,w,M] to aoi matrix

    assert(El.size() == 6);

    double e = El(1);
    double i = El(2);
    double W = El(3);
    double w = El(4);
    double M = El(5);
    double f = M2f(M, e);
    Eigen::Matrix3d R = rotz(-(w + f)) * rotx(-i) * rotz(-W);

    return R;
}

int julianDate(const int year, const int month, const int day, const int hour, const int minute, const double second, double &d1, double &d2)
{
    int ret = iauCal2jd(year, month, day, &d1, &d2);
    d2 += (hour + (minute + second / 60) / 60) / 24;
    return ret;
}

double julianDate(const int year, const int month, const int day, const int hour, const int minute, const double second)
{
    double d1 = 0.0;
    double d2 = 0.0;

    julianDate(year, month, day, hour, minute, second, d1, d2);

    return d1 + d2;
}

int julianDate(const Eigen::VectorXd &cal, double &d1, double &d2)
{
    assert(cal.size() == 6);

    int year = int(cal(0));
    int month = int(cal(1));
    int day = int(cal(2));
    int hour = int(cal(3));
    int minute = int(cal(4));
    double second = cal(5);

    return julianDate(year, month, day, hour, minute, second, d1, d2);
}

double julianDate(const Eigen::VectorXd &cal)
{
    double d1 = 0.0;
    double d2 = 0.0;

    int ret = julianDate(cal, d1, d2);
    if (ret != 0)
    {
        std::cerr << "julianDate error" << std::endl;
        return -1;
    }

    return d1 + d2;
}

double mjulianSecond(const int year, const int month, const int day, const int hour, const int minute, const double second)
{
    // relative to 2451545.0 2000-01-01 12:00:00 UTC

    double d1 = 0.0;
    double d2 = 0.0;
    int ret = julianDate(year, month, day, hour, minute, second, d1, d2);
    return ((d1 - 2451545.0) + d2) * 86400;
}

int JD2Cal(const double dj1, const double dj2, int &year, int &month, int &day, int &hour, int &minute, double &second)
{
    double fd = 0.0;
    int ret = iauJd2cal(dj1, dj2, &year, &month, &day, &fd);

    hour = int(fd * 24);
    minute = int((fd * 24 - hour) * 60);
    second = ((fd * 24 - hour) * 60 - minute) * 60;

    return ret;
}

int JD2Cal(const double JD, int &year, int &month, int &day, int &hour, int &minute, double &second)
{
    double dj1 = 2451545.0;
    double dj2 = JD - dj1;

    return JD2Cal(dj1, dj2, year, month, day, hour, minute, second);
}

Eigen::Matrix3d eci2ecef(const double utc1, const double utc2)
{
    double dut1 = 0.0;

    double uta = 0.0;
    double utb = 0.0;
    iauUtcut1(utc1, utc2, dut1, &uta, &utb);

    double tai1 = 0.0;
    double tai2 = 0.0;
    iauUtctai(uta, utb, &tai1, &tai2);

    double tta = 0.0;
    double ttb = 0.0;
    iauTaitt(tai1, tai2, &tta, &ttb);

    double xp = 0.0;
    double yp = 0.0;
    double rc2t[3][3];
    iauC2t00b(tta, ttb, uta, utb, xp, yp, rc2t);

    Eigen::Matrix3d Afi;
    Afi << rc2t[0][0], rc2t[0][1], rc2t[0][2],
        rc2t[1][0], rc2t[1][1], rc2t[1][2],
        rc2t[2][0], rc2t[2][1], rc2t[2][2];

    return Afi;
}

Eigen::Matrix3d eci2ecef(const double JD0)
{
    double utc1 = 2451545.0;
    double utc2 = JD0 - utc1;

    return eci2ecef(utc1, utc2);
}

Eigen::Vector3d twoBodyAccel(const Eigen::Vector3d &r, const double mu)
{
    Eigen::Vector3d acc;
    acc = -mu * r / pow(r.norm(), 3);
    return acc;
}

std::function<Eigen::VectorXd(double JD0, const Eigen::VectorXd &x0)> genTwoBodyModel(const double mu)
{
    return [mu](double t, const Eigen::VectorXd &x) -> Eigen::VectorXd
    {
        Eigen::VectorXd dxdt(x.size());
        dxdt.segment(0, 3) = x.segment(3, 3);
        dxdt.segment(3, 3) = twoBodyAccel(x.segment(0, 3), mu);
        return dxdt;
    };
}

std::function<Eigen::VectorXd(double JD0, const Eigen::VectorXd &x0)> genTwoBodyModel()
{
    return genTwoBodyModel(Constants::mu_E);
}

Eigen::Vector3d CentralBodyAccel(const Eigen::Vector3d &Rf, const uint64_t degree, const uint64_t order, const Eigen::MatrixXd &C, const Eigen::MatrixXd &S, const double Re, const double GM)
{
    double DU = Re;
    double TU = sqrt(DU * DU * DU / GM);
    double Req = 1;
    double mu = 1;
    Eigen::Vector3d r = Rf / DU;

    /////////////
    Eigen::Vector3d dRdr = Eigen::Vector3d::Zero();
    Eigen::Vector3d dRdphi = Eigen::Vector3d::Zero();
    Eigen::Vector3d dRdlambda = Eigen::Vector3d::Zero();

    double norm_R = r.norm();
    double rho = sqrt(r(0) * r(0) + r(1) * r(1));
    double lambda = atan2(r(1), r(0));
    double phi = asin(r(2) / norm_R);

    double F1 = norm_R * norm_R;

    ////////////
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(degree + 3, degree + 3);
    Eigen::MatrixXd scaleFactor = Eigen::MatrixXd::Zero(degree + 3, degree + 3);

    Associted_Legendre(degree + 2, r, P, scaleFactor);

    Eigen::Vector3d drdR = r / norm_R;

    double sum_dRdr = 0;

    for (uint64_t n = 2; n <= degree; n++)
    {
        double G_r = pow(F1, -(n + 2.0) / 2);
        double sum_dRdr1 = 0;
        for (uint64_t m = 0; m <= order; m++)
        {
            double cos_sin_f = C(n, m) * cos(m * lambda) + S(n, m) * sin(m * lambda);
            sum_dRdr1 += P(n, m) * cos_sin_f;
        }

        sum_dRdr += (n + 1) * pow(Req, n) * G_r * sum_dRdr1;
    }

    dRdr = -mu * sum_dRdr * drdR;

    ////////////

    Eigen::Vector3d dphidR = (1 / rho) * Eigen::Vector3d(-r(0) * r(2) / F1, -r(1) * r(2) / F1, rho * rho / F1);

    double sum_dRdphi = 0;

    for (uint64_t n = 2; n <= degree; n++)
    {
        for (uint64_t m = 0; m <= order; m++)
        {
            double g_n1 = pow(F1, -(n + 1.0) / 2);
            double tan_phi = r(2) / rho;
            sum_dRdphi += pow(Req, n) * g_n1 * (P(n, m + 1) * scaleFactor(n, m) - m * tan_phi * P(n, m)) * (C(n, m) * cos(m * lambda) + S(n, m) * sin(m * lambda));
        }
    }

    dRdphi = mu * sum_dRdphi * dphidR;

    ////////////

    Eigen::Vector3d dlambdadR = (1 / rho / rho) * Eigen::Vector3d(-r(1), r(0), 0);

    double sum_dRdlambda = 0;

    for (uint64_t n = 2; n <= degree; n++)
    {
        double G_lambda = pow(F1, -(n + 1.0) / 2);
        double sum_dRdlambda1 = 0;
        for (uint64_t m = 0; m <= order; m++)
        {
            double del_cos_sin_f = m * (-C(n, m) * sin(m * lambda) + S(n, m) * cos(m * lambda));
            sum_dRdlambda1 += P(n, m) * del_cos_sin_f;
        }

        sum_dRdlambda += pow(Req, n) * G_lambda * sum_dRdlambda1;
    }

    dRdlambda = mu * sum_dRdlambda * dlambdadR;

    ////////////

    double R = Rf.norm();
    Eigen::Vector3d acc = -GM / R / R / R * Rf + (dRdr + dRdphi + dRdlambda) * DU / TU / TU;

    return acc;
}

int Associted_Legendre(const uint64_t N, const Eigen::Vector3d &R, Eigen::MatrixXd &P, Eigen::MatrixXd &scaleFactor)
{
    double x = R(2) / R.norm();
    double y = sqrt(1 - x * x);

    Eigen::MatrixXd psi = Eigen::MatrixXd::Zero(N + 1, N + 1);
    Eigen::MatrixXd kappa = Eigen::MatrixXd::Zero(N + 1, N + 1);
    Eigen::MatrixXd zi = Eigen::MatrixXd::Zero(N + 1, N + 1);

    zi(1, 0) = sqrt(3);
    for (uint64_t n = 2; n <= N; n++)
    {
        zi(n, 0) = sqrt((2 * n + 1.0) / (2 * n));
    }

    for (uint64_t n = 0; n <= N; n++)
    {
        for (uint64_t m = 0; m < n; m++)
        {
            psi(n, m) = sqrt(((2 * n + 1.0) * (2 * n - 1.0) / (n - m) / (n + m)));
        }
    }

    for (uint64_t n = 2; n <= N; n++)
    {
        for (uint64_t m = 0; m < n - 1; m++)
        {
            kappa(n, m) = psi(n, m) / psi(n - 1, m);
        }
    }

    // calculate P
    P = Eigen::MatrixXd::Zero(N + 1, N + 1);
    P(0, 0) = 1;

    for (uint64_t n = 1; n <= N; n++)
    {
        P(n, n) = zi(n, 0) * y * P(n - 1, n - 1);
        P(n, n - 1) = psi(n, n - 1) * x * P(n - 1, n - 1);
    }

    for (uint64_t n = 2; n <= N; n++)
    {
        for (uint64_t m = 0; m < n - 1; m++)
        {
            P(n, m) = psi(n, m) * x * P(n - 1, m) - kappa(n, m) * P(n - 2, m);
        }
    }

    // calculate scaleFactor
    scaleFactor = Eigen::MatrixXd::Zero(N + 1, N + 1);
    scaleFactor(0, 0) = 0;
    scaleFactor(1, 0) = 1;
    scaleFactor(1, 1) = 0;

    for (uint64_t n = 2; n <= N; n++)
    {
        uint64_t k = n;
        for (uint64_t m = 0; m <= n; m++)
        {
            uint64_t p = m;
            if (n == m)
            {
                scaleFactor(k, k) = 0;
            }
            else if (m == 0)
            {
                scaleFactor(k, p) = sqrt((n + 1.0) * n / 2);
            }
            else
            {
                scaleFactor(k, p) = sqrt((n + m + 1.0) * (n - m));
            }
        }
    }

    return 0;
}

std::function<Eigen::VectorXd(double JD0, const Eigen::VectorXd &x0)> genCentralBodyAccelModel(const uint64_t degree, const uint64_t order, const Eigen::MatrixXd &C, const Eigen::MatrixXd &S, const double Re, const double mu)
{
    std::function<Eigen::VectorXd(double JD0, const Eigen::VectorXd &x0)> func = [degree, order, C, S, Re, mu](double JD0, const Eigen::VectorXd &x0) -> Eigen::VectorXd
    {
        Eigen::Vector3d Ri = x0.segment(0, 3);
        Eigen::Vector3d Vi = x0.segment(3, 3);

        Eigen::Matrix3d Afi = eci2ecef(JD0);
        Eigen::Vector3d Rf = Afi * Ri;

        Eigen::Vector3d af = CentralBodyAccel(Rf, degree, order, C, S, Re, mu);

        Eigen::Vector3d ai = Afi.transpose() * af;

        Eigen::VectorXd dxdt(x0.size());
        dxdt.segment(0, 3) = Vi;
        dxdt.segment(3, 3) = ai;

        return dxdt;
    };
    return func;
}

std::function<Eigen::VectorXd(double JD0, const Eigen::VectorXd &x0)> genCentralBodyAccelModel(const uint64_t degree, const uint64_t order)
{
    assert(degree >= 0 && order >= 0);
    assert(degree >= order);

    Eigen::MatrixXd C;
    Eigen::MatrixXd S;

    int ret = Constants::importEGM2008(degree, order, C, S);
    if (ret != 0)
    {
        std::cerr << "Failed to import EGM2008" << std::endl;
        exit(-1);
    }

    return genCentralBodyAccelModel(degree, order, C, S, Constants::Re_EGM2008, Constants::GM_EGM2008);
}

Eigen::Vector<double, 6> GaussPtb(const Eigen::Vector<double, 6> &El, const Eigen::Vector3d &fi, const double mu)
{
    assert(El.size() == 6);
    assert(fi.size() == 3);

    double a = El(0);
    double e = El(1);
    double i = El(2);
    double W = El(3);
    double w = El(4);
    double M = El(5);
    double f = M2f(M, e);
    double u = w + f;

    double ci = cos(i);
    double si = sin(i);
    double cW = cos(W);
    double sW = sin(W);
    double cf = cos(f);
    double sf = sin(f);
    double cu = cos(u);
    double su = sin(u);

    double e2 = e * e;
    double p = a * (1 - e2);
    double r = p / (1 + e * cf);
    double h = sqrt(mu * p);
    double n = sqrt(mu / a / a / a);

    double fix = fi(0);
    double fiy = fi(1);
    double fiz = fi(2);

    double c11 = -cW * su - sW * ci * cu;
    double c12 = cW * ci * cu - sW * su;
    double c13 = cu * si;
    double c21 = -sW * si;
    double c22 = cW * si;
    double c23 = -ci;
    double c31 = sW * ci * su - cW * cu;
    double c32 = -sW * cu - cW * ci * su;
    double c33 = -si * su;

    double fr = -(c31 * fix + c32 * fiy + c33 * fiz);
    double ft = c11 * fix + c12 * fiy + c13 * fiz;
    double fn = -(c21 * fix + c22 * fiy + c23 * fiz);

    double rhft = r * h * ft;

    double dadt = 2 * a * a / h * (e * sf * fr + p / r * fn);
    double dedt = 1 / h * (p * sf * fr + ((p + r) / r * cf + r * e) * fn);
    double didt = cu * r * h * ft;
    double dWdt = su / si * r * h * ft;
    double dwdt = 1 / h / e * (-p * cf * fr + (p + r) * sf * fn) - dWdt * ci;
    double dMdt = n + sqrt(1 - e2) / h / e * ((p * cf - 2 * e * r) * fr - (p + r) * sf * fn);

    Eigen::Vector<double, 6> dEldt;
    dEldt << dadt, dedt, didt, dWdt, dwdt, dMdt;

    return dEldt;
}

Eigen::Vector<double, 6> GaussPtb(const Eigen::Vector<double, 6> &El, const Eigen::Vector3d &fi)
{
    return GaussPtb(El, fi, Constants::mu_E);
}

std::function<Eigen::VectorXd(double t0, const Eigen::VectorXd &x0)> genGaussPtbModel()
{
    return [](double t, const Eigen::VectorXd &x) -> Eigen::VectorXd
    {
        Eigen::VectorXd dxdt(x.size());

        dxdt = GaussPtb(x, Eigen::Vector3d(0, 0, 0));

        return dxdt;
    };
}

double timeit(std::function<void()> func)
{
    // return the calculation time(seconds) of func
    double itermax = 5;
    auto start = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < itermax; i++)
    {
        func();
    }
    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1e9 / itermax;
    return duration;
}
