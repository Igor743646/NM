#pragma once

namespace {
    using ull = unsigned long long;
}

struct Problem {

    double L, R;
    double T;
    double alpha, gamma;
    double h, tau;
    double beta;

    std::function<double(double, double)> D;
    std::function<double(double)> psi;
    std::function<double(double, double)> f;
    std::function<double(double)> phiL;
    std::function<double(double)> phiR;

    ull n, k;

    Problem(double _L, double _R, double _T, double _alpha, double _gamma, double _h, double _tau, double _beta,
        std::function<double(double, double)> _D, std::function<double(double)> _psi,
        std::function<double(double, double)> _f, std::function<double(double)> _phiL,
        std::function<double(double)> _phiR) 
    {
        L = _L; R = _R; T = _T;
        alpha = _alpha; gamma = _gamma;
        h = _h; tau = _tau; beta = _beta;
        D = _D;
        psi = _psi;
        f = _f;
        phiL = _phiL;
        phiR = _phiR;

        n = (ull)((R - L) / h);
        k = (ull)(T / tau);
    }

    Problem(const Problem& p) 
    {
        L = p.L; R = p.R; T = p.T;
        alpha = p.alpha; gamma = p.gamma;
        h = p.h; tau = p.tau; beta = p.beta;
        D = p.D;
        psi = p.psi;
        f = p.f;
        phiL = p.phiL;
        phiR = p.phiR;

        n = p.n;
        k = p.k;
    }

};