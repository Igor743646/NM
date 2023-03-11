#pragma once
#include <cmath>

namespace {
    using ull = unsigned long long;
}

struct Problem {

    double L, R;
    double T;
    double alpha, gamma;
    double h, tau;
    double beta;

    std::function<double(double, double)> D;    // коэффициент диффузии при дробной производной по пространству
    std::function<double(double, double)> V;    // коэффициент сноса при производной первой степени
    std::function<double(double)> psi;          // начальное условие при t = 0, u(x, 0) = psi(x)
    std::function<double(double, double)> f;    // функция источник
    std::function<double(double)> phiL;         // граничное условие u(L, t) = phiL(t)
    std::function<double(double)> phiR;         // граничное условие u(R, t) = phiR(t)

    ull n, k;
    bool borders;

    Problem(double _L, double _R, double _T, double _alpha, double _gamma, double _h, double _tau, double _beta,
        std::function<double(double, double)> _D, std::function<double(double, double)> _V, std::function<double(double)> _psi,
        std::function<double(double, double)> _f, std::function<double(double)> _phiL,
        std::function<double(double)> _phiR) 
    {
        L = _L; R = _R; T = _T;
        alpha = _alpha; gamma = _gamma;
        h = _h; tau = _tau; beta = _beta;
        D = _D;
        V = _V;
        psi = _psi;
        f = _f;
        phiL = _phiL;
        phiR = _phiR;

        if (D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha) > gamma / alpha) {
            std::cout << "May be problem with condition\n";
            std::cout << D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha) << " " << gamma / alpha << std::endl;
        }

        n = (ull)((R - L) / h);
        k = (ull)(T / tau);
        borders = true;
    }

    Problem(double _L, double _R, double _T, double _alpha, double _gamma, double _h, double _tau, double _beta,
        std::function<double(double, double)> _D, std::function<double(double, double)> _V, std::function<double(double)> _psi,
        std::function<double(double, double)> _f) 
    {
        L = _L; R = _R; T = _T;
        alpha = _alpha; gamma = _gamma;
        h = _h; tau = _tau; beta = _beta;
        D = _D;
        V = _V;
        psi = _psi;
        f = _f;
        phiL = [](double t){ return 0.0; };
        phiR = [](double t){ return 0.0; };

        if (D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha) > gamma / alpha) {
            std::cout << "May be problem with condition\n";
            std::cout << D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha) << " " << gamma / alpha << std::endl;
        }

        n = (ull)((R - L) / h);
        k = (ull)(T / tau);
        borders = false;
    }

    Problem(const Problem& p) 
    {
        L = p.L; R = p.R; T = p.T;
        alpha = p.alpha; gamma = p.gamma;
        h = p.h; tau = p.tau; beta = p.beta;
        D = p.D;
        V = p.V;
        psi = p.psi;
        f = p.f;
        phiL = p.phiL;
        phiR = p.phiR;

        if (D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha) > gamma / alpha) {
            std::cout << "May be problem with condition\n";
            std::cout << D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha) << " " << gamma / alpha << std::endl;
        }

        n = p.n;
        k = p.k;
        borders = p.borders;
    }

    friend std::ostream& operator<<(std::ostream& out, const Problem& p) {
        return out << p.k + 1 << std::endl << p.h << std::endl << p.tau << std::endl << p.L;
    }

};