#pragma once
#include <cmath>
#include "utils.h"

namespace {
    using ull = unsigned long long;
}

class FDESbase {

protected:

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

public:

    std::vector<std::vector<double>> result;

    FDESbase(double _L, double _R, double _T, double _alpha, double _gamma, double _h, double _tau, double _beta,
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
        n = (ull)((R - L) / h);
        k = (ull)(T / tau);
        borders = true;
        result = std::vector<std::vector<double>>(k + 1, std::vector<double>(n + 1, 0.0));
    }

    FDESbase(double _L, double _R, double _T, double _alpha, double _gamma, double _h, double _tau, double _beta,
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
        n = (ull)((R - L) / h);
        k = (ull)(T / tau);
        borders = false;
        result = std::vector<std::vector<double>>(k + 1, std::vector<double>(n + 1, 0.0));
    }

    FDESbase(const FDESbase& p) 
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
        n = p.n;
        k = p.k;
        borders = p.borders;
        result = p.result;
    }

    void Info() {

        print("Info:");
        print("\tx : [", L, "; ", R, "], t : [", 0, "; ", T, "]");
        print("\th = ", h, "; tau = ", tau);
        print("\tN = ", n);
        print("\tT = ", k);

        print("\tD * (tau^gamma) / (h^alpha) = ", D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha));
        print("\tgamma / alpha = ", gamma / alpha);

        if (D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha) > gamma / alpha) {
            print("May be problem with condition");
            print(D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha), " ", gamma / alpha);
        }
        print();
    }

    inline double g(double a, ull j) {
        if (std::abs(a) - (ull)std::abs(a) < 0.000001 && j >= a + 1) {
            return 0.0;
        }
        return (j % 2 == 0 ? 1 : -1) / (a + 1) / BETA((double)j + 1.0, a - (double)j + 1.0);
    }

    inline double x_i(ull i) {
        return L + i * h;
    }

    inline double t_k(ull k) {
        return k * tau;
    }

    inline double a(double x, double t) {
        return (1.0 + beta) * (D(x, t) / 2.0) * (std::pow(tau, gamma) / std::pow(h, alpha));
    }

    inline double b(double x, double t) {
        return (1.0 - beta) * (D(x, t) / 2.0) * (std::pow(tau, gamma) / std::pow(h, alpha));
    }

    inline double c(double x, double t) {
        return V(x, t) * std::pow(tau, gamma) / 2 / h;
    }

    friend std::ostream& operator<<(std::ostream& out, const FDESbase& p) {
        out << p.k + 1 << std::endl << p.h << std::endl << p.tau << std::endl << p.L << std::endl;
        for (ull t = 0; t < p.result.size(); t++) {
            for (ull i = 0; i < p.result[t].size(); i++) {
                out << p.result[t][i] << " ";
            }
            if (t < p.result.size() - 1)
                out << std::endl;
        }
        return out;
    }

};