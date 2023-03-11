/*
    Modified Fractional Differential Equations Solver 
*/

#pragma once

#include <functional>

#include "matrix2D.h"
#include "FDESbase.h"
#include "utils.h"

namespace {
    using ull = unsigned long long;
}

class MFDES {

    Problem p;

public: 

    MFDES(const Problem& _p) : p(_p) {}

    std::vector<std::vector<double>> solve() {

        std::vector<std::vector<double>> result(p.k + 1, std::vector<double>(p.n + 1, 0.0));

        for (ull i = 0; i <= p.n; i++) {
            // print(x_i(i));
            result[0][i] = p.psi(x_i(i));
        }

        Matrix A(p.n + 1);
        
        for (ull t = 1; t <= p.k; t++) {
            
            // create d-vector
            std::vector<double> d(p.n + 1, 0.0);

            for (ull i = 0; i <= p.n; i++) {
                for (ull j = 1; j <= t; j++) {
                    d[i] += g(p.gamma, j) * result[t-j][i];
                }
                d[i] -= std::pow(p.tau, p.gamma) * p.f(x_i(i), t_k(t));
            }

            // create matrix A
            for (ull i = 0; i <= p.n; i++) {
                for (ull j = 0; j <= p.n; j++) {
                    if (i == j) {
                        A[i][j] = a(x_i(i), t_k(t)) * g(p.alpha, 1) + b(x_i(i), t_k(t)) * g(p.alpha, 1) - 1.0;
                    } else if (i + 1 < j) {
                        A[i][j] = b(x_i(i), t_k(t)) * g(p.alpha, j - i + 1);
                    } else if (i > j + 1) {
                        A[i][j] = a(x_i(i), t_k(t)) * g(p.alpha, i - j + 1);
                    } else if (i == j + 1) {
                        A[i][j] = a(x_i(i), t_k(t)) * g(p.alpha, 2) + b(x_i(i), t_k(t)) - c(x_i(i), t_k(t));
                    } else {
                        A[i][j] = a(x_i(i), t_k(t)) + b(x_i(i), t_k(t)) * g(p.alpha, 2) + c(x_i(i), t_k(t));
                    } 
                }
            }

            // borders
            if (p.borders) {
                d[0] = p.phiL(t_k(t));
                d[p.n] = p.phiR(t_k(t));
                for (ull i = 0; i <= p.n; i++) {
                    A[0][i] = (i == 0 ? 1.0 : 0.0);
                    A[p.n][i] = (i == p.n ? 1.0 : 0.0);
                }
            }

            // std::cout << A << std::endl;

            // solve system
            auto r = A.Solve(d);

            // print(r);

            for (ull i = 0; i < r.size(); i++) {
                result[t][i] = r[i];
            }

            if (p.borders) {
                result[t][0] = p.phiL(t_k(t));
                result[t][r.size() - 1] = p.phiR(t_k(t));
            }
        }

        return result;
    }

    inline double g(double a, ull j) {
        if (std::abs(a) - (ull)std::abs(a) < 0.0001 && j >= a + 1) {
            return 0.0;
        }
        return (j % 2 == 0 ? 1 : -1) / (a + 1) / BETA((double)j + 1.0, a - (double)j + 1.0);
    }

    inline double x_i(ull i) {
        return p.L + i * p.h;
    }

    inline double t_k(ull k) {
        return k * p.tau;
    }

    inline double a(double x, double t) {
        return (1.0 + p.beta) * (p.D(x, t) / 2.0) * (std::pow(p.tau, p.gamma) / std::pow(p.h, p.alpha));
    }

    inline double b(double x, double t) {
        return (1.0 - p.beta) * (p.D(x, t) / 2.0) * (std::pow(p.tau, p.gamma) / std::pow(p.h, p.alpha));
    }

    inline double c(double x, double t) {
        return p.V(x, t) * std::pow(p.tau, p.gamma) / 2 / p.h;
    }

};