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

class MFDES : public FDESbase {

public: 

    MFDES(const FDESbase& _p) : FDESbase(_p) {}

    void solve() {

        for (ull i = 0; i <= n; i++) {
            result[0][i] = psi(x_i(i));
        }

        Matrix A(n + 1);
        
        for (ull t = 1; t <= k; t++) {
            
            // create d-vector
            std::vector<double> d(n + 1, 0.0);

            for (ull i = 0; i <= n; i++) {
                for (ull j = 1; j <= t; j++) {
                    d[i] += g(gamma, j) * result[t-j][i];
                }
                d[i] -= std::pow(tau, gamma) * f(x_i(i), t_k(t));
            }

            // create matrix A
            for (ull i = 0; i <= n; i++) {
                for (ull j = 0; j <= n; j++) {
                    if (i == j) {
                        A[i][j] = a(x_i(i), t_k(t)) * g(alpha, 1) + b(x_i(i), t_k(t)) * g(alpha, 1) - 1.0;
                    } else if (i + 1 < j) {
                        A[i][j] = b(x_i(i), t_k(t)) * g(alpha, j - i + 1);
                    } else if (i > j + 1) {
                        A[i][j] = a(x_i(i), t_k(t)) * g(alpha, i - j + 1);
                    } else if (i == j + 1) {
                        A[i][j] = a(x_i(i), t_k(t)) * g(alpha, 2) + b(x_i(i), t_k(t)) - c(x_i(i), t_k(t));
                    } else {
                        A[i][j] = a(x_i(i), t_k(t)) + b(x_i(i), t_k(t)) * g(alpha, 2) + c(x_i(i), t_k(t));
                    } 
                }
            }

            // borders
            if (borders) {
                d[0] = phiL(t_k(t));
                d[n] = phiR(t_k(t));
                for (ull i = 0; i <= n; i++) {
                    A[0][i] = (i == 0 ? 1.0 : 0.0);
                    A[n][i] = (i == n ? 1.0 : 0.0);
                }
            }

            // solve system
            auto r = A.Solve(d);
            
            for (ull i = 0; i < r.size(); i++) {
                result[t][i] = r[i];
            }

            if (borders) {
                result[t][0] = phiL(t_k(t));
                result[t][r.size() - 1] = phiR(t_k(t));
            }
        }
    }
};