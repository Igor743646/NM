#include "MFDES.h"

MFDES::MFDES(const FDESbase& _p) : FDESbase(_p) {}

void MFDES::solve() {

    for (ull i = 0; i <= n; i++) {
        result[0][i] = psi(x_i(i));
    }

    Matrix A(n + 1);

    // create matrix A
    for (ull i = 0; i <= n; i++) {
        for (ull j = 0; j <= n; j++) {
            if (i == j) {
                A[i][j] = a(x_i(i)) * g(alpha, 1) + b(x_i(i)) * g(alpha, 1) - 1.0;
            } else if (i + 1 < j) {
                A[i][j] = b(x_i(i)) * g(alpha, j - i + 1);
            } else if (i > j + 1) {
                A[i][j] = a(x_i(i)) * g(alpha, i - j + 1);
            } else if (i == j + 1) {
                A[i][j] = a(x_i(i)) * g(alpha, 2) + b(x_i(i)) - c(x_i(i));
            } else {
                A[i][j] = a(x_i(i)) + b(x_i(i)) * g(alpha, 2) + c(x_i(i));
            } 
        }
    }

    //// 1-border conditions
    if (is_borders) {
        for (ull i = 0; i <= n; i++) {
            A[0][i] = 0.0;
            A[n][i] = 0.0;
        }

        A[0][0] = 1.0;
        A[n][n] = 1.0;
    }
    
    for (ull t = 1; t <= k; t++) {
        
        // create d-vector
        std::vector<double> d(n + 1, 0.0);

        for (ull i = 0; i <= n; i++) {
            for (ull j = 1; j <= t; j++) {
                d[i] += g(gamma, j) * (result[t-j][i] - result[0][i]);
            }
            d[i] -= result[0][i];
            d[i] -= std::pow(tau, gamma) * f(x_i(i), t_k(t));
        }

        //// borders
        if (is_borders) {
            d[0] = phiL(t_k(t));
            d[n] = phiR(t_k(t));
        }
        
        // solve system
        auto r = A.Solve(d);
        
        for (ull i = 0; i < r.size(); i++) {
            result[t][i] = r[i];
        }
    }
}