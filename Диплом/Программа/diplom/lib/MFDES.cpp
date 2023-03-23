#include "MFDES.h"

MFDES::MFDES(const FDESbase& _p) : FDESbase(_p) {}

void MFDES::solve() {

    for (ull i = 0; i <= n; i++) {
        result[0][i] = psi(x_i(i));
    }

    Matrix A(n + 1);
    
    for (ull t = 1; t <= k; t++) {
        
        // create d-vector
        std::vector<double> d(n + 1, 0.0);

        for (ull i = 1; i < n; i++) {
            for (ull j = 1; j <= t; j++) {
                d[i] += g(gamma, j) * result[t-j][i];
            }
            d[i] -= std::pow(tau, gamma) * f(x_i(i), t_k(t));
        }

        //// borders
        d[0] = phiL(t_k(t));
        d[n] = phiR(t_k(t));

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

        //// borders 3Т2П
        for (ull i = 0; i <= n; i++) {
            A[0][i] = 0.0;
            A[n][i] = 0.0;
        }

        A[0][0] = beta_l - 1.5 * alpha_l / h;
        A[0][1] = 2.0 * alpha_l / h;
        A[0][2] = -0.5 * alpha_l / h;
        A[n][n-2] = 0.5 * alpha_r / h;
        A[n][n-1] = -2.0 * alpha_r / h;
        A[n][n] = beta_r + 1.5 * alpha_r / h;

        // solve system
        auto r = A.Solve(d);
        
        for (ull i = 0; i < r.size(); i++) {
            result[t][i] = r[i];
        }
    }
}