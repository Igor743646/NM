/*
    Monte Carlo Fractional Differential Equations Solver 
*/

#pragma once

#include <functional>
#include <algorithm>
#include <numeric>
#include <time.h>

#include "matrix2D.h"
#include "FDESbase.h"
#include "utils.h"

namespace {
    using ull = unsigned long long;

    template<class T>
    std::vector<T> prefix_sum(const std::vector<T>& vec) {
        std::vector<T> result(vec.size(), vec[0]);

        for (ull i = 1; i < vec.size(); i++) {
            result[i] = vec[i] + result[i - 1];
        }

        return result;
    }
}

class MCFDES {

    Problem p;

    void _make_prefsum_prob(std::vector<double>& prefsum_prob) {
        std::vector<double> probabilities(2 * p.n + 2 + p.k , 0.0);

        probabilities[p.n - 1] = a(0.0, 0.0) + b(0.0, 0.0) * g(p.alpha, 2.0) - c(0.0, 0.0);
        probabilities[p.n] = p.gamma - p.alpha * (a(0.0, 0.0) + b(0.0, 0.0));
        probabilities[p.n + 1] = a(0.0, 0.0) * g(p.alpha, 2.0) + b(0.0, 0.0) + c(0.0, 0.0);

        for (ull i = 2; i <= p.n; i++) {
            probabilities[p.n + i] = a(0.0, 0.0) * g(p.alpha, i + 1);
            probabilities[p.n - i] = b(0.0, 0.0) * g(p.alpha, i + 1);
        }

        for (ull i = 1; i <= p.k; i++) {
            // std::cout << -g(p.gamma, i + 1) << std::endl;
            probabilities[2 * p.n + i] = -g(p.gamma, i + 1);
        }

        probabilities[2 * p.n + 1 + p.k] = 0.0;
        probabilities[2 * p.n + 1 + p.k] = 1.0 - std::accumulate(probabilities.begin(), probabilities.end(), 0.0);

        prefsum_prob = prefix_sum(probabilities);

        print(probabilities);
    }

    double _make_participle_count_in_cage(std::vector<ull>& participle_count_in_cage, ull count) {
        double square = 0.0;

        for (ull i = 0; i <= p.n; i++) {
            square += p.psi(x_i(i) - p.h / 2.0) + p.psi(x_i(i) + p.h / 2.0);
        }

        // print(square);

        for (ull i = 0; i <= p.n; i++) {
            participle_count_in_cage[i] = (ull)std::ceil((double)count * (p.psi(x_i(i) - p.h / 2.0) + p.psi(x_i(i) + p.h / 2.0)) / square);
        }

        // print(participle_count_in_cage);

        return square;
    }

    long long _start_x_position(std::vector<ull>& participle_count_in_cage) {
        long long result = 0;

        while (0 == participle_count_in_cage[result]) {
            result++;
        }

        if (result < participle_count_in_cage.size()) {
            participle_count_in_cage[result]--; 
        } else {
            result--;
        }
            

        // print(participle_count_in_cage);

        return result;
    }

public: 

    MCFDES(const Problem& _p) : p(_p) {}

    std::vector<std::vector<double>> solve(ull count = 1000) {

        std::vector<std::vector<double>> result(p.k + 1, std::vector<double>(p.n + 1, 0.0));

        std::vector<double> prefsum_prob;

        _make_prefsum_prob(prefsum_prob);

        // print(prefsum_prob);

        // Учёт начального и граничных условий

        for (ull i = 0; i <= p.n; i++) {
            result[0][i] = p.psi(x_i(i));
        }

        for (ull j = 1; j <= p.k; j++) {
            result[j][0] = p.phiL(t_k(j));
            result[j][p.n] = p.phiR(t_k(j));
        }

        // Симуляция

        for (ull i = 1; i < p.n; i++) {
            for (ull j = 1; j <= p.k; j++) {
                for (ull _n = 0; _n < count; _n++) {
                    long long x = i, y = j;

                    while (y > 0 && x < p.n && x > 0) {
                        double rnd = drand(0.0, 1.0);
                        long long idx = 0;

                        while ((idx < 2 * p.n + 2 + p.k) && (rnd > prefsum_prob[idx])) {
                            idx++;
                        }

                        result[j][i] += p.f(x_i(x), t_k(y)) * std::pow(p.tau, p.gamma);

                        if (idx <= 2 * p.n) { // перемещение по пространству
                            x += idx - p.n;
                            y--;
                        } else if (idx <= 2 * p.n + p.k) { // перемещение по времени
                            y -= idx - 2 * p.n + 1;
                        }

                        // if (y > 0 && x < p.n && x > 0)
                        // result[j][i] += p.f(x_i(x), t_k(y)) * std::pow(p.tau, p.gamma);

                        //result[j][i] += p.f(x_i(x), t_k(y)) * std::pow(p.tau, p.gamma);
                        // if(_n == count - 1 && i == p.n / 2 && j == p.k / 2) {
                        //     std::cout << "(" << x << ", " << y << ")\n";
                        // }
                    }

                    if (y == 0) {
                        result[j][i] += p.psi(x_i(x));
                    } else if (x == 0) {
                        result[j][i] += p.phiL(t_k(y));
                    } else if (x == p.n) {
                        result[j][i] += p.phiR(t_k(y));
                    }
                }

                result[j][i] /= (double)count;
            }
        }

        return result;
    }

    inline double drand(double dmin, double dmax) {
        return dmin + (double)rand() / RAND_MAX * (dmax - dmin);
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