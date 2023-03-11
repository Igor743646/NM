#include <iostream>
#include <cmath>

#include "utils.h"
#include "matrix2D.h"
#include "MFDES.h"
#include "MCFDES.h"
#include <fstream>
#include <chrono>

using ull = unsigned long long;

int main() {

    // std::cout << std::fixed;

    // std::cout.precision(20);

    // Problem p1(
    //     0.0, 1.0, 1.0, 1.7, 0.7, 0.05, 0.01, 0.5,
    //     [](double x, double t){ return GAMMA(3.0 - 1.7) / GAMMA(3.0) * (std::pow(x, 1.7)); },
    //     [](double x, double t){ return 0.0; },
    //     [](double x){ return 0.0; },
    //     [](double x, double t){ return GAMMA(3.0) / GAMMA(3.0 - 0.7) * (std::pow(t, 2.0 - 0.7) * std::pow(x, 2.0)) - std::pow(x*t, 2.0); },
    //     [](double t){ return 0.0; },
    //     [](double t){ return std::pow(t, 2.0); }
    // );

    Problem p1(
        -0.1, 0.1, 0.5, 1.7, 0.9, 0.01, 0.01, 0.0,
        [](double x, double t){ return 0.005; },                                // D
        [](double x, double t){ return 0.05; },                                 // V
        [](double x){ return (-0.005 <= x && x < 0.00499) ? 100.0 : 0.0; },     // psi(x)
        // [](double x){ return -(x-1)*(x+1); },
        [](double x, double t){ return std::pow(x, 2.0) * t; },                                  // f(x, t)
        [](double t){ return t; },                                            // phiL(t)
        [](double t){ return t; }                                             // phiR(t)
    );

    // std::cout << p2.D(0, 0) * std::pow(p2.tau, p2.gamma) / std::pow(p2.h, p2.alpha) << " " << p2.gamma / p2.alpha << std::endl;

    // MFDES solution1(p1);
    MFDES solution1(p1);
    MCFDES solution2(p1);

    std::chrono::time_point<std::chrono::system_clock> start_time, end_time;

    start_time = std::chrono::system_clock::now();
    auto result1 = solution1.solve();
    end_time = std::chrono::system_clock::now();

    std::cout << "MFDES (by matrix solver) time: " << std::chrono::duration<double>(end_time - start_time).count() << std::endl;

    start_time = std::chrono::system_clock::now();
    auto result2 = solution2.solve(5000);
    end_time = std::chrono::system_clock::now();

    std::cout << "MCFDES (by Monte-Carlo algo) time: " << std::chrono::duration<double>(end_time - start_time).count() << std::endl;

    /* FILE WRITING */

    std::fstream file;
    file.open("out.txt", std::ios_base::out);

    file << p1 << std::endl;
    for (ull t = 0; t < result1.size(); t++) {
        for (ull i = 0; i < result1[t].size(); i++) {
            file << result1[t][i] << " ";
        }
        file << std::endl;
    }

    file << p1 << std::endl;
    for (ull t = 0; t < result2.size(); t++) {
        for (ull i = 0; i < result2[t].size(); i++) {
            file << result2[t][i] << " ";
        }
        file << std::endl;
    }

    // file << p1 << std::endl;
    // for (ull t = 0; t < result1.size(); t++) {
    //     for (ull i = 0; i < result1[t].size(); i++) {
    //         file << result1[t][i] << " ";
    //     }
    //     file << std::endl;
    // }

    file.close();
    
    return 0;
}