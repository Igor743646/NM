#include <iostream>
#include <cmath>

#include "lib/MFDES.h"
#include "lib/MCFDES.h"
#include <fstream>
#include <chrono>

bool ill(double x, double a, double b) { // include_line_left
    return (a <= x && x < b) ? true : false;
}

bool ila(double x, double a, double b) { // include_line_all
    return (a <= x && x <= b) ? true : false;
}

int main() {

    // std::cout << std::fixed;
    // std::cout.precision(3);

    /* u(x, t) = x^2 * t^2 */
    // FDESbase p1(
    //     0.0, 1.0, 0.01, 1.7, 0.7, 0.05, 0.0001, 1.0,
    //     [](double x, double t){ return GAMMA(3.0 - 1.7) / GAMMA(3.0 - 0.7) ; },
    //     [](double x, double t){ return 0.0; },
    //     [](double x){ return 0.0; },
    //     [](double x, double t){ return GAMMA(3.0) / GAMMA(3.0 - 0.7) * (std::pow(t, 2.0 - 0.7) * std::pow(x, 2.0)) - std::pow(x, 2.0 - 1.7) * std::pow(t, 2.0); },
    //     [](double t){ return 0.0; },
    //     [](double t){ return std::pow(t, 2.0); }
    // );

    /* u(x, t) = exp(2t) * x^2 */
    // FDESbase p1(
    //     0, 10.0, 0.2, 1.8, 0.9, 0.2, 0.01, 1.0,
    //     [](double x, double t){ return std::pow(2.0, 0.9) * GAMMA(3.0 - 1.8) / GAMMA(3.0); },                                // D
    //     [](double x, double t){ return 0.0; },                                 // V
    //     [](double x){ return std::pow(x, 2.0); },     // psi(x)
    //     [](double x, double t){ return std::pow(2.0, 0.9) * std::exp(2.0 * t) * (std::pow(x, 2.0) - std::pow(x, 0.2)); },                                  // f(x, t)
    //     [](double t){ return 0.0; },                                            // phiL(t)
    //     [](double t){ return std::exp(2.0 * t) * 100.0; }                                             // phiR(t)
    // );



    /* u(x, t) = ??? */ // alpha = 2.0, gamma = 1.0
    // FDESbase p1(
    //     -0.1, 0.1, 5.0, 2.0, 1.0, 0.01, 0.05, 0.0,
    //     [](double x, double t){ return 0.0005; },                                // D
    //     [](double x, double t){ return 0.0; },                                 // V
    //     [](double x){ return -0.01 <= std::abs(x) && std::abs(x) < 0.0099999 ? 10.0 : 0.0; },     // psi(x)
    //     [](double x, double t){ return 0.0; },                                  // f(x, t)
    //     [](double t){ return 0.1 * exp(-(0.01) / (4 * 0.0005 * t)) / sqrt(std::PI * 4 * 0.0005 * t); },                                            // phiL(t)
    //     [](double t){ return 0.1 * exp(-(0.01) / (4 * 0.0005 * t)) / sqrt(std::PI * 4 * 0.0005 * t); }                                             // phiR(t)
    // );

    /* u(x, t) = ??? */ // alpha = 1.8, gamma = 1.0
    // FDESbase p1(
    //     -0.1, 0.1, 0.01, 1.8, 1.0, 0.001, 0.0001, 0.0,
    //     [](double x, double t){ return 0.0005; },                                // D
    //     [](double x, double t){ return 0.0; },                                 // V
    //     [](double x){ return -0.01 <= std::abs(x) && std::abs(x) < 0.0099999 ? 10.0 : 0.0; },     // psi(x)
    //     [](double x, double t){ return 0.0; },                                  // f(x, t)
    //     [](double t){ return 0.0; },                                            // phiL(t)
    //     [](double t){ return 0.0; }                                             // phiR(t)
    // );

    /* u(x, t) = ??? */ // alpha = 2.0, gamma = 0.9
    // FDESbase p1(
    //     -0.1, 0.1, 5.0, 2.0, 0.9, 0.01, 0.05, 0.0,
    //     [](double x, double t){ return 0.0005; },                                // D
    //     [](double x, double t){ return 0.0; },                                 // V
    //     [](double x){ return -0.01 <= std::abs(x) && std::abs(x) < 0.0099999 ? 10.0 : 0.0; },     // psi(x)
    //     [](double x, double t){ return 0.0; },                                  // f(x, t)
    //     [](double t){ return 0.0; },                                            // phiL(t)
    //     [](double t){ return 0.0; }                                             // phiR(t)
    // );

    /* u(x, t) = exp(2t) * sin(x) */ // alpha = 1.7, gamma = 0.8
    // FDESbase p1(
    //     -std::PI, std::PI, 1.0, 1.7, 0.8, 0.2, 0.005, 0.0,
    //     [](double x, double t){ return std::pow(2.0, 0.8); },                                // D
    //     [](double x, double t){ return 1.0; },                                // V
    //     [](double x){ return std::sin(x); },     // psi(x)
    //     [](double x, double t){ return -2.0 * std::pow(2.0, 0.8) * std::exp(2.0 * t) * 
    //                             std::sin(1.7 / 4.0 * std::PI) * std::cos(x + 1.7 / 4.0 * std::PI); },                      // f(x, t)
    //     [](double t){ return 0.0; },                                            // phiL(t)
    //     [](double t){ return 0.0; }                                             // phiR(t)
    // );

    /* u(x, t) = t^2 */
    // FDESbase p1(
    //     0.0, 1.0, 0.01, 1.7, 0.7, 0.05, 0.0001, 0.0,
    //     [](double x, double t){ return 0.0; },
    //     [](double x, double t){ return 0.0; },
    //     [](double x){ return 0.0; },
    //     [](double x, double t){ return GAMMA(3.0) / GAMMA(3.0 - 0.7) * std::pow(t, 2.0 - 0.7); },
    //     [](double t){ return t * t; },
    //     [](double t){ return t * t; }
    // );

    /* u(x, t) = ??? */ // alpha = 1.7, gamma = 0.8
    FDESbase p1(
        0.0, 2.0, 1.0, 1.7, 0.8, 0.025, 0.005, 0.0,
        [](double x, double t){ return 0.05; },                                // D
        [](double x, double t){ return -1.0; },                                // V
        [](double x){ return ill(x, 0.0, 1.0) ? 4.0 : (ila(x, 1.0, 2.0) ? 2.0 : 0.0); },     // psi(x)
        [](double x, double t){ return 0.0; },                      // f(x, t)
        [](double t){ return 1.0; },                                            // phiL(t)
        [](double t){ return 1.0; }                                             // phiR(t)
    );

    p1.Info();

    MFDES solution1(p1);
    MCFDES solution2(p1);

    solution2.print_probs();

    std::chrono::time_point<std::chrono::system_clock> start_time, end_time;

    start_time = std::chrono::system_clock::now();
    solution1.solve();
    end_time = std::chrono::system_clock::now();

    std::cout << "MFDES (by matrix solver) time: " << std::chrono::duration<double>(end_time - start_time).count() << std::endl;

    start_time = std::chrono::system_clock::now();
    solution2.solve_GPU(4000);
    end_time = std::chrono::system_clock::now();

    std::cout << "MCFDES (by Monte-Carlo algo) time: " << std::chrono::duration<double>(end_time - start_time).count() << std::endl;

    /* FILE WRITING */

    std::fstream file;
    file.open("out.txt", std::ios_base::out);

    file << solution1 << std::endl;

    file << solution2 << std::endl;

    file.close();
    
    return 0;
}