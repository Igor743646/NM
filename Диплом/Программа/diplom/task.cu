#include <iostream>
#include <cmath>

#include "lib/MFDES.h"
#include "lib/MCFDES.h"
#include <fstream>
#include <chrono>

int main() {


    /* u(x, t) = exp(2t) * sin(x) */ // alpha = 1.7, gamma = 0.8
    FDESbase p1(
        -0.5, 0.495, 0.999, 2.0, 0.9, 0.005, 0.001, 0.0,
        [](double x){ return 0.00068; },                                // D
        [](double x){ return 0; },                                // V
        [](double x){ return -0.005 <= std::abs(x) && std::abs(x) < 0.0049999 ? 1.0 : 0.0; },     // psi(x)
        [](double x, double t){ return 0.0; },                      // f(x, t)
        [](double t){ return 0.0; },                                            // phiL(t)
        [](double t){ return 0.0; },false                                             // phiR(t)
    );

    p1.Info();

    MFDES solution1(p1);

    std::chrono::time_point<std::chrono::system_clock> start_time, end_time;

    start_time = std::chrono::system_clock::now();
    solution1.solve();
    end_time = std::chrono::system_clock::now();

    std::cout << "MFDES (by matrix solver) time: " << std::chrono::duration<double>(end_time - start_time).count() << std::endl;

    /* FILE WRITING */

    std::fstream file;
    file.open("out.txt", std::ios_base::out);

    file << solution1 << std::endl;

    file.close();
    
    return 0;
}