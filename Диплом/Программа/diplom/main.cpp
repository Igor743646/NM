#include <iostream>
#include <cmath>

#include "utils.h"
#include "matrix2D.h"
#include "MFDES.h"

using ull = unsigned long long;

int main() {

    Problem p1(
        0.0, 1.0, 2.0, 1.7, 0.7, 0.05, 0.05, 0.5,
        [](double x, double t){ return GAMMA(3.0 - 1.7) / GAMMA(3.0) * (std::pow(x, 1.7)); },
        [](double x){ return 0.0; },
        [](double x, double t){ return GAMMA(3.0) / GAMMA(3.0 - 0.7) * (std::pow(t, 2.0 - 0.7) * std::pow(x, 2.0)) - std::pow(x*t, 2.0); },
        [](double t){ return 0.0; },
        [](double t){ return std::pow(t, 2.0); }
    );

    MFDES solution1(p1);

    auto result = solution1.solve();

    for (ull t = 0; t < result.size(); t++) {
        for (ull i = 0; i < result[t].size(); i++) {
            std::cout << result[t][i] << " ";
        }
        std::cout << std::endl;
    }
}