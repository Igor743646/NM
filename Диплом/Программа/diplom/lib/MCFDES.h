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

#define GPU_ENABLE

class MCFDES : public FDESbase {

    std::vector<double> probabilities;

public: 

    MCFDES(const FDESbase&);

    void solve(ull);

    void print_probs();

#ifdef GPU_ENABLE

    void solve_GPU(ull);

#endif

private:

    std::vector<double> _make_prob();

    void _make_prefsum_prob(std::vector<double>&);

    double drand(double, double);
};