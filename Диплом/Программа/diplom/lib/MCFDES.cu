#include "MCFDES.h"

#ifdef GPU_ENABLE

    #include "MCFDES_GPU.h"

#endif

#define DEFAULT_COUNT 1000
#define TBLOCKS 64
#define THREADS 64

MCFDES::MCFDES(const FDESbase& _p) : FDESbase(_p) {
    probabilities = _make_prob();
}

void MCFDES::solve(ull count = DEFAULT_COUNT) {

    std::vector<double> prefsum_prob;

    _make_prefsum_prob(prefsum_prob);

    // Учёт начального и граничных условий

    for (ull i = 0; i <= n; i++) {
        result[0][i] = psi(x_i(i));
    }

    for (ull j = 1; j <= k; j++) {
        result[j][0] = phiL(t_k(j));
        result[j][n] = phiR(t_k(j));
    }

    // Симуляция

    for (ull i = 1; i < n; i++) {
        for (ull j = 1; j <= k; j++) {
            for (ull _n = 0; _n < count; _n++) {
                long long x = i, y = j;

                while (y > 0 && x < n && x > 0) {
                    double rnd = drand(0.0, 1.0);
                    long long idx = 0;

                    idx = (long long)bin_search(prefsum_prob, rnd);

                    result[j][i] += f(x_i(x), t_k(y)) * std::pow(tau, gamma);

                    if (idx <= 2 * n) { // перемещение по пространству
                        x += idx - n;
                        y--;
                    } else if (idx <= 2 * n + k) { // перемещение по времени
                        y -= idx - 2 * n + 1;
                    } else {
                        break;
                    }
                }

                if (y == 0) {
                    result[j][i] += psi(x_i(x));
                } else if (x == 0) {
                    result[j][i] += phiL(t_k(y));
                } else if (x == n) {
                    result[j][i] += phiR(t_k(y));
                }
            }

            result[j][i] /= (double)count;
        }
    }
}

void MCFDES::print_probs() {
    // N
    for (ull i = 0; i <= 2 * n; i++) {
        std::cout << probabilities[i] << " ";
    }
    std::cout << std::endl;

    // T
    for (ull j = 2 * n + 1; j <= 2 * n + k; j++) {
        std::cout << probabilities[j] << " ";
    }
    std::cout << std::endl;
}

#ifdef GPU_ENABLE

void MCFDES::solve_GPU(ull count = DEFAULT_COUNT) {
    std::vector<double> prefsum_prob;

    _make_prefsum_prob(prefsum_prob);

    /* Учёт начального и граничных условий */
    for (ull i = 0; i <= n; i++) {
        result[0][i] = psi(x_i(i));
    }

    for (ull j = 1; j <= k; j++) {
        result[j][0] = phiL(t_k(j));
        result[j][n] = phiR(t_k(j));
    }

    ull grid_size = (n - 1) * k;

    /* Предпосчет граничных условий */
    double* borders_cpu = new double[n + 1 + k + k];

    for (ull i = 0; i < n + k + 1; i++) {
        if (i < n + 1) {
            borders_cpu[i] = psi(x_i(i));
        } else {
            borders_cpu[i] = phiL(t_k(i - n));
            borders_cpu[i + k] = phiR(t_k(i - n));
        }
    }

    double* borders_gpu;
    CSC(cudaMalloc(&borders_gpu, (n + 1 + k + k) * sizeof(double)));
    CSC(cudaMemcpy(borders_gpu, borders_cpu, (n + 1 + k + k) * sizeof(double), cudaMemcpyHostToDevice));

    delete[] borders_cpu;

    /* Предпосчет функции источника */
    double* source_cpu = new double[grid_size];

    for (ull i = 0; i < grid_size; i++) {
        source_cpu[i] = f(x_i((i % (n - 1)) + 1), t_k((i / (n - 1)) + 1));
    }

    double* source_gpu;
    CSC(cudaMalloc(&source_gpu, grid_size * sizeof(double)));
    CSC(cudaMemcpy(source_gpu, source_cpu, grid_size * sizeof(double), cudaMemcpyHostToDevice));

    delete[] source_cpu;

    /* Моделирование */
    double* result_gpu;
    double* prefsum_prob_gpu;
    curandState* dev_states;
    CSC(cudaMalloc(&result_gpu, grid_size * sizeof(double)));
    CSC(cudaMalloc(&prefsum_prob_gpu, prefsum_prob.size() * sizeof(double)));
    CSC(cudaMalloc(&dev_states, grid_size * sizeof(curandState)));

    CSC(cudaMemcpy(prefsum_prob_gpu, prefsum_prob.data(), prefsum_prob.size() * sizeof(double), cudaMemcpyHostToDevice));
    CSC(cudaMemset(result_gpu, 0.0, grid_size * sizeof(double)));

    setup_kernel<<<TBLOCKS, THREADS>>>(dev_states, grid_size);
    CSC(cudaGetLastError());

    solve_kernel<<<TBLOCKS, THREADS>>>(dev_states, result_gpu, prefsum_prob_gpu, n, k, count, L, h, tau, gamma, borders_gpu, source_gpu);
    CSC(cudaGetLastError());

    double* result_cpu = new double[grid_size];
    CSC(cudaMemcpy(result_cpu, result_gpu, grid_size * sizeof(double), cudaMemcpyDeviceToHost));


    for (ull i = 0; i < grid_size; i++) {
        result[(i / (n - 1)) + 1][(i % (n - 1)) + 1] = result_cpu[i]; 
    }
    
    delete[] result_cpu;

    cudaFree(result_gpu);
    cudaFree(dev_states);
    cudaFree(prefsum_prob_gpu);
    cudaFree(borders_gpu);
    cudaFree(source_gpu);
}

#endif

std::vector<double> MCFDES::_make_prob() {
    std::vector<double> probabilities(2 * n + 2 + k , 0.0);

    double a00 = a(0.0, 0.0), b00 = b(0.0, 0.0), c00 = c(0.0, 0.0); // Только для D = const

    probabilities[n - 1] = a00 + b00 * g(alpha, 2.0) - c00;
    probabilities[n] = gamma - alpha * (a00 + b00);
    probabilities[n + 1] = a00 * g(alpha, 2.0) + b00 + c00;

    for (ull i = 2; i <= n; i++) {
        probabilities[n + i] = a00 * g(alpha, i + 1);
        probabilities[n - i] = b00 * g(alpha, i + 1);
    }

    for (ull i = 1; i <= k; i++) {
        probabilities[2 * n + i] = -g(gamma, i + 1);
    }

    probabilities[2 * n + 1 + k] = 0.0;
    probabilities[2 * n + 1 + k] = 1.0 - std::accumulate(probabilities.begin(), probabilities.end(), 0.0);

    return probabilities;
}

void MCFDES::_make_prefsum_prob(std::vector<double>& prefsum_prob) {
    prefsum_prob = prefix_sum(probabilities);
}

double MCFDES::drand(double dmin, double dmax) {
    return dmin + (double)rand() / RAND_MAX * (dmax - dmin);
}