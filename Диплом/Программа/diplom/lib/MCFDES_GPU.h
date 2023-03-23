#pragma once
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <curand_kernel.h>
#include "utils.h"

#define CSC(result) do {                                                         \
    if (result != cudaSuccess) {                                                 \
        fprintf(stderr, "ERROR %d %s\n", __LINE__, cudaGetErrorString(result));  \
        exit(0);                                                                 \
    }                                                                            \
} while(0)

__global__ void setup_kernel(curandState* state, ull count_of_threads) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = gridDim.x * blockDim.x;

    while (idx < count_of_threads) {
        curand_init(1234, idx, 0, &state[idx]);
        idx += offset;
    }
}

__device__ ull bin_search_gpu(const double* v, double k, size_t sz) {
    ull l = 0, r = sz - 1, mid;

    while (r - l > 1) {

        mid = (r + l) / 2;

        if (k < v[mid]) {
            r = mid;
        }
        else if (k > v[mid]) {
            l = mid;
        }
        else {
            l = mid;
            r = mid;
        }

    }

    return k <= v[l] ? l : r;
}


__global__ void solve_kernel(curandState* state, double* result, double* prefsum_prob_gpu,
    ull n, ull k, ull count,
    double L, double h, double tau, double gamma, double* borders_gpu, double* source_gpu,
    double alpha_l, double beta_l, double alpha_r, double beta_r
) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = gridDim.x * blockDim.x;

    curandState localState = state[idx];

    double p1_l = alpha_l / (alpha_l - beta_l * h);
    double p1_r = alpha_r / (alpha_r + beta_r * h);

    while (idx < (n + 1) * k) {

        for (ull _i = 0; _i < count; _i++) {
            long long y = (idx / (n + 1)) + 1;
            long long x = idx % (n + 1);

            if (x == 0) {
                double rnd1 = curand_uniform_double(&localState);
                if (rnd1 <= p1_l) { // если частица отразилась
                    // result[idx] += source_gpu[(y - 1) * (n + 1) + x] ;
                    x++; 
                } else { // если частица поглотилась
                    result[idx] += borders_gpu[n + 1 + y] / beta_l; 
                    continue;
                }
            } else if (x == n) {
                double rnd1 = curand_uniform_double(&localState);
                if (rnd1 <= p1_r) { // если частица отразилась
                    // result[idx] += source_gpu[(y - 1) * (n + 1) + x] / (beta_l + alpha_l / h) * h;
                    x--;
                } else { // если частица поглотилась
                    result[idx] += borders_gpu[n + 1 + k + y] / beta_r; 
                    continue;
                }
            }

            while (y > 0) {
                double rnd = curand_uniform_double(&localState);
                long long idy = 0;

                idy = (long long)bin_search_gpu(prefsum_prob_gpu, rnd, 2*n+k+2);

                result[idx] += source_gpu[(y - 1) * (n + 1) + x] * pow(tau, gamma);

                if (idy <= 2 * n) { // перемещение по пространству
                    x += idy - n;
                    y--;
                } else if (idy <= 2 * n + k) { // перемещение по времени
                    y -= idy - 2 * n + 1;
                } else {
                    break;
                }

                if (y == 0 && (x >= 0) && (x <= n)) {
                    result[idx] += borders_gpu[x];
                    break;
                } else if (y == 0) {
                    break;
                }

                if (x == 0) {
                    double rnd1 = curand_uniform_double(&localState);
                    if (rnd1 <= p1_l) { // если частица отразилась
                        // result[idx] += source_gpu[(y - 1) * (n + 1) + x] * pow(tau, gamma);
                        x++; 
                    } else { // если частица поглотилась
                        result[idx] += borders_gpu[n + 1 + y] / beta_l; 
                        break;
                    }
                } else if (x == n) {
                    double rnd1 = curand_uniform_double(&localState);
                    if (rnd1 <= p1_r) { // если частица отразилась
                        // result[idx] += source_gpu[(y - 1) * (n + 1) + x] / (beta_l + alpha_l / h) * h;
                        x--;
                    } else { // если частица поглотилась
                        result[idx] += borders_gpu[n + 1 + k + y] / beta_r; 
                        break;
                    }
                }
            }
            
            // if (y == 0 && (x >= 0) && (x <= n)) { // если частица дошла до начального момента времени
            //     result[idx] += borders_gpu[x];
            // } else if (x == 0 && y > 0) { // если частица попала на левую границу
            //     double rnd = curand_uniform_double(&localState);

            //     if (rnd < p1_l) { // если частица отразилась

            //     } else {

            //     }
            //     result[idx] += borders_gpu[n + 1 + y];
            // } else if (x == n && y > 0) {
            //     result[idx] += borders_gpu[n + 1 + k + y];
            // }
        }
        
        result[idx] /= (double)count;
        
        idx += offset;
    }
}