#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <random>
#include <chrono>

#include "cuda_runtime.h"
#include "curand_kernel.h"
#include "device_launch_parameters.h"
#include <iostream>

#define CSC(result) do { \
	if (result != cudaSuccess) { \
		fprintf(stderr, "ERROR %d %s\n", __LINE__, cudaGetErrorString(result)); \
		exit(0); \
	} \
} while(0)

struct Tensor {

	double*** cage;
	size_t n1, n2, n3;

	Tensor(size_t _n1, size_t _n2, size_t _n3) : n1(_n1), n2(_n2), n3(_n3) {
		cage = new double**[n1 + 1];

		for (size_t k = 0; k <= n1; k++) {
			cage[k] = new double* [n2 + 1];

			for (size_t j = 0; j <= n2; j++) {
				cage[k][j] = new double[n3 + 1]{0.0};
			}
		}
	}

	void print() {
		for (size_t k = 0; k <= n1; k++) {
			for (size_t j = 0; j <= n2; j++) {
				for (size_t i = 0; i <= n3; i++) {
					printf("%lf ", cage[k][j][i]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}

	void fprint(std::string filename, double Lx, double hx, double Ly, double hy, size_t N) {
		FILE* out = fopen(filename.c_str(), "w");

		for (size_t k = 0; k <= n1; k++) {
			for (size_t i = 0; i <= n3; i++) {
				for (size_t j = 0; j <= n2; j++) {
					fprintf(out, "%d %d %lf\n", i, j, cage[k][j][i] / (double)N);
					//fprintf(out, "%lf ", cage[k][j][i]);
				}
			}
			fprintf(out, "\n");
		}

		//for (size_t k = 0; k <= n1; k++) {
		//	for (size_t i = 0; i <= n3; i++) {
		//		for (size_t j = 0; j <= n2; j++) {
		//			fprintf(out, "%lf %lf %lf\n", Lx + hx * (double)i, Ly + hy * (double)j, cage[k][j][i] / (double)N);
		//		}
		//		//fprintf(out, "\n");
		//	}
		//	fprintf(out, "\n");
		//}

		fclose(out);
	}

	~Tensor() {
		delete[] cage;
	}
};

__host__ double drand() { // 0.0 ... 1.0
	return ((double)rand()) / (double)RAND_MAX;
}

void MonteCarlo_CPU(Tensor& cage, double a, double b, double c, double d, double hx, size_t nx, double hy, size_t ny, double tau, size_t nt, size_t N) {
	
	double probabilities[5] = {
		1.0 - (2.0 * tau * a / hx / hx) - (2.0 * tau * b / hy / hy),
		(tau * a / hx / hx) + (tau * c / hx / 2.0),
		(tau * a / hx / hx) - (tau * c / hx / 2.0),
		(tau * b / hx / hx) + (tau * d / hx / 2.0),
		(tau * b / hx / hx) - (tau * d / hx / 2.0),
	};

	double sum_probabilities[5] = { probabilities[0] };

	for (size_t i = 1; i < 5; i++) {
		sum_probabilities[i] = sum_probabilities[i - 1] + probabilities[i];
	}

	/*for (size_t p = 0; p < 5; p++) {
		printf("%lf ", sum_probabilities[p]);
	}
	printf("\n");*/

	for (size_t particle = 0; particle < N; particle++) {
		size_t cur_x = nx / 2, cur_y = ny / 2;
		cage.cage[0][cur_y][cur_x] += 1.0;

		for (size_t cur_k = 1; cur_k <= nt; cur_k++) {
			double rnd = drand();

			int i = 0;
			while (i < 5 && rnd > sum_probabilities[i]) {
				i++;
			}

			switch (i) {
			case 0:
				break;
			case 1:
				cur_x--;
				break;
			case 2:
				cur_x++;
				break;
			case 3:
				cur_y--;
				break;
			case 4:
				cur_y++;
				break;
			}

			if (cur_x < 0 || cur_x > nx || cur_y < 0 || cur_y > ny) {
				break;
			}

			cage.cage[cur_k][cur_y][cur_x] += 1.0;
		}
	}
}

__global__ void setup_kernel(curandState* state, size_t count_of_threads) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int offset = gridDim.x * blockDim.x;

	while (idx < count_of_threads) {
		curand_init(1234, idx, 0, &state[idx]);
		idx += offset;
	}
}

__global__ void MonteCarlo_GPU(
	curandState* state,
	unsigned int *cage,
	double a, double b, double c, double d, 
	double hx, size_t nx, double hy, size_t ny, double tau, size_t nt, size_t N)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int offset = gridDim.x * blockDim.x;

	double probabilities[5] = {
		1.0 - (2.0 * tau * a / hx / hx) - (2.0 * tau * b / hy / hy),
		(tau * a / hx / hx) + (tau * c / hx / 2.0),
		(tau * a / hx / hx) - (tau * c / hx / 2.0),
		(tau * b / hx / hx) + (tau * d / hx / 2.0),
		(tau * b / hx / hx) - (tau * d / hx / 2.0),
	};

	double sum_probabilities[5] = { probabilities[0] };

	for (size_t i = 1; i < 5; i++) {
		sum_probabilities[i] = sum_probabilities[i - 1] + probabilities[i];
	}

	while (idx < N) {
		curandState localState = state[idx];

		size_t cur_x = nx / 2, cur_y = ny / 2;
		
		atomicAdd(cage + (cur_y * (nx + 1) + cur_x), 1);

		for (size_t cur_k = 1; cur_k <= nt; cur_k++) {
			double rnd = curand_uniform_double(&localState);

			int i = 0;
			while (i < 5 && rnd > sum_probabilities[i]) {
				i++;
			}

			switch (i) {
			case 0:
				break;
			case 1:
				cur_x--;
				break;
			case 2:
				cur_x++;
				break;
			case 3:
				cur_y--;
				break;
			case 4:
				cur_y++;
				break;
			}

			if (cur_x < 0 || cur_x > nx || cur_y < 0 || cur_y > ny) {
				break;
			}

			//cage[cur_k][cur_y][cur_x] += 1.0;
			atomicAdd(cage + (cur_k * (nx + 1) * (ny + 1) + cur_y * (nx + 1) + cur_x), 1);
			//atomicAdd(cage[cur_k][cur_y] + cur_x, 1.0);
		}

		idx += offset;
	}
}

int main()
{
	srand(1234);

	double a = 0.718, b = 0.03, c = 0.0, d = 0.0;
	double Lx = -1.0, Rx = 1.0, Ly = -1.0, Ry = 1.0, T = 0.5;
	size_t nx = 40, ny = 40, nt = 300;

	double hx = (Rx - Lx) / (double)nx;
	double hy = (Ry - Ly) / (double)ny;
	double tau = T / (double)nt;
	double sigma = (tau * a / hx / hx) + (tau * b / hy / hy);
	size_t N = 500;

	printf("sigma = %lf\n", sigma);
	if (sigma > 0.5) {
		printf("Condition not met (%lf > 0.5)\n", sigma);
		return 0;
	}

	for (size_t N = 200000; N < 200001; N++) {
		Tensor cage(nt, ny, nx);

		//cage.print();

		{
			unsigned int* cage_gpu;
			curandState* devStates;
			CSC(cudaMalloc(&cage_gpu, (nx + 1) * (ny + 1) * (nt + 1) * sizeof(unsigned int)));
			CSC(cudaMalloc(&devStates, N * sizeof(curandState)));

			CSC(cudaMemset(cage_gpu, 0, (nx + 1) * (ny + 1) * (nt + 1) * sizeof(unsigned int)));

			setup_kernel << <512, 512 >> > (devStates, N);
			CSC(cudaGetLastError());

			cudaEvent_t start, stop;
			float time;
			cudaEventCreate(&start);
			cudaEventCreate(&stop);
			cudaEventRecord(start, 0);

			MonteCarlo_GPU << <512, 512 >> > (devStates, cage_gpu, a, b, c, d, hx, nx, hy, ny, tau, nt, N);
			CSC(cudaGetLastError());

			cudaEventRecord(stop, 0);
			cudaEventSynchronize(stop);
			cudaEventElapsedTime(&time, start, stop);

			unsigned int* cage_cpu = (unsigned int*)malloc((nx + 1) * (ny + 1) * (nt + 1) * sizeof(unsigned int));
			CSC(cudaMemcpy(cage_cpu, cage_gpu, (nx + 1) * (ny + 1) * (nt + 1) * sizeof(unsigned int), cudaMemcpyDeviceToHost));

			for (size_t k = 0; k <= nt; k++) {
				for (size_t i = 0; i <= nx; i++) {
					for (size_t j = 0; j <= ny; j++) {
						cage.cage[k][j][i] = (double)cage_cpu[k * (ny + 1) * (nx + 1) + j * (nx + 1) + i];
					}
				}
			}

			printf("%llu %f\n", N, time);

			cudaFree(cage_gpu);
			cudaFree(devStates);
			free(cage_cpu);
		}

		cage.fprint("output.txt", Lx, hx, Ly, hy, N);
	}

	//for (size_t N = 10; N < 500000; N *= 2) {
	//	Tensor cage(nt, ny, nx);

	//	//cage.print();
	//	{
	//		auto start1 = std::chrono::system_clock::now();
	//		MonteCarlo_CPU(cage, a, b, c, d, hx, nx, hy, ny, tau, nt, N);
	//		auto end1 = std::chrono::system_clock::now();
	//		auto time = std::chrono::duration<float, std::ratio<1, 1000>>(end1 - start1).count();

	//		printf("%llu %f\n", N, time);
	//	}
	//	//cage.fprint("output.txt", Lx, hx, Ly, hy, N);
	//}
	

    return 0;
}