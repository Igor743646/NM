# CC = g++
# CFLAGS = --std=c++14 -O3
CC = nvcc
CFLAGS = -Wno-deprecated-gpu-targets --extended-lambda -gencode arch=compute_50,code=sm_50 --std=c++14 -O3
SOURSES = *.cu lib/*.cpp lib/*.cu
BIN = main.out

all:
	$(CC) $(CFLAGS) -o $(BIN) $(SOURSES)

run:
	./main.out