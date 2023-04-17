# CC = g++
# CFLAGS = --std=c++14 -O3
CC = nvcc
CFLAGS = -Wno-deprecated-gpu-targets --extended-lambda -gencode arch=compute_50,code=sm_50 --std=c++14 -O3
SOURSES = main.cu lib/*.cpp lib/*.cu
BIN = main.out
TASK = task.out

all:
	$(CC) $(CFLAGS) -o $(BIN) $(SOURSES)
run:
	./main.out
taskk:
	$(CC) $(CFLAGS) -o $(TASK) task.cu lib/*.cpp lib/*.cu