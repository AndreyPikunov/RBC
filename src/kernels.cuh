#ifndef RBC_KERNELS_CUH
#define RBC_KERNELS_CUH

#include <cstddef>
#include <cstdio>

#include "DEFINITIONS.cuh"

__global__ void add(REAL *answer, REAL *a, REAL *b, size_t size);

__global__ void mul(REAL *answer, REAL *a, REAL mult, size_t size);

__global__ void fill(REAL *a, REAL value, size_t size);


inline __device__ size_t get_tid() {
    return blockIdx.x * blockDim.x + threadIdx.x;
}

inline __device__ size_t get_stride() {
    return blockDim.x * gridDim.x;
}

#endif // RBC_KERNELS_CUH