#include "kernels.cuh"

__global__ void add(REAL *answer, REAL *a, REAL *b, size_t size) {
    size_t tid = get_tid();
    size_t stride = get_stride();
    for ( ; tid<size; tid+=stride)
        answer[tid] = a[tid] + b[tid];
}

__global__ void mul(REAL *answer, REAL *a, REAL mult, size_t size) {
    size_t tid = get_tid();
    size_t stride = get_stride();
    for ( ; tid<size; tid+=stride)
        answer[tid] = mult * a[tid];
}

__global__ void fill(REAL *a, REAL value, size_t size) {
    size_t tid = get_tid();
    size_t stride = get_stride();
    for ( ; tid<size; tid+=stride) {
        a[tid] = value;
    }
}