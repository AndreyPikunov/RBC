#ifndef RBC_FUNCTIONS_CUH
#define RBC_FUNCTIONS_CUH

#include <cstddef>
#include <cmath>
#include <cstdio>
#include "kernels.cuh"
#include "getters.cuh"
#include "Lock.cuh"


__global__ void rhs(REAL *v_out, REAL *v_in, size_t N, REAL *f_x, REAL *f_y, REAL mass = 1.0, REAL betta = 1.0);

__global__ void calculate_area(Lock lock, REAL *area, REAL *var, size_t N);

inline __host__ __device__
REAL length(REAL x1, REAL y1, REAL x2, REAL y2) {
    return std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

inline __host__ __device__
int sign(REAL x) {
    return (x > 0) - (x < 0);
}

#endif // RBC_FUNCTIONS_CUH