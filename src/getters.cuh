//
// Created by andrey on 06/05/18.
//

#ifndef RBC_GETTERS_CUH
#define RBC_GETTERS_CUH

#include <cmath>
#include "DEFINITIONS.cuh"


/******************************************/
///       X     ||     Y     ||     Vx    ||    Vy
/// |  |  |  |  ||  |  |  |  ||  |  |  |  ||  |  |  |  ||
/// 0            N           2N           3N           4N
// (N+i)%N for avoiding (negative % N)

inline __host__ __device__
REAL& get_x(int i, REAL *array, int N) {
    return array[(N+i)%N];
}

inline __host__ __device__
REAL& get_y(int i, REAL *array, int N) {
    return array[N+(N+i)%N];
}

inline __host__ __device__
REAL& get_vx(int i, REAL *array, int N) {
    return array[2*N+(N+i)%N];
}

inline __host__ __device__
REAL& get_vy(int i, REAL *array, int N) {
    return array[3*N+(N+i)%N];
}

#endif // RBC_GETTERS_CUH