#ifndef RBC_FORCES_CUH
#define RBC_FORCES_CUH

#include <cstddef> //size_t
#include <cstdio>
#include "kernels.cuh"
#include "getters.cuh"
#include "functions.cuh"
#include "DEFINITIONS.cuh"

__global__ void calculate_fl(REAL *var, size_t N, REAL *fl_x, REAL *fl_y, REAL l0 = 0.66, REAL kl = 0.3);

__global__ void calculate_fs(REAL *var, REAL *area, size_t N, REAL *fs_x, REAL *fs_y, REAL s0 = 314.0, REAL ks = 2000.0);

__global__ void calculate_fb(REAL *var, size_t N, REAL *fb_x, REAL *fb_y, REAL kb = 0.05);

//__global__ void fr(Vesicle *ves, double lr = 1.0, double kr = 1.0, Segment *segm, Vesicle *ves_out);

#endif // RBC_FORCES_CUH