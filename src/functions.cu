#include "functions.cuh"

/**************** R I G H T  H A N D  S I D E ****************************************************/
__global__ void rhs(REAL *v_out, REAL *v_in, size_t N, REAL *f_x, REAL *f_y, REAL mass, REAL betta) {

    size_t tid = get_tid();
    size_t stride = get_stride();

    for ( ; tid<N; tid+=stride) {
        // COORDINATES
        get_x(tid, v_out, N) = get_vx(tid, v_in, N);
        get_y(tid, v_out, N) = get_vy(tid, v_in, N);
        // VELOCITIES
        get_vx(tid, v_out, N) = f_x[tid] / mass - betta * get_vx(tid, v_in, N);
        get_vy(tid, v_out, N) = f_y[tid] / mass - betta * get_vy(tid, v_in, N);
    }
}

/*********** A R E A  F O R   4 N - V E C T O R  : X Y VX VY *************/
// will not work for other structures
__global__ void calculate_area(Lock lock, REAL *area, REAL *var, size_t N) {

    __shared__ REAL cache[THREADS]; //threads per block
    size_t tid = get_tid();
    size_t stride = get_stride();

    REAL temp = 0;
    for ( ; tid<N; tid+=stride) {
        temp += 0.5 * ( get_x(tid, var, N) * get_y(tid+1, var, N) - get_x(tid+1, var, N) * get_y(tid, var, N) );
    }

    cache[threadIdx.x] = temp;

    __syncthreads();

    int i = blockDim.x/2; // TODO general case, now it MUST BE PowOf2
    while (i != 0) {
        if (threadIdx.x < i)
            cache[threadIdx.x] += cache[threadIdx.x + i];
        __syncthreads();
        i /= 2;
    }

    if (threadIdx.x == 0) {
        lock.lock();
        *area += cache[0];
        //printf("calculate_area : treadIdx.x = %d || blockIdx.x = %d\n", threadIdx.x, blockIdx.x);
        //printf("area = %f\n",*area);
        lock.unlock();
    }
}