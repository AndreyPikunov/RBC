//
// Created by andrey on 07/05/18.
//

#ifndef RBC_LOCK_CUH
#define RBC_LOCK_CUH

/****** M U T E X  S T R U C T U R E **********/
// from the book "CUDA by Example"
struct Lock {
    int *mutex;
    Lock( void ) {
        int state = 0;
        cudaMalloc((void **) &mutex,
                   sizeof(int));
        cudaMemcpy(mutex, &state, sizeof(int),
                   cudaMemcpyHostToDevice);
    }

    ~Lock( void ) {
        cudaFree( mutex );
    }

    __device__ void lock( void ) {
        while (atomicCAS(mutex, 0, 1) != 0);
    }

    __device__ void unlock( void ) {
        atomicExch(mutex, 0);
    }
};


#endif //RBC_LOCK_CUH
