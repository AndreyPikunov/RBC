//
// Created by andrey on 02/05/18.
//

#ifndef RBC_VESICLE_H
#define RBC_VESICLE_H

//#include "Walls.h"
//#include "kernels.cuh"
#include <cstddef>
#include <cassert>

template <typename Type>
struct Vesicle {

    __host__ __device__
    Vesicle(int N, Type *variables_) : N(N), variables(new Type[4*N]) {
        for (int i=0; i<4*N; ++i)
            this->variables[i] = variables_[i];
    };

    //Vesicle(Walls &w); //TODO

    __host__ __device__
    ~Vesicle() {
        delete [] variables;
    };

    __host__ __device__
    size_t get_size() {
        return sizeof(Vesicle) + N * sizeof(Type);
    }

    __host__ __device__
    Type& get_x(int i) {
        return variables[i%N];
    }

    __host__ __device__
    Type& get_y(int i) {
        return variables[N + i%N];
    }

    __host__ __device__
    Type& get_vx(int i) {
        return variables[2*N + i%N];
    }

    __host__ __device__
    Type& get_vy(int i) {
        return variables[3 * N + i % N];
    }

    __host__ __device__
    Vesicle<Type> operator+(Vesicle<Type> const &other) {
        assert(this->N==other.N);
        Vesicle<Type> v(this->N, this->variables);
        for (int i=0; i<4*v.N; ++i)
            v.variables[i] += other.variables[i];
        //auto *answer = new double[v1.N];
        //cudaMalloc( (void**)&answer , v1.N * sizeof(double) );

        return v;
    }

    template <typename Multiplier>
    __host__ __device__
    Vesicle<Type> operator*(Multiplier d) {
        Vesicle<Type> v_(this->N, this->variables);
        for (int i=0; i<4*v_.N; ++i)
            v_.variables[i] *= d;
        return v_;
    }

    int N;
    Type *variables;

};

//template <typename Type>                       <-        DOES NOT WORK!!!
//template <typename Multiplier>
//Vesicle<Type> Vesicle<Type>::operator*(Multiplier d, Vesicle<Type> const &v) {
    //return v * d;}



#endif //RBC_VESICLE_H
