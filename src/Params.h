//
// Created by andrey on 02/05/18.
//

#ifndef RBC_PARAMS_H
#define RBC_PARAMS_H

#include "DEFINITIONS.cuh"

struct Params {

    explicit Params(REAL t = 0.01, REAL m = 0.9, REAL d = 0.3, REAL i = 2.0) :
            tol(t), marg(m), decr(d), incr(i) {};

    REAL tol; // tolerance
    REAL marg; // margin
    REAL decr; // decrement
    REAL incr; // increment

};


#endif //RBC_PARAMS_H
