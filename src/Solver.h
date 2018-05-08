//
// Created by andrey on 02/05/18.
//

#ifndef RBC_SOLVER_H
#define RBC_SOLVER_H

#include "Vesicle.cuh"
#include "Walls.h"

class Solver {

public:
    virtual Vesicle do_step(Vesicle const & vesicle,
                    Walls const & walls,
                    double t_step,
                    void (*FuncName)) = 0;

    virtual ~Solver();

};


#endif //RBC_SOLVER_H
