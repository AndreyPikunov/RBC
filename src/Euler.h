//
// Created by andrey on 02/05/18.
//

#ifndef RBC_EULER_H
#define RBC_EULER_H

#include "Solver.h"


class Euler final : Solver {

    Vesicle do_step(Vesicle const &vesicle,
                    Walls const &walls,
                    double t_step,
                    void (*FuncName));

};


#endif //RBC_EULER_H
