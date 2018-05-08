//
// Created by andrey on 02/05/18.
//

#ifndef RBC_ADAPTIVESOLVER_H
#define RBC_ADAPTIVESOLVER_H

#include "Solver.h"
#include "Params.h"

class AdaptiveSolver : Solver {

public:
    double adapt_t_step(Vesicle const & vesicle,
                        Walls const & walls,
                        double t_step,
                        Params params,
                        void (*FuncName));

    virtual ~AdaptiveSolver();

};


#endif //RBC_ADAPTIVESOLVER_H
