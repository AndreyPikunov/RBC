//
// Created by andrey on 02/05/18.
//

#include <algorithm>

#include "AdaptiveSolver.h"

AdaptiveSolver::~AdaptiveSolver() = default;

double AdaptiveSolver::adapt_t_step(Vesicle const &vesicle,
                                    Walls const &walls,
                                    double t_step,
                                    Params params,
                                    void *FuncName) {

    /*
    // Try one step with t_step
    Vesicle v1 = do_step(vesicle, walls, t_step, FuncName);
    // Try two steps with t_step/2
    Vesicle v2 = do_step(do_step(vesicle, walls, 0.5*t_step, FuncName),
                         walls, 0.5*t_step, FuncName);

    double delta = calculate_delta(v1,v2);

    double t_adapted = 0;

    if (delta != 0)
        t_adapted = t_step * params.marg * std::min( std::max(params.tol/delta, params.decr) , params.incr );
    else
        t_adapted = t_step * params.marg * params.incr;

    return t_adapted;
     */

}