//
// Created by andrey on 02/05/18.
//

#ifndef RBC_CONSTANTS_H
#define RBC_CONSTANTS_H

#include "DEFINITIONS.cuh"

struct Constants {

    explicit Constants(REAL mass = 1.0, REAL betta = 1.0,
                       REAL l0 = 0.66, REAL kl = 0.3,
                       REAL s0 = 314.0, REAL ks = 200.0,
                       REAL kb = 0.05,
                       REAL lr = 1.0, REAL kr = 1.0,
                       REAL lp = 1.0, REAL kp = 1.0) :
            mass(mass), betta(betta),
            l0(l0), kl(kl),
            s0(s0), ks(ks),
            kb(kb), lr(lr),
            kr(kr),
            lp(lp), kp(kp) { }

    REAL mass;
    REAL betta;

    REAL l0;
    REAL kl;

    REAL s0;
    REAL ks;

    REAL kb;

    REAL lr;
    REAL kr;

    REAL lp;
    REAL kp;
};


#endif //RBC_CONSTANTS_H
