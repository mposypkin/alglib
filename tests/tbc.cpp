/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   tbc.cpp
 * Author: posypkin
 *
 * Created on November 19, 2018, 11:57 AM
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"

using namespace alglib;
void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr) 
{
    //
    // this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
    // and its derivatives df/dx0 and df/dx1
    //
    func = 100*pow(x[0]+3,4) + pow(x[1]-3,4);
    grad[0] = 400*pow(x[0]+3,3);
    grad[1] = 4*pow(x[1]-3,3);
}

int main(int argc, char **argv)
{
    //
    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
    // subject to bound constraints -1<=x<=+1, -1<=y<=+1, using BLEIC optimizer.
    //
    real_1d_array x = "[0,0]";
    real_1d_array bndl = "[-1,-1]";
    real_1d_array bndu = "[+1,+1]";
    minbleicstate state;
    minbleicreport rep;

    //
    // These variables define stopping conditions for the optimizer.
    //
    // We use very simple condition - |g|<=epsg
    //
    double epsg = 0.000001;
    double epsf = 0;
    double epsx = 0;
    ae_int_t maxits = 0;

    //
    // Now we are ready to actually optimize something:
    // * first we create optimizer
    // * we add boundary constraints
    // * we tune stopping conditions
    // * and, finally, optimize and obtain results...
    //
    minbleiccreate(x, state);
    minbleicsetbc(state, bndl, bndu);
    minbleicsetcond(state, epsg, epsf, epsx, maxits);
    alglib::minbleicoptimize(state, function1_grad);
    minbleicresults(state, x, rep);

    //
    // ...and evaluate these results
    //
    printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
    printf("%s\n", x.tostring(2).c_str()); // EXPECTED: [-1,1]
    return 0;
}


