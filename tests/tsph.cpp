#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"

using namespace alglib;
void  function1_fvec(const real_1d_array &x, real_1d_array &fi, void *ptr)
{
    //
    // this callback calculates
    // f0(x0,x1) = 100*(x0+3)^4,
    // f1(x0,x1) = (x1-3)^4
    //
    fi[0] = x[0]*x[0] + x[1]*x[1] + (x[2] - 2) * (x[2] - 2) - 9;
    fi[1] = x[0]*x[0] + x[1]*x[1] + (x[2] + 2) * (x[2] + 2) - 9;
    //fi[0] = fi[0] * fi[0]; 
    //fi[0] = fi[1] * fi[1]; 
}
void  function1_jac(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr)
{
    //
    // this callback calculates
    // f0(x0,x1) = 100*(x0+3)^4,
    // f1(x0,x1) = (x1-3)^4
    // and Jacobian matrix J = [dfi/dxj]
    //
    fi[0] = x[0]*x[0] + x[1]*x[1] + (x[2] - 2) * (x[2] - 2) - 9;
    fi[1] = x[0]*x[0] + x[1]*x[1] + (x[2] + 2) * (x[2] + 2) - 9;
    //fi[0] = fi[0] * fi[0]; 
    //fi[0] = fi[1] * fi[1]; 
    jac[0][0] = 2 * x[0];
    jac[0][1] = 2 * x[1];
    jac[0][2] = 2 * x[2];
    jac[1][0] = 2 * x[0];
    jac[1][1] = 2 * x[1];
    jac[1][2] = 2 * x[2];
}

int main(int argc, char **argv)
{
    //
    // This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
    //
    //     f0(x0,x1) = 10*(x0+3)^2
    //     f1(x0,x1) = (x1-3)^2
    //
    // using "VJ" mode of the Levenberg-Marquardt optimizer.
    //
    // Optimization algorithm uses:
    // * function vector f[] = {f1,f2}
    // * Jacobian matrix J = {dfi/dxj}.
    //
    real_1d_array x = "[2,2,0]";
    real_1d_array bndl = "[1,1,-1]";
    real_1d_array bndu = "[3,3,1]";
    double epsx = 0.0000000001;
    ae_int_t maxits = 0;
    minlmstate state;
    minlmreport rep;


    minlmcreatev(2, x, 0.0001, state);
    //minlmcreatevj(2, x, state);

    minlmsetcond(state, epsx, maxits);
    minlmsetbc(state, bndl, bndu);
    alglib::minlmoptimize(state, function1_fvec);
//    alglib::minlmoptimize(state, function1_fvec, function1_jac);
    minlmresults(state, x, rep);


    double f0 = x[0]*x[0] + x[1]*x[1] + (x[2] - 2) * (x[2] - 2) - 9;
    double f1 = x[0]*x[0] + x[1]*x[1] + (x[2] + 2) * (x[2] + 2) - 9;


    printf("%s\n", x.tostring(8).c_str()); // EXPECTED: [-3,+3]
    printf("%lf, %lf\n", f0, f1);
    return 0;
}


