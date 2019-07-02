#include <math.h>
#include "PROPOSAL/math/MathMethods.h"

namespace PROPOSAL {

    double NewtonRaphson(std::function<double(double)> func, std::function<double(double)> dfunc, double x1, double x2,
                         double xacc) {

        /*
         * Method adapted from rtsafe method from:
         * Numerical Recipes in C (2Nd Ed.): The Art of Scientific Computing, 1992
         * Press, William H. and Teukolsky, Saul A. and Vetterling, William T. and Flannery, Brian P.
         */

        double df, dx, dxold, f, fh, fl;
        double temp, xh, xl, rts;
        const int MAXIT = 100;

        fl = func(x1);
        fh = func(x2);

        if (fl * fh > 0.0) {
            log_error("Root must be bracketed in NetwonRaphson method!");
            return 0;
        }

        if (fl == 0.0) {
            return x1;
        }
        if (fh == 0.0) {
            return x2;
        }

        if (fl < 0.0) {
            //swap orientation of search
            xl = x1;
            xh = x2;
        } else {
            xh = x1;
            xl = x2;
        }

        rts = 0.5 * (x1 + x2); //initial guess for root
        dxold = fabs(x2 - x1);
        dx = dxold;

        f = func(rts);
        df = dfunc(rts);

        for (int j = 1; j < MAXIT; j++) {
            if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) || (fabs(2.0 * f) > fabs(dxold * df)))
            {
                //use bisection if Newton method produces an x that is out of range
                //or if Newton method is not decreasing fast enough
                dxold = dx;
                dx = 0.5 * (xh - xl);
                rts = xl + dx;
                if (xl == rts) {
                    return rts;
                }
            }
            else{
                //use Newton step
                dxold = dx;
                dx = f / df;
                temp = rts;
                rts -= dx;
                if (temp == rts) {
                    return rts;
                }
            }

            if (fabs(dx) < xacc) {
                //convergence criterion reached
                return rts;
            }
            f = func(rts);
            df = dfunc(rts);
            // update bracket
            if (f < 0.0) {
                xl = rts;
            } else {
                xh = rts;
            };
        }
        log_error("Maxmum number of iteration exeeded in NewtonRaphson");
        return 0.0;
    }

}