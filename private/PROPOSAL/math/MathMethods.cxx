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

    std::vector<SplineCoefficients> CalculateSpline(std::vector<double> x, std::vector<double> y){
        //Algorithm from https://en.wikipedia.org/wiki/Spline_(mathematics)
        int n = x.size() - 1;
        if(x.size()!=y.size()){
            log_error("CalculateSpline: x and y (abscissa and ordinate) must have same dimension");
        }
        //1
        std::vector<double> a(n+1);
        for(int i = 0; i<n+1; i++){
            a[i] = y[i];
        }
        //2
        std::vector<double> b(n);
        std::vector<double> d(n);
        //3
        std::vector<double> h(n);
        for(int i = 0; i<n; i++){
            h[i] = x[i+1] - x[i];
        }
        //4
        std::vector<double> alpha(n);
        for(int i = 1; i<n; i++){
            alpha[i] = 3./h[i] * (a[i+1] - a[i]) - 3./h[i-1] * (a[i] - a[i-1]);
        }
        //5
        std::vector<double> c(n+1);
        std::vector<double> l(n+1);
        std::vector<double> mu(n+1);
        std::vector<double> z(n+1);
        //6
        l[0] = 1;
        mu[0] = 0;
        z[0] = 0;
        //7
        for(int i = 1; i<n; i++){
            l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i-1] * z[i-1])/l[i];
        }
        //8
        l[n] = 1;
        z[n] = 0;
        c[n] = 0;
        //9
        for(int j = n-1; j>=0; j--){
            c[j] = z[j] - mu[j] * c[j+1];
            b[j] = (a[j+1] - a[j])/h[j] - (h[j] * (c[j+1] + 2 * c[j]))/3;
            d[j] = (c[j+1] - c[j]) / (3 * h[j]);
        }
        //10ff
        std::vector<SplineCoefficients> return_values(n);
        for(int i = 0; i<n; i++){
            return_values[i] = SplineCoefficients(a[i], b[i], c[i], d[i], x[i]);
        }

        return return_values;

    }


}