#include <math.h>
#include <iostream>
#include "PROPOSAL/math/MathMethods.h"

#include "PROPOSAL/json.hpp"
#include <fstream>
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Logging.h"

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
            throw MathException("Root must be bracketed in NewtonRaphson method!");
            // log_error("Root must be bracketed in NewtonRaphson method!");
            // return 0;
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
        log_error("Maximum number of iteration exeeded in NewtonRaphson");
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

    std::pair<std::vector<double>, std::vector<double>> ParseSplineCoordinates(const std::string& path_to_coordinates){
        nlohmann::json json_coords;
        try
        {
            std::string expanded_coords_path = Helper::ResolvePath(path_to_coordinates, true);
            std::ifstream infilestream(expanded_coords_path);
            infilestream >> json_coords;
        }
        catch (const nlohmann::json::parse_error& e)
        {
            log_fatal("Unable parse \"%s\" as json file", path_to_coordinates.c_str());
        }

        std::vector<double> x_vec;
        std::vector<double> y_vec;

        if (json_coords.find("x") != json_coords.end())
        {
            if (json_coords["x"].is_array())
            {
                if (json_coords["x"].size() >= 2)
                {
                    x_vec.resize(json_coords["x"].get<std::vector<double>>().size());
                    for (size_t i=0; i < json_coords["x"].get<std::vector<double>>().size(); i++)
                    {
                        if (json_coords["x"][i].is_number())
                        {
                            x_vec[i] = json_coords["x"][i].get<double>();
                        }
                        else
                        {
                            log_fatal("Dens_spline: x coordinate is not a double");
                        }
                    }
                }
                else
                {
                    log_fatal("Dens_spline: at least 2 coordinates must be available to create a spline");
                }
            }
            else
            {
                log_fatal("Dens_spline: x object is not an array");
            }
        }
        else
        {
            log_fatal("Dens_spline: no x coordinates are found");
        }

        if (json_coords.find("y") != json_coords.end())
        {
            if (json_coords["y"].is_array())
            {
                if (json_coords["y"].size() >= 2)
                {
                    y_vec.resize(json_coords["y"].get<std::vector<double>>().size());
                    for (size_t i=0; i < json_coords["y"].get<std::vector<double>>().size(); i++)
                    {
                        if (json_coords["y"][i].is_number())
                        {
                            y_vec[i] = json_coords["y"][i].get<double>();
                        }
                        else
                        {
                            log_fatal("Dens_spline: y coordinate is not a double");
                        }
                    }
                }
                else
                {
                    log_fatal("Dens_spline: at least 2 coordinates must be available to create a spline");
                }
            }
            else
            {
                log_fatal("Dens_spline: y object is not an array");
            }
        }
        else
        {
            log_fatal("Dens_spline: no y coordinates are found");
        }

        std::pair<std::vector<double>, std::vector<double>> return_val;
        return_val.first  = x_vec;
        return_val.second = y_vec;
        return return_val;

    }

}
