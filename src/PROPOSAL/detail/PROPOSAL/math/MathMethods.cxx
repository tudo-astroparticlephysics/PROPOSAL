
#include <cmath>
#include <utility>
#include <fstream>
#include <iostream>
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"
#include <nlohmann/json.hpp>
#include <limits>

namespace PROPOSAL {

double normalppf(double p) {
    if (p <= 0 || p >= 1) {
        Logging::Get("proposal.math")->critical(
            "The percent point function can just handle values between 0 and "
            "1.");
    }
    double a_arr[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                      -2.759285104469687e+02, 1.383577518672690e+02,
                      -3.066479806614716e+01, 2.506628277459239e+00};
    double b_arr[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                      -1.556989798598866e+02, 6.680131188771972e+01,
                      -1.328068155288572e+01};
    double c_arr[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                      -2.400758277161838e+00, -2.549732539343734e+00,
                      4.374664141464968e+00,  2.938163982698783e+00};
    double d_arr[] = {7.784695709041462e-03, 3.224671290700398e-01,
                      2.445134137142996e+00, 3.754408661907416e+00};
    double p_low = 0.02425;
    double p_high = 1 - p_low;

    double q, r, x, e, u;
    if (p < p_low) {
        // Rational approximation for lower region.
        q = std::sqrt(-2 * std::log(p));
        x = (((((c_arr[0] * q + c_arr[1]) * q + c_arr[2]) * q + c_arr[3]) * q +
              c_arr[4]) *
                 q +
             c_arr[5]) /
            ((((d_arr[0] * q + d_arr[1]) * q + d_arr[2]) * q + d_arr[3]) * q +
             1);
    } else if (p_low <= p && p <= p_high) {
        // Rational approximation for central region.
        q = p - 0.5;
        r = q * q;
        x = (((((a_arr[0] * r + a_arr[1]) * r + a_arr[2]) * r + a_arr[3]) * r +
              a_arr[4]) *
                 r +
             a_arr[5]) *
            q /
            (((((b_arr[0] * r + b_arr[1]) * r + b_arr[2]) * r + b_arr[3]) * r +
              b_arr[4]) *
                 r +
             1);
    } else {
        // Rational approximation for upper region.
        q = std::sqrt(-2 * std::log(1 - p));
        x = -(((((c_arr[0] * q + c_arr[1]) * q + c_arr[2]) * q + c_arr[3]) * q +
               c_arr[4]) *
                  q +
              c_arr[5]) /
            ((((d_arr[0] * q + d_arr[1]) * q + d_arr[2]) * q + d_arr[3]) * q +
             1);
    }

    // Refining the result:
    // One iteration of Halley’s rational method (third order)
    // gives full machine precision.
    e = 0.5 * std::erfc(-x / SQRT2) - p;
    u = e * std::sqrt(2 * PI) * std::exp(0.5 * x * x);
    x = x - u / (1 + x * u * 0.5);
    return x;
}

/// dilog implementation from Alexander Voigt (https://orcid.org/0000-0001-8963-6512)
/// From the repository https://github.com/Expander/polylogarithm (version 7.0.0)
/// Distributed under MIT license

double dilog(double x)
{
    const double PI = 3.1415926535897932;
    const double P[] = {
            0.9999999999999999502e+0,
            -2.6883926818565423430e+0,
            2.6477222699473109692e+0,
            -1.1538559607887416355e+0,
            2.0886077795020607837e-1,
            -1.0859777134152463084e-2
    };
    const double Q[] = {
            1.0000000000000000000e+0,
            -2.9383926818565635485e+0,
            3.2712093293018635389e+0,
            -1.7076702173954289421e+0,
            4.1596017228400603836e-1,
            -3.9801343754084482956e-2,
            8.2743668974466659035e-4
    };

    double y = 0, r = 0, s = 1;

    // transform to [0, 1/2]
    if (x < -1) {
        const double l = std::log(1 - x);
        y = 1/(1 - x);
        r = -PI*PI/6 + l*(0.5*l - std::log(-x));
        s = 1;
    } else if (x == -1) {
        return -PI*PI/12;
    } else if (x < 0) {
        const double l = std::log1p(-x);
        y = x/(x - 1);
        r = -0.5*l*l;
        s = -1;
    } else if (x == 0) {
        return x;
    } else if (x < 0.5) {
        y = x;
        r = 0;
        s = 1;
    } else if (x < 1) {
        y = 1 - x;
        r = PI*PI/6 - std::log(x)*std::log1p(-x);
        s = -1;
    } else if (x == 1) {
        return PI*PI/6;
    } else if (x < 2) {
        const double l = std::log(x);
        y = 1 - 1/x;
        r = PI*PI/6 - l*(std::log(y) + 0.5*l);
        s = 1;
    } else {
        const double l = std::log(x);
        y = 1/x;
        r = PI*PI/3 - 0.5*l*l;
        s = -1;
    }

    const double y2 = y*y;
    const double y4 = y2*y2;
    const double p = P[0] + y * P[1] + y2 * (P[2] + y * P[3]) +
                     y4 * (P[4] + y * P[5]);
    const double q = Q[0] + y * Q[1] + y2 * (Q[2] + y * Q[3]) +
                     y4 * (Q[4] + y * Q[5] + y2 * Q[6]);

    return r + s*y*p/q;
}

double NewtonRaphson(std::function<double(double)> func,
                     std::function<double(double)> dfunc,
                     double x1,
                     double x2,
                     double xinit,
                     int MAX_STEPS,
                     double xacc) {
    /*
     * Method adapted from rtsafe method from:
     * Numerical Recipes in C (2Nd Ed.): The Art of Scientific Computing, 1992
     * Press, William H. and Teukolsky, Saul A. and Vetterling, William T. and
     * Flannery, Brian P.
     */

    double df, dx, dxold, f, fh, fl;
    double temp, xh, xl, rts;

    fl = func(x1);
    fh = func(x2);

    if (fl * fh > 0.0) {
        throw MathException("Root must be bracketed in NewtonRaphson method!");
        // Logging::Get("proposal.math")->error("Root must be bracketed in NewtonRaphson method!");
        // return 0;
    }

    if (fl == 0.0) {
        return x1;
    }
    if (fh == 0.0) {
        return x2;
    }

    if (fl < 0.0) {
        // swap orientation of search
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
    }

    rts = xinit;  // initial guess for root
    dxold = std::abs(x2 - x1);
    dx = dxold;

    f = func(rts);
    df = dfunc(rts);

    for (int j = 1; j < MAX_STEPS; j++) {
        if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) ||
            (std::abs(2.0 * f) > std::abs(dxold * df))) {
            // use bisection if Newton method produces an x that is out of range
            // or if Newton method is not decreasing fast enough
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if (xl == rts) {
                return rts;
            }
        } else {
            // use Newton step
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if (temp == rts) {
                return rts;
            }
        }

        if (std::abs(dx) < xacc) {
            // convergence criterion reached
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
    Logging::Get("proposal.math")->warn("Maximum number of iteration exeeded in NewtonRaphson");
    return rts;
}

std::pair<double, double> Bisection(std::function<double(double)> f, double x1,
                                    double x2, double xacc, double MAX_ITER) {
    if (f(x1) * f(x2) > 0.)
        throw MathException("Root must be bracketed in Bisection method!");

    if (x1 > x2) {
        std::swap(x1, x2);
    }

    double center;

    for (int i = 0; i<=MAX_ITER; i++) {
        center = (x1 + x2) / 2;
        if (sgn(f(center)) == sgn(f(x1)))
            x1 = center;
        else
            x2 = center;
        if (std::abs(x2 - x1) < xacc)
            return std::make_pair(x1, x2);
    }
    Logging::Get("proposal.math")->warn("Maximum number of iteration exceeded in Bisection");
    return std::make_pair(x1, x2);
}

std::vector<SplineCoefficients> CalculateSpline(std::vector<double> x,
                                                std::vector<double> y) {
    // Algorithm from https://en.wikipedia.org/wiki/Spline_(mathematics)
    int n = x.size() - 1;
    if (x.size() != y.size()) {
        Logging::Get("proposal.math")->error(
            "CalculateSpline: x and y (abscissa and ordinate) must have same "
            "dimension");
    }
    // 1
    std::vector<double> a(n + 1);
    for (int i = 0; i < n + 1; i++) {
        a[i] = y[i];
    }
    // 2
    std::vector<double> b(n);
    std::vector<double> d(n);
    // 3
    std::vector<double> h(n);
    for (int i = 0; i < n; i++) {
        h[i] = x[i + 1] - x[i];
    }
    // 4
    std::vector<double> alpha(n);
    for (int i = 1; i < n; i++) {
        alpha[i] =
            3. / h[i] * (a[i + 1] - a[i]) - 3. / h[i - 1] * (a[i] - a[i - 1]);
    }
    // 5
    std::vector<double> c(n + 1);
    std::vector<double> l(n + 1);
    std::vector<double> mu(n + 1);
    std::vector<double> z(n + 1);
    // 6
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;
    // 7
    for (int i = 1; i < n; i++) {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    // 8
    l[n] = 1;
    z[n] = 0;
    c[n] = 0;
    // 9
    for (int j = n - 1; j >= 0; j--) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - (h[j] * (c[j + 1] + 2 * c[j])) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }
    // 10ff
    std::vector<SplineCoefficients> return_values(n);
    for (int i = 0; i < n; i++) {
        return_values[i] = SplineCoefficients(a[i], b[i], c[i], d[i], x[i]);
    }

    return return_values;
}

std::pair<std::vector<double>, std::vector<double>> ParseSplineCoordinates(
    const std::string& path_to_coordinates) {
    nlohmann::json json_coords;
    try {
        std::ifstream infilestream(path_to_coordinates);
        infilestream >> json_coords;
    } catch (const nlohmann::json::parse_error& e) {
        Logging::Get("proposal.math")->critical("Unable parse \"%s\" as json file",
                  path_to_coordinates.c_str());
    }

    std::vector<double> x_vec;
    std::vector<double> y_vec;

    if (json_coords.find("x") != json_coords.end()) {
        if (json_coords["x"].is_array()) {
            if (json_coords["x"].size() >= 2) {
                x_vec.resize(
                    json_coords["x"].get<std::vector<double>>().size());
                for (size_t i = 0;
                     i < json_coords["x"].get<std::vector<double>>().size();
                     i++) {
                    if (json_coords["x"][i].is_number()) {
                        x_vec[i] = json_coords["x"][i].get<double>();
                    } else {
                        Logging::Get("proposal.math")->critical("Dens_spline: x coordinate is not a double");
                    }
                }
            } else {
                Logging::Get("proposal.math")->critical(
                    "Dens_spline: at least 2 coordinates must be available to "
                    "create a spline");
            }
        } else {
            Logging::Get("proposal.math")->critical("Dens_spline: x object is not an array");
        }
    } else {
        Logging::Get("proposal.math")->critical("Dens_spline: no x coordinates are found");
    }

    if (json_coords.find("y") != json_coords.end()) {
        if (json_coords["y"].is_array()) {
            if (json_coords["y"].size() >= 2) {
                y_vec.resize(
                    json_coords["y"].get<std::vector<double>>().size());
                for (size_t i = 0;
                     i < json_coords["y"].get<std::vector<double>>().size();
                     i++) {
                    if (json_coords["y"][i].is_number()) {
                        y_vec[i] = json_coords["y"][i].get<double>();
                    } else {
                        Logging::Get("proposal.math")->critical("Dens_spline: y coordinate is not a double");
                    }
                }
            } else {
                Logging::Get("proposal.math")->critical(
                    "Dens_spline: at least 2 coordinates must be available to "
                    "create a spline");
            }
        } else {
            Logging::Get("proposal.math")->critical("Dens_spline: y object is not an array");
        }
    } else {
        Logging::Get("proposal.math")->critical("Dens_spline: no y coordinates are found");
    }

    std::pair<std::vector<double>, std::vector<double>> return_val;
    return_val.first = x_vec;
    return_val.second = y_vec;
    return return_val;
}

std::pair<double, double> welfords_online_algorithm(double& new_Value, unsigned int& iter, double& mean, double& cov) {
    double new_mean = mean + (new_Value - mean) / static_cast<double>(iter);
    double new_cov = cov + ((new_Value - mean) * (new_Value - new_mean) - cov) / static_cast<double>(iter);

    return std::make_pair (new_mean, new_cov);
}

double SampleFromGaussian(double mean, double sigma, double rnd, double min, double max) {
    if(sigma == 0){
        return mean;
    }
    auto xlo = 0.5 + std::erf((min - mean) / (SQRT2 * sigma)) / 2;
    auto xhi = 0.5 + std::erf((max - mean) / (SQRT2 * sigma)) / 2;
    auto rndtmp = xlo + (xhi - xlo) * rnd;
    return sigma * normalppf(rndtmp)  + mean;
}

double SampleFromExponential(double p, double lambda) {
    return - log1p(-p) / lambda;
}

}  // namespace PROPOSAL
