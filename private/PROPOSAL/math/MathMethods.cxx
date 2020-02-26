
#include <cmath>
#include <utility>
#include <fstream>
#include <iostream>
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/json.hpp"

namespace PROPOSAL {

double inverseErrorFunction(double p) {
    if (p <= 0 || p >= 1) {
        log_fatal(
            "The inverse Error function can just handle values between 0 and "
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

// ------------------------------------------------------------------------- //
double dilog(double x) {
    double C_arr[] = {
        0.42996693560813697,  0.40975987533077105,  -0.01858843665014592,
        0.00145751084062268,  -0.00014304184442340, 0.00001588415541880,
        -0.00000190784959387, 0.00000024195180854,  -0.00000003193341274,
        0.00000000434545063,  -0.00000000060578480, 0.00000000008612098,
        -0.00000000001244332, 0.00000000000182256,  -0.00000000000027007,
        0.00000000000004042,  -0.00000000000000610, 0.00000000000000093,
        -0.00000000000000014, 0.00000000000000002};

    double HF = 0.5;

    if (x == 1)
        return PI * PI / 6.0;
    else if (x == -1)
        return -PI * PI / 12.0;
    else {
        double _T = -x;
        double Y, S, A, B0, B1, B2, _H, _ALFA;

        if (_T <= -2.0) {
            Y = -1.0 / (1.0 + _T);
            S = 1.0;
            B1 = std::log(-_T);
            B2 = std::log(1.0 + 1.0 / _T);
            A = -PI * PI / 3.0 + HF * (B1 * B1 - B2 * B2);
        } else if (_T < -1.0) {
            Y = -1.0 - _T;
            S = -1.0;
            A = std::log(-_T);
            A = -PI * PI / 6.0 + A * (A + std::log(1 + 1 / _T));
        } else if (_T <= -0.5) {
            Y = -(1 + _T) / _T;
            S = 1.0;
            A = std::log(-_T);
            A = -PI * PI / 6.0 + A * (-HF * A + std::log(1.0 + _T));
        } else if (_T < 0.0) {
            Y = -_T / (1 + _T);
            S = -1.0;
            A = std::log(1 + _T);
            A = HF * A * A;
        } else if (_T <= 1.0) {
            Y = _T;
            S = 1.0;
            A = 0.0;
        } else {
            Y = 1.0 / _T;
            S = -1.0;
            A = std::log(_T);
            A = PI * PI / 6 + HF * A * A;
        }

        _H = Y + Y - 1;
        _ALFA = _H + _H;
        B1 = 0.0;
        B2 = 0.0;

        for (int i = 19; i >= 0.0; --i) {
            B0 = C_arr[i] + _ALFA * B1 - B2;
            B2 = B1;
            B1 = B0;
        }

        return -(S * (B0 - _H * B2) + A);
    }

    return x;
}

// double trilog(double x)
// {
//     if (x == 0.)
//         return 0.;
//     else if (x == 1.)
//         return ZETA3;
//     else if (x == -1.)
//         return -0.75*ZETA3;
//     else if (std::abs(x) < 1.)
//     {
//         int max_iter = 100;
//         double uncertainty_limit = 1e-6;
//         double power_series = 0.25*x;
//         double x_k = x;
//         double addend;
//         for (int k=2; k <= max_iter; k++)
//         {
//             x_k = x_k*x;
//             addend = x_k / (k*k*k * (k+1)*(k+1));
//             if (k > 5 && std::abs(addend/power_series) < uncertainty_limit)
//                 break;
//             if (k == max_iter)
//                 log_fatal("reach max iterations in trilog");
//             power_series += addend;
//         }
//         return -4 + (3 - 3./x)*std::log(1 - x) + (2 + 1./x)*dilog(x) + power_series;
//     }
//     else if (x < -1)
//     {
//         double ln_x = std::log(-x);
//         return trilog(1./x) - ln_x*ln_x*ln_x/6 - PI*PI/6*ln_x;
//     }
//     else // (x > 1)
//     {
//         double ln_x = std::log(x);
//         return -trilog(1./x) - trilog(1 - 1./x) + ZETA3 + ln_x*ln_x*ln_x/6 - 0.5*ln_x*ln_x*std::log(x - 1) + PI*PI/6*ln_x;
//     }
// }

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
    log_warn("Maximum number of iteration exeeded in NewtonRaphson");
    return rts;
}

std::vector<SplineCoefficients> CalculateSpline(std::vector<double> x,
                                                std::vector<double> y) {
    // Algorithm from https://en.wikipedia.org/wiki/Spline_(mathematics)
    int n = x.size() - 1;
    if (x.size() != y.size()) {
        log_error(
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
        std::string expanded_coords_path =
            Helper::ResolvePath(path_to_coordinates, true);
        std::ifstream infilestream(expanded_coords_path);
        infilestream >> json_coords;
    } catch (const nlohmann::json::parse_error& e) {
        log_fatal("Unable parse \"%s\" as json file",
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
                        log_fatal("Dens_spline: x coordinate is not a double");
                    }
                }
            } else {
                log_fatal(
                    "Dens_spline: at least 2 coordinates must be available to "
                    "create a spline");
            }
        } else {
            log_fatal("Dens_spline: x object is not an array");
        }
    } else {
        log_fatal("Dens_spline: no x coordinates are found");
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
                        log_fatal("Dens_spline: y coordinate is not a double");
                    }
                }
            } else {
                log_fatal(
                    "Dens_spline: at least 2 coordinates must be available to "
                    "create a spline");
            }
        } else {
            log_fatal("Dens_spline: y object is not an array");
        }
    } else {
        log_fatal("Dens_spline: no y coordinates are found");
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

}  // namespace PROPOSAL
