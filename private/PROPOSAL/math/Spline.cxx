#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/Function.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/Spline.h"
#include "PROPOSAL/math/TableWriter.h"

using namespace PROPOSAL;

Spline::Spline(std::vector<double> x, std::vector<double> y) : x_(x), y_(y) {
    if (x.size() != y.size())
        log_error(
            "CalculateSpline: x and y (abscissa and ordinate) must have same "
            "dimension");
}

Spline::Spline(std::vector<Polynom> splines, std::vector<double> subintervall)
    : splines_(splines),
      subintervall_(subintervall),
      n_subintervalls_(subintervall_.size()) {}

Spline::Spline(std::string spline_path) {}

Spline::Spline(const Spline& spline)
    : splines_(spline.splines_),
      subintervall_(spline.subintervall_),
      n_subintervalls_(spline.n_subintervalls_) {}

double Spline::evaluate(double x) {
    for (unsigned int i = 0; i < n_subintervalls_; ++i) {
        if (subintervall_[i] <= x && x <= subintervall_[i + 1])
            return splines_[i].evaluate(x);
    }

    if (x > subintervall_.back())
        return splines_[n_subintervalls_ - 1].evaluate(x);
    else
        return splines_[0].evaluate(x);
}

void Spline::Derivative() {
    for (unsigned int i = 0; i < n_subintervalls_; ++i) {
        splines_[i] = splines_[i].GetDerivative();
    }
}

void Spline::Antiderivative(double c) {
    double aux1 = c;
    double aux2 = 0;
    for (unsigned int i = 0; i < n_subintervalls_; ++i) {
        aux2 = splines_[i].GetAntiderivative(0).evaluate(subintervall_[i]);
        splines_[i] = splines_[i].GetAntiderivative(aux1 - aux2);
        aux1 += splines_[i].evaluate(subintervall_[i + 1]) -
                splines_[i].evaluate(subintervall_[i]);
    }
}

std::pair<double, double> Spline::GetDomain() {
    std::pair<double, double> domain(subintervall_.front(),
                                     subintervall_.back());  // value init

    return domain;
}
////----------------------------------------------------------------------------//
////------------------------------- Linear Spline
///------------------------------//
////----------------------------------------------------------------------------//

// Linear_Spline::Linear_Spline(std::vector<double> x, std::vector<double>
// y) :
//    Spline(x, y)
//{
//    calculate_splines(x_, y_);
//}

// Linear_Spline::Linear_Spline(std::vector<Polynom> splines,
// std::vector<double> subintervall) :
//    Spline(splines, subintervall)
//{

//}

// Linear_Spline::Linear_Spline(std::string path) :
//    Spline(path)
//{
//}

// Linear_Spline::Linear_Spline(const Linear_Spline& spline) :
//    Spline( spline )
//{
//}

// void Linear_Spline::calculate_splines(std::vector<double> x,
// std::vector<double> y)
//{
//}

// std::vector<linear_spline_container> Linear_Spline::GetSplineContainer()
//{
//    std::vector<linear_spline_container> container;

//    linear_spline_container lin_container;
//    for (unsigned int i = 0; i < splines_.size() ; ++i) {
//        std::vector<double> coeff = splines_[i].GetCoefficient();
//        lin_container.h_i = subintervall_[i];
//        lin_container.a_0 = coeff[0];
//        lin_container.a_1 = coeff[1];
//        container.push_back(lin_container);
//    }

//    return container;
//}

//----------------------------------------------------------------------------//
//-------------------------------- Cubic Spline
//------------------------------//
//----------------------------------------------------------------------------//

Cubic_Spline::Cubic_Spline(std::vector<double> x, std::vector<double> y)
    : Spline(x, y) {
    calculate_splines(x_, y_);
}

Cubic_Spline::Cubic_Spline(std::vector<Polynom> splines,
                           std::vector<double> subintervall)
    : Spline(splines, subintervall) {}

// Cubic_Spline::Cubic_Spline(std::string path) :
//     Spline(path)
// {
// }

Cubic_Spline::Cubic_Spline(const Cubic_Spline& spline) : Spline(spline) {}

void Cubic_Spline::calculate_splines(std::vector<double> x,
                                     std::vector<double> y) {
    // Algorithm from https://en.wikipedia.org/wiki/Spline_(mathematics)
    int n = x.size() - 1;
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
    std::vector<double> coeff;
    for (int i = 0; i < n; i++) {
        coeff = {a[i], b[i], c[i], d[i]};
        auto p = Polynom(coeff);
        p.shift(-x[i]);
        splines_.push_back(p);
        subintervall_.push_back(x[i]);
    }
    subintervall_.push_back(x.back());
    n_subintervalls_ = n;
}

std::vector<cubic_spline_container> Cubic_Spline::GetSplineContainer() {
    std::vector<cubic_spline_container> container;

    cubic_spline_container cub_container;
    for (unsigned int i = 0; i < splines_.size(); ++i) {
        std::vector<double> coeff = splines_[i].GetCoefficient();
        cub_container.h_i = subintervall_[i];
        cub_container.a_0 = coeff[0];
        cub_container.a_1 = coeff[1];
        cub_container.a_2 = coeff[2];
        cub_container.a_3 = coeff[3];
        container.push_back(cub_container);
    }

    return container;
}

//----------------------------------------------------------------------------//
//--------------------------------Save and
// Load-------------------------------//
//----------------------------------------------------------------------------//

// bool Linear_Spline::save(std::string path, bool binary)
// {

//         Table_write writer(path, binary);

//         writer.write(n_subintervalls_);

//         std::vector<linear_spline_container> spline_container =
//         GetSplineContainer(); for (auto i : spline_container) {
//             writer.write(i);
//         }

//         writer.close();

//         return 1;
// }

Cubic_Spline::Cubic_Spline(std::string path, bool binary) : Spline(path) {
    Table_read reader(path, binary);

    reader.read(n_subintervalls_);

    std::vector<cubic_spline_container> spline_container;
    cubic_spline_container container;

    for (unsigned int i = 0; i < n_subintervalls_; ++i) {
        reader.read(container);
        spline_container.push_back(container);
    }
    reader.close();

    for (auto i : spline_container) {
        subintervall_.push_back(i.h_i);

        std::vector<double> coeff = {i.a_0, i.a_1, i.a_2, i.a_3};
        Polynom p(coeff);
        splines_.push_back(p);
    }
    n_subintervalls_ = subintervall_.size();
}

bool Cubic_Spline::save(std::string path, bool binary) {
    Table_write writer(path, binary);

    writer.write(n_subintervalls_);

    std::vector<cubic_spline_container> spline_container = GetSplineContainer();
    for (auto i : spline_container) {
        writer.write(i);
    }
    writer.close();

    return 1;
}

