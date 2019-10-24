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

Spline::Spline(const Spline& spline)
    : splines_(spline.splines_),
      subintervall_(spline.subintervall_),
      n_subintervalls_(spline.n_subintervalls_) {}

Spline::Spline(std::string path, bool binary) {
    Table_read reader(path, binary);
    reader.jump(2);

    reader.read(*this);
}

bool Spline::operator==(const Spline& spline) const {
    if(splines_ != spline.splines_)
        return false;
    if( n_subintervalls_ != spline.n_subintervalls_ )
        return false;
    if( subintervall_ != spline.subintervall_ )
        return false;
    return true;
}

bool Spline::operator!=(const Spline& spline) const {
    return !(*this == spline);
}

double Spline::evaluate(double x) {
    for (unsigned int i = 0; i < n_subintervalls_; ++i) {
        if (subintervall_[i] <= x && x <= subintervall_[i + 1])
            return splines_[i].evaluate(x);
    }

    // logger warining
    if (x > subintervall_.back())
        return splines_[n_subintervalls_ - 1].evaluate(x);
    else
        return splines_[0].evaluate(x);
}

void Spline::Derivative() {
    for (auto spline : splines_)
        spline = spline.GetDerivative();
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
                                     subintervall_.back());
    return domain;
}

std::vector<spline_container> Spline::GetSplineContainer() const {
    std::vector<spline_container> containers;

    spline_container container;
    for (unsigned int i = 0; i < splines_.size(); ++i) {
        container.coeff = splines_[i].GetCoefficient();
        container.domain.first = subintervall_[i];
        container.domain.second = subintervall_[i + 1];
        containers.push_back(container);
    }

    return containers;
}

bool Spline::save(std::string path, bool binary) {
    std::stringstream sstr;
    sstr << *this;

    Table_write writer(path, binary);
    writer.write(sstr.str());
    writer.close();

    return 1;
}

namespace PROPOSAL {

std::ostream& operator<<(std::ostream& os, const Spline& s) {
    os << "# n_splines" << std::endl;
    os << "# xmin xmax degree a_0 ... a_i" << std::endl;
    os << s.n_subintervalls_ << std::endl;
    for (spline_container i : s.GetSplineContainer())
        os << i << std::endl;
    return os;
}

std::istream& operator>>(std::istream& is, Spline& s) {
    spline_container container;

    is >> s.n_subintervalls_;
    for (unsigned int i = 0; i < s.n_subintervalls_; ++i) {
        is >> container;
        s.subintervall_.push_back(container.domain.first);
        s.splines_.push_back(Polynom(container.coeff));
        container.coeff.clear();
    }
    s.subintervall_.push_back(container.domain.second);
    return is;
}

std::ostream& operator<<(std::ostream& stream, spline_container& s) {
    stream << s.domain.first << " " << s.domain.second << " "
           << s.coeff.size() - 1 << " ";
    for (auto a_i : s.coeff)
        stream << a_i << " ";
    return stream;
}

std::istream& operator>>(std::istream& is, spline_container& s) {
    unsigned int degree;
    is >> s.domain.first >> s.domain.second >> degree;
    double a_i;
    for (unsigned int i = 0; i < degree + 1; ++i) {
        is >> a_i;
        s.coeff.push_back(a_i);
    }
    return is;
}

}  // namespace PROPOSAL
//----------------------------------------------------------------------------//
//------------------------------- Linear Spline
//------------------------------//
//----------------------------------------------------------------------------//

Linear_Spline::Linear_Spline(std::vector<double> x, std::vector<double> y)
    : Spline(x, y) {
    calculate_splines(x_, y_);
}

Linear_Spline::Linear_Spline(std::vector<Polynom> splines,
                             std::vector<double> subintervall)
    : Spline(splines, subintervall) {}

Linear_Spline::Linear_Spline(std::string spline_path, bool binary)
    : Spline(spline_path, binary) {}

Linear_Spline::Linear_Spline(const Spline& spline) : Spline(spline) {}

void Linear_Spline::calculate_splines(std::vector<double> x,
                                      std::vector<double> y) {
    unsigned int n = x.size() - 1;
    double a_1;
    double a_0;

    for (unsigned int i = 0; i < n; ++i) {
        a_1 = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        a_0 = y[i] - a_1 * x[i];
        Polynom p = Polynom({a_0, a_1});

        splines_.push_back(p);
    }

    n_subintervalls_ = n;
    subintervall_ = x;
}

//----------------------------------------------------------------------------//
//------------------------------- Cubic Spline
//-------------------------------//
//----------------------------------------------------------------------------//

Cubic_Spline::Cubic_Spline(std::vector<double> x, std::vector<double> y)
    : Spline(x, y) {
    calculate_splines(x_, y_);
}

Cubic_Spline::Cubic_Spline(std::vector<Polynom> splines,
                           std::vector<double> subintervall)
    : Spline(splines, subintervall) {}

Cubic_Spline::Cubic_Spline(std::string spline_path, bool binary)
    : Spline(spline_path, binary) {}

Cubic_Spline::Cubic_Spline(const Spline& spline) : Spline(spline) {}

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
