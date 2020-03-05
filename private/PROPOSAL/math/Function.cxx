#include "PROPOSAL/math/Function.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include "PROPOSAL/Constants.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Polynom      %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace PROPOSAL;

Polynom::Polynom(std::vector<double> coefficients) : N_(coefficients.size()) {
    coeff_ = new double[N_];

    std::copy(coefficients.begin(), coefficients.end(), coeff_);
}

Polynom::Polynom(const Polynom& poly) : N_(poly.N_), coeff_(poly.coeff_) {}

bool Polynom::operator==(const Polynom& polynom) const {
    if (N_ != polynom.N_)
        return false;
    for (int i = 0; i < N_; ++i) {
        if(coeff_[i] != polynom.coeff_[i])
            return false;
    }
    return true;
}

bool Polynom::operator!=(const Polynom& polynom) const {
    return !(*this == polynom);
}

double Polynom::evaluate(double x) {
    double aux = coeff_[N_ - 1];

    for (int i = N_ - 2; i >= 0; --i)
        aux = aux * x + coeff_[i];

    return aux;
}

void Polynom::shift(double x) {
    // Shaw and Traub method for the Taylor shift
    // https://planetcalc.com/7726/#fnref1:shaw

    if (std::fabs(x) > GEOMETRY_PRECISION) {
        int n = N_ - 1;
        double** t = new double*[N_];
        for (int count = 0; count < N_; ++count)
            t[count] = new double[N_];

        for (int i = 0; i < n; ++i) {
            t[i][0] = coeff_[n - i - 1] * std::pow(x, n - i - 1);
            t[i][i + 1] = coeff_[n] * std::pow(x, n);
        }

        for (int j = 0; j <= n - 1; ++j) {
            for (int i = j + 1; i <= n; ++i) {
                t[i][j + 1] = t[i - 1][j] + t[i - 1][j + 1];
            }
        }

        for (int i = 0; i <= n - 1; ++i) {
            coeff_[i] = t[n][i + 1] / std::pow(x, i);
        }
    }
}

Polynom Polynom::GetDerivative() {
    std::vector<double> derivative_coeff_;

    for (auto i = 1; i < N_; ++i)
        derivative_coeff_.push_back(coeff_[i] * i);

    return Polynom(derivative_coeff_);
}

Polynom Polynom::GetAntiderivative(double constant) {
    std::vector<double> derivative_coeff_{constant};

    for (auto i = 0; i < N_; ++i)
        derivative_coeff_.push_back(coeff_[i] / (i + 1));

    return Polynom(derivative_coeff_);
}

std::vector<double> Polynom::GetCoefficient() const {
    std::vector<double> v(coeff_, coeff_ + N_);
    return v;
}

std::function<double(double)> Polynom::GetFunction() {
    return (std::function<double(double)>)std::bind(&Polynom::evaluate, this,
                                                    std::placeholders::_1);
}

namespace PROPOSAL {
std::ostream& operator<<(std::ostream& os, const Polynom& p) {
    os << "p(x) =";
    for (int i = 0; i < p.N_; ++i) {
        if (p.coeff_[i] != 0) {
            if (!std::signbit(p.coeff_[i]))
                os << "+";
            os << p.coeff_[i] << "*x^{" << i << "}";
        }
    }
    return os;
}
}  // namespace PROPOSAL
