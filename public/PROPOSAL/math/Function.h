#pragma once
#include <functional>
#include <vector>

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Polynom      %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace PROPOSAL {

class Polynom {
   public:
    Polynom(std::vector<double> coefficients);
    Polynom(const Polynom&);

    // Polynom& operator=(const Polynom& poly);

    Polynom* clone() const { return new Polynom(*this); };

    double evaluate(double x);
    void shift(double x);

    Polynom GetDerivative();
    Polynom GetAntiderivative(double constant);
    std::vector<double> GetCoefficient();

    friend std::ostream& operator<<(std::ostream& os, const Polynom& p);

    std::function<double(double)> GetFunction();

   protected:
    int N;
    double* coeff;
};

}  // namespace PROPOSAL
