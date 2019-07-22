#include <functional>
#include <algorithm>
#include <iostream>
#include "PROPOSAL/math/Function.h"


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Polynom      %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace PROPOSAL;

Polynom::Polynom(std::vector<double> coefficients) :
    N( coefficients.size() )
{
    coeff = new double[N];

    std::copy(coefficients.begin(), coefficients.end(), coeff);
}

Polynom::Polynom(const Polynom& poly):
    N(poly.N),
    coeff(poly.coeff)
{}

double Polynom::evaluate(double x)
{
    double aux {0};

    // Hier sollte demnächst das Hornerschema stehen!
    for (int i = 0; i < N; ++i) 
        aux += coeff[i] * std::pow(x, i);

    return aux;
}

void Polynom::shift(double x)
{
    // Shaw and Traub method for the Taylor shift
    // https://planetcalc.com/7726/#fnref1:shaw
	
    int n = N-1;
    double **t = new double*[N]; 
    for (int count = 0; count < N; ++count)
        t[count] = new double[N]; 

    for (int i = 0; i < n; ++i) 
    {
        t[i][0] = coeff[n-i-1] * std::pow(x, n-i-1);
        t[i][i+1] = coeff[n] * std::pow(x, n);
    }

    for (int j = 0; j <= n-1; ++j) 
    {
        for (int i = j+1; i <= n; ++i) 
        {
            t[i][j+1] = t[i-1][j] + t[i-1][j+1];   
        }
    }

    for (int i = 0; i <= n-1; ++i) 
    {
        coeff[i] = t[n][i+1] / std::pow(x, i);
    }
}

Polynom Polynom::GetDerivative()
{
    std::vector<double> derivative_coeff;

    for (auto i = 1; i < N; ++i)
        derivative_coeff.push_back(coeff[i] * i);

    return Polynom(derivative_coeff);
}

Polynom Polynom::GetAntiderivative(double constant)
{
    std::vector<double> derivative_coeff {constant};

    for (auto i = 0; i < N; ++i)
        derivative_coeff.push_back(coeff[i] / (i+1));

    return Polynom(derivative_coeff);
}
        
std::vector<double> Polynom::GetCoefficient()
{
    std::vector<double> v(coeff, coeff + N);
    return v;
}

        
std::function<double(double)> Polynom::GetFunction()
{
    return (std::function<double(double)>)std::bind(&Polynom::evaluate, 
                                                    this, 
                                                    std::placeholders::_1);
}

namespace PROPOSAL {
    std::ostream& operator<<(std::ostream& os, const Polynom& p)
    {
        bool flag = false;
        for (int i = 0; i < p.N; ++i) {
            if(p.coeff[i] != 0)
            {
                if(flag==true)
                    os << " + " << p.coeff[i] << " * x^{" << i << "}";
                else {
                    os << "p(x) = " << p.coeff[i] << " * x^{" << i << "}";
                    flag = true;
                }
            }
        }
        return os;
    }
} // namespace PROPOSAL 
