#include "PROPOSAL/density_distr/density_polynomial.h" 
#include <functional>
#include <algorithm>
#include <iostream>


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Polynom      %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Polynom::Polynom(std::vector<double> coefficients)
{
    N = coefficients.size();
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

    for (int i = 0; i < N; ++i) 
        aux += coeff[i] * std::pow(x, i);

    return aux;
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

        
std::function<double(double)> Polynom::GetFunction()
{
    return (std::function<double(double)>)std::bind(&Polynom::evaluate, 
                                                    this, 
                                                    std::placeholders::_1);
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Density_polynomial::Density_polynomial(const Axis& axis, const Polynom& polynom):
    Density_distr( axis ),
    polynom_( polynom )
{
    density_distribution = polynom_.GetFunction();
    antiderived_density_distribution = polynom_.GetAntiderivative(0).GetFunction();
}

Density_polynomial::Density_polynomial(const Density_polynomial& dens):
    Density_distr( dens ),
    polynom_( dens.polynom_ )
{
    density_distribution = polynom_.GetFunction();
    antiderived_density_distribution = polynom_.GetAntiderivative(0).GetFunction();
}


double Density_polynomial::Correct(Vector3D xi, Vector3D direction, double res) const 
{
    double f(double l)
    {
        return antiderived_density_distribution(GetDepth(xi) + l * delta) 
            - antiderived_density_distribution(GetDepth(xi))
            - res;
    }
    std::function<double(double)> = f;

    return res;
}

double Density_polynomial::Integrate(Vector3D xi, Vector3D direction, double l) const
{
    double delta = GetEffectiveDistance(xi, direction);

    return antiderived_density_distribution(GetDepth(xi) + l * delta);
}

double Density_polynomial::Calculate(Vector3D xi, Vector3D direction, double distance) const
{
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}
