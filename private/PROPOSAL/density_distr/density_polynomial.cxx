#include "PROPOSAL/density_distr/density_polynomial.h" 
#include "PROPOSAL/math/MathMethods.h" 
#include <functional>
#include <algorithm>
#include <iostream>


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Polynom      %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace PROPOSAL;

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

namespace PROPOSAL {
    std::ostream& operator<<(std::ostream& os, const Polynom& p)
    {
        bool flag = false;
        for (int i = 0; i < p.N; ++i) {
            if(p.coeff[i] != 0)
            {
                if(flag==true)
                    os << "+ " << p.coeff[i] << " * x^{" << i << "}";
                else {
                    os << "p(x) = " << p.coeff[i] << " * x^{" << i << "}";
                    flag = true;
                }
            }
        }
        return os;
    }
} // namespace PROPOSAL 

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Density_polynomial::Density_polynomial(const Axis& axis, const Polynom& polynom):
    Density_distr( axis ),
    polynom_( polynom ),
    Polynom_( polynom_.GetAntiderivative(0) ),
    density_distribution( polynom_.GetFunction() ),
    antiderived_density_distribution( Polynom_.GetFunction() )
{ }

Density_polynomial::Density_polynomial(const Density_polynomial& dens):
    Density_distr( dens ),
    polynom_( dens.polynom_ ),
    Polynom_( dens.Polynom_ ),
    density_distribution( polynom_.GetFunction() ),
    antiderived_density_distribution( Polynom_.GetFunction() )
{ }

double Density_polynomial::Helper_function(Vector3D xi, 
                                           Vector3D direction, 
                                           double res, 
                                           double l) const 
{
    // std::cout << "Integrate(xi, direction, l): " 
    //           << Integrate(xi, direction, l)
    //           << " Integrate(xi, direction, 0): "
    //           << Integrate(xi, direction, 0) 
    //           << " res:" << ees
    //           << std::endl;
    // std::cout << "Helper function(" 
    //           << l
    //           << "): "
    //           << Integrate(xi, direction, 0) - Integrate(xi, direction, l) + res
    //           << std::endl;

    return Integrate(xi, direction, 0) - Integrate(xi, direction, l) + res;
}

double Density_polynomial::helper_function(Vector3D xi, 
                                           Vector3D direction, 
                                           double res, 
                                           double l) const 
{
    double delta = axis_->GetEffectiveDistance(xi, direction);

    return density_distribution(axis_->GetDepth(xi) + l * delta);
}

double Density_polynomial::Correct(Vector3D xi, 
                                   Vector3D direction,
                                   double res) const 
{
    std::function<double(double)> F = std::bind(&Density_polynomial::Helper_function, 
                                                this, 
                                                xi, 
                                                direction, 
                                                res, 
                                                std::placeholders::_1);
    
    std::function<double(double)> dF = std::bind(&Density_polynomial::helper_function, 
                                                this, 
                                                xi, 
                                                direction, 
                                                res,
                                                std::placeholders::_1);

    // check if direction * axis larger or less than zero 
    // direction * fAxis_

    try {
        res = NewtonRaphson(F, dF, 0, 1e15, 1.e-6);
    } catch (MathException& e) {
        double depth = GetCorrection(xi + res * direction);
        if ( depth < 0 ) 
            DensityException("Densities smaller than zero are non-physical. Check coefficients from polynom.");
        else
            throw e;
    }

    return res;
}

double Density_polynomial::Integrate(Vector3D xi, 
                                     Vector3D direction, 
                                     double l) const
{
    double delta = axis_->GetEffectiveDistance(xi, direction);

    // std::cout << "Integrate( "<<axis_->GetDepth(xi) + l * delta << "): " 
    //           << antiderived_density_distribution(axis_->GetDepth(xi) + l * delta) 
    //           << std::endl;

    return antiderived_density_distribution(axis_->GetDepth(xi) + l * delta);
}

double Density_polynomial::Calculate(Vector3D xi, 
                                     Vector3D direction, 
                                     double distance) const
{
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_polynomial::GetCorrection(Vector3D xi) const
{
    return density_distribution(axis_->GetDepth(xi));
}

