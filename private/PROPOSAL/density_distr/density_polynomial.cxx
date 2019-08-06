#include "PROPOSAL/density_distr/density_polynomial.h" 
#include "PROPOSAL/math/MathMethods.h" 
#include "PROPOSAL/math/Function.h"
#include <functional>
#include <algorithm>
#include <iostream>

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
    // std::cout << "Helper function(" 
    //           << l
    //           << "): "
    //           << Integrate(xi, direction, 0) - Integrate(xi, direction, l) + res
    //           << std::endl;

    return Integrate(xi, direction, l) - Integrate(xi, direction, 0) + res;
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
        res = NewtonRaphson(F, dF, 0, 1e7, 5e3);
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

