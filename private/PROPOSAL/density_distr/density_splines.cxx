#include "PROPOSAL/density_distr/density_splines.h" 
#include "PROPOSAL/math/MathMethods.h" 
#include <functional>
#include <algorithm>
#include <iostream>

Density_splines::Density_splines(const Axis& axis, 
                                 const std::vector<SplineCoefficients>& splines) :
    Density_distr(axis)
{
    std::pair<double, std::vector<double>> spline_container;
    
    for (auto const& spline : splines) {
        spline_container = spline.GetSpline();
        Polynom poly(spline_container.second);
        dens_polynom_.push_back(new Density_polynomial(axis, poly));
        definition_area_.push_back(spline_container.first);
    }
}

Density_splines::Density_splines(const Density_splines& dens_splines):
    Density_distr( dens_splines ),
    dens_polynom_( dens_splines.dens_polynom_ ),
    definition_area_( dens_splines.definition_area_ )
{
}

double Density_splines::Helper_function(Vector3D xi, 
                                           Vector3D direction, 
                                           double res, 
                                           double l) const 
{
    return Integrate(xi, direction, l) - Integrate(xi, direction, 0) - res;
}

double Density_splines::helper_function(Vector3D xi, 
                                           Vector3D direction, 
                                           double res, 
                                           double l) const 
{
    double depth = axis_->GetDepth(xi) + l * 
                   axis_->GetEffectiveDistance(xi, direction);
    
    unsigned int loss_in_nth_spline = 0;
    for (unsigned int i = 0; i < definition_area_.size() ; ++i) 
    {
        if(depth < definition_area_[i])
            break;
        loss_in_nth_spline = i;
    }

    return dens_polynom_[loss_in_nth_spline]->GetCorrection(xi + l * direction);
}

double Density_splines::Correct(Vector3D xi, 
                                Vector3D direction,
                                double res) const 
{
    std::function<double(double)> F = std::bind(&Density_splines::Helper_function, 
                                                this, 
                                                xi, 
                                                direction, 
                                                res, 
                                                std::placeholders::_1);
    
    std::function<double(double)> dF = std::bind(&Density_splines::helper_function, 
                                                this, 
                                                xi, 
                                                direction, 
                                                res,
                                                std::placeholders::_1);

    res = NewtonRaphson(F, dF, 0, 1e15, 1.e-6);

    return res;
}

double Density_splines::Integrate(Vector3D xi, 
                                  Vector3D direction, 
                                  double l) const
{
    double depth = axis_->GetDepth(xi) + l * 
                   axis_->GetEffectiveDistance(xi, direction);

    unsigned int loss_in_nth_spline = 0;
    for (unsigned int i = 0; i < definition_area_.size() ; ++i) 
    {
        if(depth < definition_area_[i])
            break;
        loss_in_nth_spline = i;
    }

    return dens_polynom_[loss_in_nth_spline]->Integrate(xi, direction, l + definition_area_[loss_in_nth_spline]);
}

double Density_splines::Calculate(Vector3D xi, 
                                     Vector3D direction, 
                                     double distance) const
{
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_splines::GetCorrection(Vector3D xi) const
{
    double depth = axis_->GetDepth(xi);
    unsigned int loss_in_nth_spline = 0;

    for (unsigned int i = 0; i < definition_area_.size() ; ++i) 
    {
        if(depth < definition_area_[i])
            break;
        loss_in_nth_spline = i;
    }
    return dens_polynom_[loss_in_nth_spline]->GetCorrection(xi);
}

