#pragma once
#include <functional>
#include <iostream>
#include <vector>
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/Spline.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/medium/density_distr/density_polynomial.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class Density_splines : public Density_distr {
   public:
    Density_splines(const Axis&, const Spline&);
    Density_splines(const Density_splines&);

    Density_splines* clone() const { return new Density_splines(*this); };

    double Correct(Vector3D xi,
                   Vector3D direction,
                   double res,
                   double distance_to_border) const override;
    double Integrate(Vector3D xi, Vector3D direction, double l) const override;
    double Evaluate(Vector3D xi, Vector3D direction, double l) const;
    double Calculate(Vector3D xi,
                     Vector3D direction,
                     double distance) const override;
    // double GetCorrection(Vector3D xi) const override;

    double Helper_function(Vector3D xi,
                           Vector3D direction,
                           double res,
                           double l) const;
    double helper_function(Vector3D xi,
                           Vector3D direction,
                           double res,
                           double l) const;

   protected:
    Spline* spline_;
    Spline* integrated_spline_;

    std::function<double(double)> density_distribution;
    std::function<double(double)> antiderived_density_distribution;
};

