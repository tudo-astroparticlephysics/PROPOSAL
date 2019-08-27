#pragma once
#include <functional>
#include <iostream>
#include <vector>
#include "PROPOSAL/math/Function.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/medium/density_distr/density_distr.h"

class Density_polynomial : public Density_distr {
   public:
    Density_polynomial(const Axis&, const Polynom&);
    Density_polynomial(const Density_polynomial&);

    Density_polynomial* clone() const { return new Density_polynomial(*this); };

    double Correct(Vector3D xi,
                   Vector3D direction,
                   double res,
                   double distance_to_border) const override;
    double Integrate(Vector3D xi, Vector3D direction, double l) const override;
    double Evaluate(Vector3D xi, Vector3D direction, double l) const;
    double Calculate(Vector3D xi,
                     Vector3D direction,
                     double distance) const override;

    double Helper_function(Vector3D xi,
                           Vector3D direction,
                           double res,
                           double l) const;
    double helper_function(Vector3D xi,
                           Vector3D direction,
                           double res,
                           double l) const;

   protected:
    Polynom polynom_;
    Polynom Polynom_;

    std::function<double(double)> density_distribution;
    std::function<double(double)> antiderived_density_distribution;
};

