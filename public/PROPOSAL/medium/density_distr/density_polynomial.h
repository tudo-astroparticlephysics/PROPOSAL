#pragma once
#include <functional>
#include <iostream>
#include <vector>
#include "PROPOSAL/math/Function.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/medium/density_distr/density_distr.h"

namespace PROPOSAL {
class Density_polynomial : public Density_distr {
   public:
    Density_polynomial(const Axis&, const Polynom&);
    Density_polynomial(const Density_polynomial&);
    ~Density_polynomial();

    bool compare(const Density_distr& dens_distr) const override;

    Density_distr* clone() const override {
        return new Density_polynomial(*this);
    };

    double Correct(const Vector3D& xi,
                   const Vector3D& direction,
                   double res,
                   double distance_to_border) const override;
    double Integrate(const Vector3D& xi, const Vector3D& direction, double l) const override;
    double Evaluate(const Vector3D& xi) const override;
    double Calculate(const Vector3D& xi,
                     const Vector3D& direction,
                     double distance) const override;

    double Helper_function(const Vector3D& xi,
                           const Vector3D& direction,
                           double res,
                           double l) const;
    double helper_function(const Vector3D& xi,
                           const Vector3D& direction,
                           double res,
                           double l) const;

   protected:
    Polynom polynom_;
    Polynom Polynom_;

    std::function<double(double)> density_distribution;
    std::function<double(double)> antiderived_density_distribution;
};
}  // namespace PROPOSAL
