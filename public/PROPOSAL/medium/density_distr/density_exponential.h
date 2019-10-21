#pragma once

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/medium/density_distr/density_distr.h"

namespace PROPOSAL {
class Density_exponential : public Density_distr {
   public:
    Density_exponential(const Axis& axis, double sigma);

    ~Density_exponential(){};

    bool compare(const Density_distr& dens_distr) const override;

    Density_distr* clone() const override {
        return new Density_exponential(*this);
    };

    double Integrate(const Vector3D& xi,
                     const Vector3D& direction,
                     double res) const override;
    double Correct(const Vector3D& xi,
                   const Vector3D& direction,
                   double res,
                   double distance_to_border) const override;

    double Calculate(const Vector3D& xi,
                     const Vector3D& direction,
                     double distance) const override;
    double Evaluate(const Vector3D& xi) const override;

    double GetDepth(const Vector3D& xi) const;
    double GetEffectiveDistance(const Vector3D& xi, const Vector3D& direction) const;

   private:
    double sigma_;
};
}  // namespace PROPOSAL
