#pragma once

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/medium/density_distr/density_distr.h"

namespace PROPOSAL {
class Density_exponential : public Density_distr {
   public:
    Density_exponential(const Axis& axis, double sigma);

    ~Density_exponential(){};

    Density_distr* clone() const override {
        return new Density_exponential(*this);
    };

    double Integrate(Vector3D xi,
                     Vector3D direction,
                     double res) const override;
    double Correct(Vector3D xi,
                   Vector3D direction,
                   double res,
                   double distance_to_border) const override;

    double Calculate(Vector3D xi,
                     Vector3D direction,
                     double distance) const override;
    double Evaluate(Vector3D xi) const override;

    double GetDepth(Vector3D xi) const;
    double GetEffectiveDistance(Vector3D xi, Vector3D direction) const;

   private:
    double sigma_;
};
}  // namespace PROPOSAL
