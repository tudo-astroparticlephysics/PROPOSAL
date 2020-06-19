#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/secondaries/photopairproduction/PhotopairProduction.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"

using PROPOSAL::Components::Component;

namespace PROPOSAL {
namespace secondaries {
    class PhotoTsai : public PhotopairProduction ,
                                 public RegisteredInDefault<PhotoTsai> {
        Integral integral;

        double FunctionToIntegral(
            double energy, double x, double theta, const Component&);

    public:
        static constexpr int n_rnd = 3;

        PhotoTsai() = default;
        PhotoTsai(ParticleDef, Medium) {};

        double CalculateRho(double, double) override;
        tuple<Vector3D, Vector3D> CalculateDirections(Vector3D, double, double,
            const Component&, vector<double>) override;
        tuple<double, double> CalculateEnergy(double, double, double ) override;

        size_t RequiredRandomNumbers() final { return n_rnd; }
        vector<Loss::secondary_t> CalculateSecondaries(
            double, Loss::secondary_t, const Component&, vector<double>) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
