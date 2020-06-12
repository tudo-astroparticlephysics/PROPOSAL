#include "PROPOSAL/secondaries/photopairproduction/PhotopairProduction.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/Constants.h"

using PROPOSAL::Components::Component;

namespace PROPOSAL {
namespace secondaries {
    class PhotoTsai : public PhotopairProduction {
        Integral integral;

        double FunctionToIntegral(double energy, double x, double theta, const Component&);

    public:
        static constexpr int n_rnd = 3;

        PhotoTsai() = default;

        double CalculateRho(double, double) override;
        tuple<Vector3D, Vector3D> CalculateDirections(
            Vector3D, double, double, const Component&, array<double, 3>) override;
        tuple<double, double> CalculateEnergy(double, double) override;
    };
} // namespace secondaries
} // namespace PROPOSAL
