
#include <fstream>

#include "gtest/gtest.h"

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/math/Spherical3D.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include "PROPOSAL/scattering/multiple_scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/multiple_scattering/Highland.h"
#include "PROPOSAL/scattering/multiple_scattering/HighlandIntegral.h"
#include "PROPOSAL/scattering/multiple_scattering/Moliere.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"
#include "PROPOSAL/crosssection/CrossSection.h"

#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"

#include "PROPOSALTestUtilities/TestFilesHandling.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus")
    {
        return MuMinusDef();
    } else if (name == "TauMinus")
    {
        return TauMinusDef();
    } else
    {
        return EMinusDef();
    }
}

auto GetCrossSections(const ParticleDef& p_def, const Medium& med, std::shared_ptr<const EnergyCutSettings> cuts, bool interpolate) {
    // old TestFiles were created using the old StandardCrossSections that are recreated here
    crosssection_list_t cross;

    auto brems = crosssection::BremsKelnerKokoulinPetrukhin{ true, p_def, med };
    cross.push_back(make_crosssection(brems, p_def, med, cuts, interpolate));
    auto epair = crosssection::EpairKelnerKokoulinPetrukhin{ true, p_def, med };
    cross.push_back(make_crosssection(epair, p_def, med, cuts, interpolate));
    auto ioniz = crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cuts)};
    cross.push_back(make_crosssection(ioniz, p_def, med, cuts, interpolate));
    auto photo = crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikheyev>() };
    cross.push_back(make_crosssection(photo, p_def, med, cuts, interpolate));

    return cross;
}

TEST(Comparison, Comparison_equal)
{
    ParticleDef mu = MuMinusDef();
    auto water = Water();

    multiple_scattering::Parametrization* moliere1 = new multiple_scattering::Moliere(mu, water);
    multiple_scattering::Moliere moliere2(mu, water);

    EXPECT_TRUE(*moliere1 == moliere2);

    multiple_scattering::Parametrization* high1 = new multiple_scattering::Highland(mu, water);
    multiple_scattering::Highland high2(mu, water);

    EXPECT_TRUE(*high1 == high2);

    // TODO: Add ScatteringHighlandIntegral as soon as it gets a compare operator
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef mu = MuMinusDef();
    auto water = Water();

    multiple_scattering::Moliere moliere1(mu, water);
    multiple_scattering::Moliere moliere2 = moliere1;

    EXPECT_TRUE(moliere1 == moliere2);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef mu = MuMinusDef();
    auto water = Water();

    multiple_scattering::Highland moliere1(mu, water);
    multiple_scattering::Highland moliere2(moliere1);
    EXPECT_TRUE(moliere1 == moliere2);
}


// Tests for virtual Scattering class

class ScatterDummy : public multiple_scattering::Parametrization {
public:
    ScatterDummy(const ParticleDef& p_def) : Parametrization(p_def.mass){}
    double GetMass(){return mass;}

    void SetOffset(std::pair<double, double> offset){theta_offset_ = offset.first; phi_offset_ = offset.second;}
    void SetScatteringAngle(std::pair<double, double> scatter){theta_scatter_ = scatter.first; phi_scatter_ = scatter.second;}

    std::unique_ptr<Parametrization> clone() const final
    {
        return std::unique_ptr<Parametrization>(
                std::make_unique<ScatterDummy>(*this));
    }

    multiple_scattering::ScatteringOffset CalculateRandomAngle(double, double,
                                double, const std::array<double, 4>&) final {
        multiple_scattering::ScatteringOffset random_angles;
        //Set offset
        random_angles.sx = std::sin(theta_offset_) * std::cos(phi_offset_);
        random_angles.sy = std::sin(theta_offset_) * std::sin(phi_offset_);
        //Set Scattering
        random_angles.tx = std::sin(theta_scatter_) * std::cos(phi_scatter_);
        random_angles.ty = std::sin(theta_scatter_) * std::sin(phi_scatter_);

        return random_angles;
    };

    bool compare(const Parametrization&) const override {return false;};
    void print(std::ostream&) const override {};

    double theta_offset_; //polar angle for offset
    double phi_offset_; //azimuthal angle for offset

    double theta_scatter_; //polar angle for scattering
    double phi_scatter_; //azimithal angle for scattering
};

TEST(Scattering, Constructor){
    ScatterDummy* dummy1 = new ScatterDummy(MuMinusDef());
    ScatterDummy* dummy2 = new ScatterDummy(TauMinusDef());

    ASSERT_DOUBLE_EQ(dummy1->GetMass(), MMU);
    ASSERT_DOUBLE_EQ(dummy2->GetMass(), MTAU);
}

TEST(Scattering, Scatter){
    std::array<std::pair<double, double>, 5> pairs{{{0,0}, {PI/4, PI/4},
                                                    {PI/2, PI/2.}, {PI/4, -PI/2.}}} ;

    std::vector<Cartesian3D> direction_list;
    direction_list.emplace_back(1.,0.,0.);
    direction_list.emplace_back(0.,-1.,0.);
    direction_list.emplace_back(1./SQRT2,0.,-1/SQRT2);

    ScatterDummy* dummy3 = new ScatterDummy(MuMinusDef());
    auto position_init = Cartesian3D(0, 0, 0);

    for(const std::pair<double, double> & offset: pairs) {
        for (const std::pair<double, double> &scatter: pairs) {
            for(auto direction_init: direction_list) {
                dummy3->SetOffset(offset);
                dummy3->SetScatteringAngle(scatter);

                auto coords = dummy3->CalculateRandomAngle(1, 10., 1., {0, 0, 0, 0});
                auto directions = multiple_scattering::ScatterInitialDirection(
                        direction_init, coords);
                // Expect directions to be normalized
                ASSERT_DOUBLE_EQ(std::get<0>(directions).magnitude(), 1.);
                ASSERT_DOUBLE_EQ(std::get<1>(directions).magnitude(), 1.);
                // Expect input polar angles to be equal to angles between input direction and output directions
                ASSERT_NEAR(std::acos(std::min((direction_init * std::get<0>(directions)), 1.)), offset.first,
                            offset.first * 1e-10);
                ASSERT_NEAR(std::acos(std::min((direction_init * std::get<1>(directions)), 1.)), scatter.first,
                            scatter.first * 1e-10);
            }
        }
    }

}

TEST(Scattering, BorderCases){
    auto medium = StandardRock();
    auto position_init  = Cartesian3D(0, 0, 0);
    auto direction_init = Cartesian3D(1, 0, 0);

    std::array<multiple_scattering::Parametrization*, 2> scatter_list = {new multiple_scattering::Moliere(MuMinusDef(), medium),
                                               new multiple_scattering::Highland(MuMinusDef(), medium)};

    for(multiple_scattering::Parametrization* scatter : scatter_list) {
        // Expect no change of direction for displacement of almost zero
        auto coords = scatter->CalculateRandomAngle(1e-20, 1e4, 1e3, {0.1, 0.2, 0.3, 0.4});
        auto directions = multiple_scattering::ScatterInitialDirection(direction_init, coords);

        EXPECT_NEAR((std::get<0>(directions) - direction_init).magnitude(), 0, 1e-10);
        EXPECT_NEAR((std::get<1>(directions) - direction_init).magnitude(), 0, 1e-10);

        // Expect change of direction
        coords = scatter->CalculateRandomAngle(1000, 1e4, 1e3, {0.1, 0.2, 0.3, 0.4});
        directions = multiple_scattering::ScatterInitialDirection(direction_init, coords);

        EXPECT_FALSE(std::get<0>(directions) == direction_init);
        EXPECT_FALSE(std::get<1>(directions) == direction_init);
    }
}

TEST(Scattering, FirstMomentum){
    RandomGenerator::Get().SetSeed(24601);
    auto medium = StandardRock();
    auto position_init  = Cartesian3D(0, 0, 0);
    auto direction_init = Cartesian3D(0, 0, 1);
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, false);

    int statistics = 1e3;
    auto cross = GetCrossSections(MuMinusDef(), medium, cuts, true);

    std::array<std::unique_ptr<multiple_scattering::Parametrization>, 3> scatter_list = {make_multiple_scattering("moliere", MuMinusDef(), medium),
                                                               make_multiple_scattering("highland", MuMinusDef(), medium),
                                                               make_multiple_scattering("highlandintegral", MuMinusDef(), medium, cross)};
    Cartesian3D scatter_sum;
    Cartesian3D offset_sum;

    for(auto const& scatter: scatter_list){
        scatter_sum.SetCoordinates({0., 0., 0.});
        offset_sum.SetCoordinates({0., 0., 0.});

        for (int n=1; n<=statistics; ++n) {
            auto coords = scatter->CalculateRandomAngle(
                    1e3, 1e4, 1e3,
                    {RandomGenerator::Get().RandomDouble(), RandomGenerator::Get().RandomDouble(),
                     RandomGenerator::Get().RandomDouble(), RandomGenerator::Get().RandomDouble()}
            );

            auto sampled_vectors = multiple_scattering::ScatterInitialDirection(
                    direction_init, coords);

            offset_sum = offset_sum + ( std::get<0>(sampled_vectors) - offset_sum) * (1./n);
            scatter_sum = scatter_sum + ( std::get<1>(sampled_vectors) - scatter_sum)* (1./n);
        }

        EXPECT_NEAR((offset_sum.GetX() - direction_init.GetX()), 0., 1e-3);
        EXPECT_NEAR((offset_sum.GetY() - direction_init.GetY()), 0., 1e-3);

        EXPECT_NEAR((scatter_sum.GetX() - scatter_sum.GetX()), 0., 1e-3);
        EXPECT_NEAR((scatter_sum.GetY() - scatter_sum.GetY()), 0., 1e-3);
    }
}

TEST(Scattering, SecondMomentum){
    RandomGenerator::Get().SetSeed(24601);
    auto medium = StandardRock();
    auto position_init  = Cartesian3D(0, 0, 0);
    auto direction_init = Cartesian3D(0, 0, 1);
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, false);

    int statistics = 1e3;

    double E_i = 1e14;
    std::array<double, 6> final_energies = {1e13, 1e11, 1e9, 1e7, 1e5, 1e3};
    auto cross = GetCrossSections(MuMinusDef(), medium, cuts, true);

    std::array<std::unique_ptr<multiple_scattering::Parametrization>, 3> scatter_list = {make_multiple_scattering("moliere", MuMinusDef(), medium),
                                                               make_multiple_scattering("highland", MuMinusDef(), medium),
                                                               make_multiple_scattering("highlandintegral", MuMinusDef(), medium, cross)};
    double scatter_sum;
    double offset_sum;
    double displacement;
    double old_displacement;
    auto displacement_calculator = DisplacementBuilder(cross, std::false_type());
    std::array<double, 2> variances = {0,0};
    std::array<double, 2> old_variances;

    for(auto const& scatter: scatter_list){
        old_variances = {0, 0};
        old_displacement = 0;
        for(double E_f: final_energies) {
            scatter_sum = 0;
            offset_sum = 0;
            displacement = displacement_calculator.SolveTrackIntegral(E_i, E_f);
            EXPECT_GT(displacement, old_displacement);
            for (int n = 1; n <= statistics; ++n) {
                auto coords = scatter->CalculateRandomAngle(
                        displacement, E_i, E_f,
                        {RandomGenerator::Get().RandomDouble(), RandomGenerator::Get().RandomDouble(),
                         RandomGenerator::Get().RandomDouble(), RandomGenerator::Get().RandomDouble()}
                );
                auto sampled_vectors = multiple_scattering::ScatterInitialDirection(
                        direction_init, coords);
                offset_sum = offset_sum + ((std::get<0>(sampled_vectors) - direction_init).magnitude() - offset_sum) * (1. / n);
                scatter_sum = scatter_sum + ((std::get<1>(sampled_vectors) - direction_init).magnitude() - scatter_sum) * (1. / n);
            }

            variances = {offset_sum, scatter_sum};
            EXPECT_GT(variances, old_variances);

            old_displacement = displacement;
            old_variances = variances;
        }
    }
}

TEST(Scattering, compare_integral_interpolant) {
    RandomGenerator::Get().SetSeed(24601);
    auto medium = StandardRock();
    auto position_init  = Cartesian3D(0, 0, 0);
    auto direction_init = Cartesian3D(0, 0, 1);

    std::vector<ParticleDef> particles = {EMinusDef(), MuMinusDef()};
    auto cut = std::make_shared<EnergyCutSettings>(INF, 1, false);

    for (auto p : particles) {
        auto cross = GetCrossSections(p, Ice(), cut, true);
        auto scatter_integral = make_multiple_scattering("highlandintegral", p, medium, cross, false);
        auto scatter_interpol = make_multiple_scattering("highlandintegral", p, medium, cross, true);
        auto energies = std::array<double, 5>{1e6, 1e7, 1e8, 1e9, 1e10};
        for (auto E_i : energies) {
            auto rnd = std::array<double, 4>{RandomGenerator::Get().RandomDouble(),
                                             RandomGenerator::Get().RandomDouble(),
                                             RandomGenerator::Get().RandomDouble(),
                                             RandomGenerator::Get().RandomDouble()};
            auto coords_integral = scatter_integral->CalculateRandomAngle(1e4, E_i, 1e5, rnd);
            auto coords_interpol = scatter_interpol->CalculateRandomAngle(1e4, E_i, 1e5, rnd);
            auto vec_integral = multiple_scattering::ScatterInitialDirection(
                    direction_init, coords_integral);
            auto vec_interpol = multiple_scattering::ScatterInitialDirection(
                    direction_init, coords_interpol);
            EXPECT_NEAR(std::get<0>(vec_integral).GetX(), std::get<0>(vec_interpol).GetX(),
                        std::abs(std::get<0>(vec_integral).GetX() * 1e-3));
            EXPECT_NEAR(std::get<0>(vec_integral).GetY(), std::get<0>(vec_interpol).GetY(),
                        std::abs(std::get<0>(vec_integral).GetY() * 1e-3));
            EXPECT_NEAR(std::get<0>(vec_integral).GetZ(), std::get<0>(vec_interpol).GetZ(),
                        std::abs(std::get<0>(vec_integral).GetZ() * 1e-3));
            EXPECT_NEAR(std::get<1>(vec_integral).GetX(), std::get<1>(vec_interpol).GetX(),
                        std::abs(std::get<1>(vec_integral).GetX() * 1e-3));
            EXPECT_NEAR(std::get<1>(vec_integral).GetY(), std::get<1>(vec_interpol).GetY(),
                        std::abs(std::get<1>(vec_integral).GetY() * 1e-3));
            EXPECT_NEAR(std::get<1>(vec_integral).GetZ(), std::get<1>(vec_interpol).GetZ(),
                        std::abs(std::get<1>(vec_integral).GetZ() * 1e-3));
        }
    }
}

TEST(Scattering, ScatterReproducibilityTest)
{
    auto in = getTestFiles("Scattering_scatter.txt");

    std::string particleName;
    std::string mediumName;
    std::string parametrization;

    double energy_init, energy_final, distance;
    double rnd1, rnd2, rnd3, rnd4;
    double energy_previous = -1;
    double ecut, vcut;
    Cartesian3D position_init  = Cartesian3D(0, 0, 0);
    Cartesian3D direction_init = Cartesian3D(1, 0, 0);
    Cartesian3D position_out;
    Cartesian3D direction_out;
    double x_f, y_f, z_f;
    double radius_f, phi_f, theta_f;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    double error    = 1e-3;

    while (in >> particleName >> mediumName >> parametrization >> ecut >> vcut >> energy_init >> energy_final >> distance >> rnd1 >> rnd2 >> rnd3 >> rnd4 >> x_f >> y_f >> z_f >> radius_f >> phi_f >> theta_f)
    {
        energy_previous = -1;

        ParticleDef particle_def = getParticleDef(particleName);

        std::shared_ptr<const Medium> medium = CreateMedium(mediumName);

        //reprouce old behaviour
        if(ecut==-1){
            ecut = std::numeric_limits<double>::infinity();
        }
        if(vcut==-1){
            vcut = 1;
        }

        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, false);

        crosssection_list_t cross;

        std::unique_ptr<multiple_scattering::Parametrization> scattering = NULL;
        if (parametrization == "NoScattering")
        {
            continue; // not implemented anymore
        } else if (parametrization == "HighlandIntegral") {
            cross = GetCrossSections(particle_def, *medium, ecuts, false);
        }

        scattering = make_multiple_scattering(parametrization, particle_def, *medium, cross, false);

        // There has been a correction in the LPM effect parametrization for
        // bremsstrahlung which influences the scattering angles for electrons
        // see commit 7be271c3eeafc8b7093340168c2cd739392ee6c4
        if (particleName == "EMinus" and parametrization == "HighlandIntegral")
            continue;

        while (energy_previous < energy_init)
        {
            energy_previous = energy_init;
            if(energy_final > particle_def.mass) {
                std::array<double, 4> rnd{rnd1, rnd2, rnd3, rnd4};
                auto coords = scattering->CalculateRandomAngle(distance * medium->GetMassDensity(),
                                                      energy_init,
                                                      energy_final,
                                                      rnd);
                auto directions = multiple_scattering::ScatterInitialDirection(
                        direction_init, coords);
                position_out = position_init + distance * std::get<0>(directions);
                direction_out = std::get<1>(directions);

                EXPECT_NEAR(position_out.GetX(), x_f, std::abs(error * x_f));
                EXPECT_NEAR(position_out.GetY(), y_f, std::abs(error * y_f));
                EXPECT_NEAR(position_out.GetZ(), z_f, std::abs(error * z_f));

                auto direction_out_spherical = Spherical3D(direction_out);
                EXPECT_NEAR(direction_out_spherical.GetRadius(), radius_f, std::abs(error * radius_f));
                EXPECT_NEAR(direction_out_spherical.GetAzimuth(), phi_f, std::abs(error * phi_f));
                EXPECT_NEAR(direction_out_spherical.GetZenith(), theta_f, std::abs(error * theta_f));
            }
            in >> particleName >> mediumName >> parametrization >> ecut >> vcut >> energy_init >> energy_final >>
                distance >> rnd1 >> rnd2 >> rnd3 >> rnd4 >> x_f >> y_f >> z_f >> radius_f >> phi_f >> theta_f;

            //reprouce old behaviour
            if(ecut==-1){
                ecut = std::numeric_limits<double>::infinity();
            }
            if(vcut==-1){
                vcut = 1;
            }

        }
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
