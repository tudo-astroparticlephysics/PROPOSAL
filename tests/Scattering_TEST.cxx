
#include <fstream>

#include "gtest/gtest.h"

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Vector3D.h"

#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"

#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"

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
    crosssection_list_t<ParticleDef, Medium> cross;

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

    Scattering* moliere1 = new ScatteringMoliere(mu, water);
    ScatteringMoliere moliere2(mu, water);

    EXPECT_TRUE(*moliere1 == moliere2);

    Scattering* high1 = new ScatteringHighland(mu, water);
    ScatteringHighland high2(mu, water);

    EXPECT_TRUE(*high1 == high2);

    // TODO: Add ScatteringHighlandIntegral as soon as it gets a compare operator
}

TEST(Comparison, Comparison_not_equal)
{
    ParticleDef mu  = MuMinusDef();
    ParticleDef tau = TauMinusDef();
    auto water = Water();
    auto ice = Ice();

    ScatteringMoliere moliere1(mu, water);
    ScatteringMoliere moliere2(tau, water);
    ScatteringMoliere moliere3(mu, ice);

    EXPECT_TRUE(moliere1 != moliere2);
    EXPECT_TRUE(moliere1 != moliere3);

    ScatteringHighland high1(mu, water);
    ScatteringHighland high2(tau, water);
    ScatteringHighland high3(mu, ice);

    EXPECT_TRUE(high1 != high2);
    EXPECT_TRUE(high1 != high3);

    // TODO: Add ScatteringHighlandIntegral as soon as it gets a compare operator
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef mu = MuMinusDef();
    auto water = Water();

    ScatteringMoliere moliere1(mu, water);
    ScatteringMoliere moliere2 = moliere1;
    EXPECT_TRUE(moliere1 == moliere2);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef mu = MuMinusDef();
    auto water = Water();

    ScatteringHighland moliere1(mu, water);
    ScatteringHighland moliere2(moliere1);
    EXPECT_TRUE(moliere1 == moliere2);
}


// Tests for virtual Scattering class

class ScatterDummy : public Scattering{
public:
    ScatterDummy(const ParticleDef& p_def) : Scattering(p_def){}
    double GetMass(){return mass;}

    void SetOffset(std::pair<double, double> offset){theta_offset_ = offset.first; phi_offset_ = offset.second;}
    void SetScatteringAngle(std::pair<double, double> scatter){theta_scatter_ = scatter.first; phi_scatter_ = scatter.second;}

private:
    RandomAngles CalculateRandomAngle(double grammage,
                                      double ei,
                                      double ef,
                                      const std::array<double, 4>& rnd){
        Scattering::RandomAngles random_angles;
        //Set offset
        random_angles.sx = std::sin(theta_offset_) * std::cos(phi_offset_);
        random_angles.sy = std::sin(theta_offset_) * std::sin(phi_offset_);
        //Set Scattering
        random_angles.tx = std::sin(theta_scatter_) * std::cos(phi_scatter_);
        random_angles.ty = std::sin(theta_scatter_) * std::sin(phi_scatter_);

        return random_angles;
    };

    bool compare(const Scattering&) const{return false;};
    void print(std::ostream&) const{};

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

    std::vector<Vector3D> direction_list;
    direction_list.emplace_back(1.,0.,0.);
    direction_list.emplace_back(0.,-1.,0.);
    direction_list.emplace_back(1./SQRT2,0.,-1/SQRT2);

    ScatterDummy* dummy3 = new ScatterDummy(MuMinusDef());
    Vector3D position_init  = Vector3D(0, 0, 0);

    for(const std::pair<double, double> & offset: pairs) {
        for (const std::pair<double, double> &scatter: pairs) {
            for(Vector3D &direction_init: direction_list) {
                direction_init.CalculateSphericalCoordinates();
                dummy3->SetOffset(offset);
                dummy3->SetScatteringAngle(scatter);
                auto directions = dummy3->Scatter(1, 10., 1., direction_init, {0, 0, 0, 0});
                // Expect directions to be normalized
                ASSERT_DOUBLE_EQ(std::get<0>(directions).magnitude(), 1.);
                ASSERT_DOUBLE_EQ(std::get<1>(directions).magnitude(), 1.);
                // Expect input polar angles to be equal to angles between input direction and output directions
                ASSERT_NEAR(std::acos(std::min(scalar_product(direction_init, std::get<0>(directions)), 1.)), offset.first,
                            offset.first * 1e-10);
                ASSERT_NEAR(std::acos(std::min(scalar_product(direction_init, std::get<1>(directions)), 1.)), scatter.first,
                            scatter.first * 1e-10);
            }
        }
    }

}

TEST(Scattering, BorderCases){
    auto medium = StandardRock();
    Vector3D position_init  = Vector3D(0, 0, 0);
    Vector3D direction_init = Vector3D(1, 0, 0);
    direction_init.CalculateSphericalCoordinates();

    std::array<Scattering*, 2> scatter_list = {new ScatteringMoliere(MuMinusDef(), medium),
                                               new ScatteringHighland(MuMinusDef(), medium)};

    for(Scattering* scatter : scatter_list) {
        // Expect no change of direction for displacement of almost zero
        EXPECT_NEAR((std::get<0>(scatter->Scatter(1e-20, 1e4, 1e3, direction_init, {0.1, 0.2, 0.3, 0.4})) - direction_init).magnitude(),
                    0, 1e-10);
        EXPECT_NEAR((std::get<1>(scatter->Scatter(1e-20, 1e4, 1e3, direction_init, {0.1, 0.2, 0.3, 0.4})) - direction_init).magnitude(),
                    0, 1e-10);

        // Expect change of direction
        EXPECT_FALSE(std::get<0>(scatter->Scatter(1000, 1e4, 1e3, direction_init, {0.1, 0.2, 0.3, 0.4})) ==
                direction_init);
        EXPECT_FALSE(std::get<1>(scatter->Scatter(1000, 1e4, 1e3, direction_init, {0.1, 0.2, 0.3, 0.4})) ==
                direction_init);

    }
}

TEST(Scattering, FirstMomentum){
    RandomGenerator::Get().SetSeed(24601);
    auto medium = StandardRock();
    Vector3D position_init  = Vector3D(0, 0, 0);
    Vector3D direction_init = Vector3D(0, 0, 1);
    direction_init.CalculateSphericalCoordinates();
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, false);

    int statistics = 1e7;
    auto cross = GetCrossSections(MuMinusDef(), medium, cuts, true);

    std::array<std::unique_ptr<Scattering>, 3> scatter_list = {make_scattering("moliere", MuMinusDef(), medium),
                                                               make_scattering("highland", MuMinusDef(), medium),
                                                               make_scattering("highland_integral", MuMinusDef(), medium, cross)};
    Vector3D scatter_sum;
    Vector3D offset_sum;

    for(auto const& scatter: scatter_list){
        Vector3D scatter_sum = Vector3D(0., 0., 0.);
        Vector3D offset_sum = Vector3D(0., 0., 0.);

        for (int n=1; n<=statistics; ++n) {
            auto sampled_vectors = scatter->Scatter(
                    1e3, 1e4, 1e3, direction_init,
                    {RandomGenerator::Get().RandomDouble(), RandomGenerator::Get().RandomDouble(),
                     RandomGenerator::Get().RandomDouble(), RandomGenerator::Get().RandomDouble()}
            );
            offset_sum = offset_sum + ( std::get<0>(sampled_vectors) - offset_sum) * (1./n);
            scatter_sum = scatter_sum + ( std::get<1>(sampled_vectors) - scatter_sum)* (1./n);
        }

        EXPECT_NEAR((offset_sum.GetX() - direction_init.GetX()), 0., 1e-5);
        EXPECT_NEAR((offset_sum.GetY() - direction_init.GetY()), 0., 1e-5);

        EXPECT_NEAR((scatter_sum.GetX() - scatter_sum.GetX()), 0., 1e-5);
        EXPECT_NEAR((scatter_sum.GetY() - scatter_sum.GetY()), 0., 1e-5);
    }
}

TEST(Scattering, SecondMomentum){
    RandomGenerator::Get().SetSeed(24601);
    auto medium = StandardRock();
    Vector3D position_init  = Vector3D(0, 0, 0);
    Vector3D direction_init = Vector3D(0, 0, 1);
    direction_init.CalculateSphericalCoordinates();
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, false);

    int statistics = 1e5;

    double E_i = 1e14;
    std::array<double, 10> final_energies = {1e13, 1e12, 1e11, 1e10, 1e9, 1e8, 1e7, 1e6, 1e5, 1e4};
    auto cross = GetCrossSections(MuMinusDef(), medium, cuts, true);

    std::array<std::unique_ptr<Scattering>, 3> scatter_list = {make_scattering("moliere", MuMinusDef(), medium),
                                                               make_scattering("highland", MuMinusDef(), medium),
                                                               make_scattering("highland_integral", MuMinusDef(), medium, cross)};
    double scatter_sum;
    double offset_sum;
    double displacement;
    double old_displacement;
    auto displacement_calculator = DisplacementBuilder<UtilityIntegral>(cross);
    std::array<double, 2> variances = {0,0};
    std::array<double, 2> old_variances;

    for(auto const& scatter: scatter_list){
        old_variances = {0, 0};
        old_displacement = 0;
        for(double E_f: final_energies) {
            scatter_sum = 0;
            offset_sum = 0;
            displacement = displacement_calculator.SolveTrackIntegral(E_i, E_f);
            ASSERT_TRUE(displacement > old_displacement);
            for (int n = 1; n <= statistics; ++n) {
                auto sampled_vectors = scatter->Scatter(
                        displacement, E_i, E_f, direction_init,
                        {RandomGenerator::Get().RandomDouble(), RandomGenerator::Get().RandomDouble(),
                         RandomGenerator::Get().RandomDouble(), RandomGenerator::Get().RandomDouble()}
                );
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
    Vector3D position_init  = Vector3D(0, 0, 0);
    Vector3D direction_init = Vector3D(0, 0, 1);
    direction_init.CalculateSphericalCoordinates();

    std::vector<ParticleDef> particles = {EMinusDef(), MuMinusDef()};
    auto cut1 = std::make_shared<EnergyCutSettings>(500, 0.05, false);
    auto cut2 = std::make_shared<EnergyCutSettings>(INF, 1, false);
    std::vector<std::shared_ptr<EnergyCutSettings>> cuts = {cut1, cut2};

    for (auto p : particles) {
        for (auto cut : cuts) {

            auto cross = GetCrossSections(p, Ice(), cut, true);

            auto scatter_integral = make_scattering("highland_integral", p, medium, cross, false);
            auto scatter_interpol = make_scattering("highland_integral", p, medium, cross, true);

            auto energies = std::array<double, 5>{1e6, 1e7, 1e8, 1e9, 1e10};

            for (auto E_i : energies) {
                auto rnd = std::array<double, 4>{RandomGenerator::Get().RandomDouble(),
                                                 RandomGenerator::Get().RandomDouble(),
                                                 RandomGenerator::Get().RandomDouble(),
                                                 RandomGenerator::Get().RandomDouble()};
                auto vec_integral = scatter_integral->Scatter(1e4, E_i, 1e5, direction_init, rnd);
                auto vec_interpol = scatter_interpol->Scatter(1e4, E_i, 1e5, direction_init, rnd);

                EXPECT_NEAR(std::get<0>(vec_integral).GetX(), std::get<0>(vec_interpol).GetX(),
                            std::abs(std::get<0>(vec_integral).GetX() * 1e-5));
                EXPECT_NEAR(std::get<0>(vec_integral).GetY(), std::get<0>(vec_interpol).GetY(),
                            std::abs(std::get<0>(vec_integral).GetY() * 1e-5));
                EXPECT_NEAR(std::get<0>(vec_integral).GetZ(), std::get<0>(vec_interpol).GetZ(),
                            std::abs(std::get<0>(vec_integral).GetZ() * 1e-5));

                EXPECT_NEAR(std::get<1>(vec_integral).GetX(), std::get<1>(vec_interpol).GetX(),
                            std::abs(std::get<1>(vec_integral).GetX() * 1e-5));
                EXPECT_NEAR(std::get<1>(vec_integral).GetY(), std::get<1>(vec_interpol).GetY(),
                            std::abs(std::get<1>(vec_integral).GetY() * 1e-5));
                EXPECT_NEAR(std::get<1>(vec_integral).GetZ(), std::get<1>(vec_interpol).GetZ(),
                            std::abs(std::get<1>(vec_integral).GetZ() * 1e-5));
            }
        }
    }
}

TEST(Scattering, ScatterReproducibilityTest)
{
    std::string filename = "bin/TestFiles/Scattering_scatter.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string particleName;
    std::string mediumName;
    std::string parametrization;

    double energy_init, energy_final, distance;
    double rnd1, rnd2, rnd3, rnd4;
    double energy_previous = -1;
    double ecut, vcut;
    Vector3D position_init  = Vector3D(0, 0, 0);
    Vector3D direction_init = Vector3D(1, 0, 0);
    Vector3D position_out;
    Vector3D direction_out;
    direction_init.CalculateSphericalCoordinates();
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

        crosssection_list_t<ParticleDef, Medium> cross;

        std::unique_ptr<Scattering> scattering = NULL;
        if (parametrization == "NoScattering")
        {
            continue; // not implemented anymore
        } else if (parametrization == "HighlandIntegral") {
            parametrization = "highland_integral";
            cross = GetCrossSections(particle_def, *medium, ecuts, false);
        }

        scattering = make_scattering(parametrization, particle_def, *medium, cross, false);


        while (energy_previous < energy_init)
        {
            energy_previous = energy_init;
            if(energy_final > particle_def.mass) {
                std::array<double, 4> rnd{rnd1, rnd2, rnd3, rnd4};
                auto directions = scattering->Scatter(distance * medium->GetMassDensity(),
                                                      energy_init,
                                                      energy_final,
                                                      direction_init,
                                                      rnd);
                position_out = position_init + distance * std::get<0>(directions);
                direction_out = std::get<1>(directions);

                EXPECT_NEAR(position_out.GetX(), x_f, std::abs(error * x_f));
                EXPECT_NEAR(position_out.GetY(), y_f, std::abs(error * y_f));
                EXPECT_NEAR(position_out.GetZ(), z_f, std::abs(error * z_f));

                EXPECT_NEAR(direction_out.GetRadius(), radius_f, std::abs(error * radius_f));
                EXPECT_NEAR(direction_out.GetPhi(), phi_f, std::abs(error * phi_f));
                EXPECT_NEAR(direction_out.GetTheta(), theta_f, std::abs(error * theta_f));
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
