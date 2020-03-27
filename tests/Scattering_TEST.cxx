
#include <fstream>

#include "gtest/gtest.h"

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/math/RandomGenerator.h"

#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"

#include "PROPOSAL/crossection/factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crossection/factories/IonizationFactory.h"
#include "PROPOSAL/crossection/factories/EpairProductionFactory.h"
#include "PROPOSAL/crossection/factories/PhotonuclearFactory.h"

#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus")
    {
        return MuMinusDef::Get();
    } else if (name == "TauMinus")
    {
        return TauMinusDef::Get();
    } else
    {
        return EMinusDef::Get();
    }
}

TEST(Comparison, Comparison_equal)
{
    ParticleDef mu = MuMinusDef::Get();
    std::shared_ptr<const Medium> water(new Water(1.0));

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
    ParticleDef mu  = MuMinusDef::Get();
    ParticleDef tau = TauMinusDef::Get();
    std::shared_ptr<const Medium> water(new Water(1.0));
    std::shared_ptr<const Medium> ice(new Ice());


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
    ParticleDef mu = MuMinusDef::Get();
    std::shared_ptr<const Medium> water(new Water(1.0));

    ScatteringMoliere moliere1(mu, water);
    ScatteringMoliere moliere2 = moliere1;
    EXPECT_TRUE(moliere1 == moliere2);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef mu = MuMinusDef::Get();
    std::shared_ptr<const Medium> water(new Water(1.0));

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
    RandomAngles CalculateRandomAngle(double dr,
                                      double ei,
                                      double ef,
                                      const Vector3D& pos,
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
    ScatterDummy* dummy1 = new ScatterDummy(getParticleDef("MuMinus"));
    ScatterDummy* dummy2 = new ScatterDummy(getParticleDef("TauMinus"));

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

    ScatterDummy* dummy3 = new ScatterDummy(getParticleDef("MuMinus"));
    Vector3D position_init  = Vector3D(0, 0, 0);

    for(const std::pair<double, double> & offset: pairs) {
        for (const std::pair<double, double> &scatter: pairs) {
            for(Vector3D &direction_init: direction_list) {
                direction_init.CalculateSphericalCoordinates();
                dummy3->SetOffset(offset);
                dummy3->SetScatteringAngle(scatter);
                auto directions = dummy3->Scatter(1, 10., 1., position_init, direction_init, {0, 0, 0, 0,});
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
    std::shared_ptr<const Medium> medium = CreateMedium("StandardRock");
    Vector3D position_init  = Vector3D(0, 0, 0);
    Vector3D direction_init = Vector3D(1, 0, 0);
    direction_init.CalculateSphericalCoordinates();

    std::array<Scattering*, 2> scatter_list = {new ScatteringMoliere(getParticleDef("MuMinus"), medium),
                                               new ScatteringHighland(getParticleDef("MuMinus"), medium)};

    for(Scattering* scatter : scatter_list) {
        // Expect no change of direction for displacement of almost zero
        EXPECT_NEAR((std::get<0>(scatter->Scatter(1e-20, 1e4, 1e3, position_init, direction_init, {0.1, 0.2, 0.3, 0.4})) - direction_init).magnitude(),
                    0, 1e-10);
        EXPECT_NEAR((std::get<1>(scatter->Scatter(1e-20, 1e4, 1e3, position_init, direction_init, {0.1, 0.2, 0.3, 0.4})) - direction_init).magnitude(),
                    0, 1e-10);

        // Expect change of direction
        EXPECT_FALSE(std::get<0>(scatter->Scatter(1000, 1e4, 1e3, position_init, direction_init, {0.1, 0.2, 0.3, 0.4})) ==
                direction_init);
        EXPECT_FALSE(std::get<1>(scatter->Scatter(1000, 1e4, 1e3, position_init, direction_init, {0.1, 0.2, 0.3, 0.4})) ==
                direction_init);

    }
}

TEST(Scattering, FirstMomentum){
    RandomGenerator::Get().SetSeed(24601);
    std::shared_ptr<const Medium> medium = CreateMedium("StandardRock");
    Vector3D position_init  = Vector3D(0, 0, 0);
    Vector3D direction_init = Vector3D(0, 0, 1);
    direction_init.CalculateSphericalCoordinates();
    int statistics = 1e7;

    auto cross_dummy = std::make_shared<CrossSectionBuilder>("scattering_dummy");
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-4 * energy;});

    std::array<Scattering*, 3> scatter_list = {new ScatteringMoliere(getParticleDef("MuMinus"), medium),
                                               new ScatteringHighland(getParticleDef("MuMinus"), medium),
                                               new ScatteringHighlandIntegral<UtilityIntegral>(getParticleDef("MuMinus"), medium, CrossSectionList{cross_dummy})
                                                 };
    Vector3D scatter_sum;
    Vector3D offset_sum;

    for(Scattering* scatter: scatter_list){
        Vector3D scatter_sum = Vector3D(0., 0., 0.);
        Vector3D offset_sum = Vector3D(0., 0., 0.);

        for (int n=1; n<=statistics; ++n) {
            auto sampled_vectors = scatter->Scatter(
                    1e3, 1e4, 1e3, position_init, direction_init,
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



/*
TEST(Scattering, Scatter)
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
    bool first_line = true;

    while (in.good())
    {
        if (first_line)
        {
            in >> particleName >> mediumName >> parametrization >> ecut >> vcut >> energy_init >> energy_final >>
                distance >> rnd1 >> rnd2 >> rnd3 >> rnd4 >> x_f >> y_f >> z_f >> radius_f >> phi_f >> theta_f;

            first_line = false;
        }

        energy_previous = -1;

        ParticleDef particle_def = getParticleDef(particleName);

        std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);

        //reprouce old behaviour
        if(ecut==-1){
            ecut = std::numeric_limits<double>::infinity();
        }
        if(vcut==-1){
            vcut = 1;
        }

        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, false);

        Scattering* scattering = NULL;
        std::cout << parametrization << std::endl;
        if (parametrization == "HighlandIntegral")
        {
            auto brems = BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def, medium, ecuts, BremsstrahlungFactory::Definition(), std::make_shared<InterpolationDef>());
            auto ioniz = IonizationFactory::Get().CreateIonization(particle_def, medium, ecuts, IonizationFactory::Definition(), std::make_shared<InterpolationDef>());
            auto epair = EpairProductionFactory::Get().CreateEpairProduction(particle_def, medium, ecuts, EpairProductionFactory::Definition(), std::make_shared<InterpolationDef>());
            auto photo = PhotonuclearFactory::Get().CreatePhotonuclear(particle_def, medium, ecuts, PhotonuclearFactory::Definition(), std::make_shared<InterpolationDef>());
            auto cross = new CrossSectionList{std::shared_ptr<CrossSection>(brems), std::shared_ptr<CrossSection>(ioniz), std::shared_ptr<CrossSection>(epair), std::shared_ptr<CrossSection>(photo)};
            scattering = ScatteringFactory::Get().CreateScattering(parametrization, particle_def, medium, std::make_shared<InterpolationDef>(), std::unique_ptr<CrossSectionList>(cross));
        }
        else if(parametrization != "NoScattering")
        {
            scattering = ScatteringFactory::Get().CreateScattering(parametrization, particle_def, medium);

        }


        while (energy_previous < energy_init)
        {
            energy_previous = energy_init;
            if(parametrization!="NoScattering") {
                std::array<double, 4> rnd{rnd1, rnd2, rnd3, rnd4};
                auto directions = scattering->Scatter(distance,
                                                      energy_init,
                                                      energy_final,
                                                      position_init,
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
        }

        delete scattering;
    }
}
*/

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
