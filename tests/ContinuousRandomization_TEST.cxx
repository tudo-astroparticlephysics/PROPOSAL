
#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/crossection/CrossSectionBuilder.h"

#include "PROPOSAL/crossection/factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crossection/factories/IonizationFactory.h"
#include "PROPOSAL/crossection/factories/EpairProductionFactory.h"
#include "PROPOSAL/crossection/factories/PhotonuclearFactory.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name) {
    if (name == "MuMinus") {
        return MuMinusDef::Get();
    } else if (name == "TauMinus") {
        return TauMinusDef::Get();
    } else {
        return EMinusDef::Get();
    }
}


TEST(ContinuousRandomization, Constructor) {
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", getParticleDef("MuMinus"));
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross_dummy->SetdE2dx_function([](double energy)->double {return std::exp(-14 + 2 * std::log(energy));});


    ContRand* dummy1 = new ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy});
    ContRand* dummy2 = new ContRandBuilder<UtilityInterpolant>(CrossSectionList{cross_dummy});

    auto dummy3 = ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy});
    auto dummy4 = ContRandBuilder<UtilityInterpolant>(CrossSectionList{cross_dummy});

}

TEST(ContinuousRandomization, FirstMomentum){
    // The first momentum of the randomized energy is supposed to be equal to the final_energy
    // when our distribution is neither cut by our upper limit (initial_energy) or our lower limit (mass/low)
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", getParticleDef("MuMinus"));
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross_dummy->SetdE2dx_function([](double energy)->double {return std::exp(-14 + 2 * std::log(energy));});

    RandomGenerator::Get().SetSeed(24601);
    auto contrand = ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy});
    auto energies = std::array<std::pair<double, double>, 5>{{{5e4, 1e5}, {5e6, 1e7}, {5e8, 1e9}, {5e10, 1e11}, {5e12, 1e13}}};
    int statistics = 1e4;
    double average, randomized;
    for(auto E : energies){
        average = 0;
        for(int n=1; n<=statistics; n++){
            randomized = contrand.EnergyRandomize(E.second, E.first, RandomGenerator::Get().RandomDouble());
            average = average + (randomized - average) / n;
        }
        EXPECT_FALSE(average == E.first);
        EXPECT_NEAR(average, E.first, E.first * 1e-2);
    }
}

TEST(ContinuousRandomization, IdenticalEnergies){
    // Expect the randomized energy to be equal to the final energy when the final energy is equal to the
    // initial energy (e.g. the continuous loss is zero)
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", getParticleDef("MuMinus"));
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross_dummy->SetdE2dx_function([](double energy)->double {return std::exp(-14 + 2 * std::log(energy));});
    auto contrand = ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy});
    double energy = 1e5;
    double energy_randomized = contrand.EnergyRandomize(energy, energy, RandomGenerator::Get().RandomDouble());
    EXPECT_DOUBLE_EQ(energy_randomized, energy);
}

TEST(ContinuousRandomization, PhysicalProperties){
    // Expecting an increasing variance in randomized_energies when the energy
    // difference between final_energy and an fixed initial energy is increasing
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", getParticleDef("MuMinus"));
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross_dummy->SetdE2dx_function([](double energy)->double {return std::exp(-14 + 2 * std::log(energy));});
    auto contrand = ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy});

    double E_i = 1e12;
    std::array<double, 5> final_energies = {6e11, 5.5e11, 5e11, 4.5e11, 4e11};

    int statistics = 1e4;
    auto average = std::pair<double, double>{0., 0.};
    double old_variance = 0;
    for(auto E_f : final_energies){
        average = {0., 0.};
        for(unsigned int n=1; n<=statistics; n++){
            double sampled = contrand.EnergyRandomize(E_i, E_f, RandomGenerator::Get().RandomDouble());
            average = welfords_online_algorithm(sampled, n, average.first, average.second);
        }
        EXPECT_TRUE(old_variance < average.second);
        old_variance = average.second;
    }
}

TEST(ContinuousRandomization, Constraints){
    // The randomized energy should never be above the particle mass or above the initial_energy
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", getParticleDef("MuMinus"));
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross_dummy->SetdE2dx_function([](double energy)->double {return std::exp(-14 + 2 * std::log(energy));});
    auto contrand = ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy});
    double energy = 1e5;
    int statistics = 1e3;
    double randomized;

    for(int n=1; n<statistics; n++){
        randomized = contrand.EnergyRandomize(1e4, MMU, RandomGenerator::Get().RandomDouble());
        EXPECT_TRUE(randomized >= MMU);
        EXPECT_TRUE(randomized <= 1e4);

        randomized = contrand.EnergyRandomize(110, MMU, RandomGenerator::Get().RandomDouble());
        EXPECT_TRUE(randomized <= 110);
        EXPECT_TRUE(randomized > MMU);
    }

}

TEST(ContinuousRandomization, compare_integral_interpolant) {
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy", getParticleDef("MuMinus"));
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross_dummy->SetdE2dx_function([](double energy)->double {return std::exp(-14 + 2 * std::log(energy));});

    auto contrand_integral = ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy});
    auto contrand_interpol = ContRandBuilder<UtilityInterpolant>(CrossSectionList{cross_dummy});

    auto energies = std::array<double, 5> {1e6, 1e7, 1e8, 1e9, 1e10};

    for(auto E_i : energies){
        double rnd = RandomGenerator::Get().RandomDouble();
        double randomized_integral = contrand_integral.EnergyRandomize(E_i, 1e5, rnd);
        double randomized_interpol = contrand_interpol.EnergyRandomize(E_i, 1e5, rnd);
        EXPECT_NEAR(randomized_integral, randomized_interpol, 1e-5 * E_i);
    }
}

TEST(ContinuousRandomization, Randomize_interpol) {
    std::string filename = "bin/TestFiles/continous_randomization.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";


    char firstLine[256];
    in.getline(firstLine, 256);
    double initial_energy;
    double final_energy;
    double randomized_energy;
    double randomized_energy_new;
    double vcut;
    double ecut;
    std::string mediumName;
    std::string particleName;
    double rnd;
    double energy_old;
    bool first = true;

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in.good()) {
        if (first)
            in >> rnd >> particleName >> mediumName >> ecut >> vcut >>
                initial_energy >> final_energy >> randomized_energy;

        first = false;
        energy_old = -1;

        std::shared_ptr<const Medium> medium  = CreateMedium(mediumName);
        ParticleDef particle_def = getParticleDef(particleName);

        //reprouce old behaviour
        if(ecut==-1){
            ecut = std::numeric_limits<double>::infinity();
        }
        if(vcut==-1){
            vcut = 1;
        }

        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, false);

        auto brems = BremsstrahlungFactory::Get().CreateBremsstrahlung(particle_def, medium, ecuts, BremsstrahlungFactory::Definition(), std::make_shared<InterpolationDef>());
        auto ioniz = IonizationFactory::Get().CreateIonization(particle_def, medium, ecuts, IonizationFactory::Definition(), std::make_shared<InterpolationDef>());
        auto epair = EpairProductionFactory::Get().CreateEpairProduction(particle_def, medium, ecuts, EpairProductionFactory::Definition(), std::make_shared<InterpolationDef>());
        auto photo = PhotonuclearFactory::Get().CreatePhotonuclear(particle_def, medium, ecuts, PhotonuclearFactory::Definition(), std::make_shared<InterpolationDef>());
        auto cross = new CrossSectionList{std::shared_ptr<CrossSection>(brems), std::shared_ptr<CrossSection>(ioniz), std::shared_ptr<CrossSection>(epair), std::shared_ptr<CrossSection>(photo)};

        ContRandBuilder<UtilityIntegral> cont(*cross);

        while (energy_old < initial_energy) {
            energy_old = initial_energy;

            randomized_energy_new =
                cont.EnergyRandomize(initial_energy, final_energy, rnd);

            ASSERT_NEAR(randomized_energy_new, randomized_energy,
                        1e-1 * randomized_energy);

            in >> rnd >> particleName >> mediumName >> ecut >> vcut >>
                initial_energy >> final_energy >> randomized_energy;
        }
    }
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
