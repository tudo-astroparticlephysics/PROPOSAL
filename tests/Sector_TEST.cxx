
#include "gtest/gtest.h"
#include <algorithm>

#include "PROPOSAL/PROPOSAL.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"

#include <string>
using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus") {
        return MuMinusDef();
    } else if (name == "TauMinus") {
        return TauMinusDef();
    } else {
        return EMinusDef();
    }
}

auto GetCrossSections(const ParticleDef& p_def, const Medium& med, std::shared_ptr<const EnergyCutSettings> cuts, bool interpolate) {
    // old TestFiles were created using the old StandardCrossSections that are recreated here
    crosssection_list_t cross;

    auto brems = crosssection::BremsKelnerKokoulinPetrukhin{ true, p_def, med };
    cross.push_back(make_crosssection(brems, p_def, med, cuts, interpolate));
    auto photo = crosssection::PhotoAbramowiczLevinLevyMaor97 { make_unique<crosssection::ShadowButkevichMikheyev>() };
    cross.push_back(make_crosssection(photo, p_def, med, cuts, interpolate));
    auto epair = crosssection::EpairKelnerKokoulinPetrukhin{ true, p_def, med };
    cross.push_back(make_crosssection(epair, p_def, med, cuts, interpolate));
    auto ioniz = crosssection::IonizBetheBlochRossi{ EnergyCutSettings(*cuts)};
    cross.push_back(make_crosssection(ioniz, p_def, med, cuts, interpolate));


    return cross;
}

TEST(Sector, Continuous)
{
    auto in = getTestFiles("Sector_ContinousLoss.txt");

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy;

    std::shared_ptr<ParticleDef> particle = std::make_shared<MuMinusDef>();
    std::shared_ptr<const Medium> medium = CreateMedium("ice");
    auto cuts = std::make_shared<EnergyCutSettings>(500, 1, false);

    auto cross = GetCrossSections(*particle, *medium, cuts, true);
    auto interaction = make_interaction(cross, true);
    auto decay = make_decay(cross, *particle, true);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        //reproduce old behaviour
        if (ecut == -1)
            ecut = INF;
        if (vcut == -1)
            vcut = 1;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {

            cuts = std::make_shared<EnergyCutSettings>(ecut, vcut, false);
            medium = CreateMedium(mediumName);
            particle = std::make_shared<ParticleDef>(getParticleDef(particleName));

            auto cross = GetCrossSections(*particle, *medium, cuts, true);
            interaction = make_interaction(cross, true);
            decay = make_decay(cross, *particle, true);
        }

        while (ss >> energy) {
            double rndd = RandomGenerator::Get().RandomDouble();
            double decay_energy = decay->EnergyDecay(initial_energy, rndd, medium->GetMassDensity());
            double rndi = RandomGenerator::Get().RandomDouble();
            double inter_energy
                = interaction->EnergyInteraction(initial_energy, rndi);
            double energy_calc = std::max(decay_energy, inter_energy);
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-10 * energy_calc));
            initial_energy = energy;
        }
    }
}


TEST(Sector, Stochastic)
{
    auto in = getTestFiles("Sector_StochasticLoss.txt");

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double energy, initial_energy, rnd;
    int interaction_type;

    std::shared_ptr<ParticleDef> particle = std::make_shared<MuMinusDef>();
    std::shared_ptr<const Medium> medium = CreateMedium("ice");
    auto cuts = std::make_shared<EnergyCutSettings>(500, 1, false);

    auto cross = GetCrossSections(*particle, *medium, cuts, true);
    auto interaction = make_interaction(cross, true);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        //reproduce old behaviour
        if (ecut == -1)
            ecut = INF;
        if (vcut == -1)
            vcut = 1;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            cuts = std::make_shared<EnergyCutSettings>(ecut, vcut, false);
            medium = CreateMedium(mediumName);
            particle = std::make_shared<ParticleDef>(getParticleDef(particleName));

            cross = GetCrossSections(*particle, *medium, cuts, true);
            interaction = make_interaction(cross, true);
        }

        while (ss >> energy >> interaction_type >> rnd) {
            // This section in necessary since we need to reproduce
            // the old behaviour of Sector::MakeStochasticLoss

            auto rates = interaction->Rates(initial_energy);

            double rnd1 = RandomGenerator::Get().RandomDouble();
            double rnd2 = RandomGenerator::Get().RandomDouble();
            double rnd3 = RandomGenerator::Get().RandomDouble();

            InteractionType type;
            double v_loss;

            // 1: determine the interaction type
            double sum = 0.;
            for (auto c : cross) {
                sum += c->CalculatedNdx(initial_energy);
            }
            sum *= rnd1;
            double rates_sum = 0;
            for (auto c : cross) {
                double rate_c = c->CalculatedNdx(initial_energy);
                rates_sum += rate_c;
                if (rates_sum > sum) {
                    type = c->GetInteractionType(); // save interaction_type
                    // 2: c is the right cross section, find the component
                    double rates_comp_sum = 0;
                    if (type == InteractionType::Ioniz) {
                        v_loss = c->CalculateStochasticLoss(medium->GetHash(), initial_energy, rnd2 * rate_c);
                    } else {
                        for (auto comp : medium->GetComponents()) {
                            double rate_for_comp = c->CalculatedNdx(initial_energy, comp.GetHash());
                            rates_comp_sum += rate_for_comp;
                            if (rates_comp_sum > rate_c * rnd3) {
                                // 3: comp is the right component. calculate loss
                                v_loss = c->CalculateStochasticLoss(comp.GetHash(), initial_energy,
                                                                    rnd2 * rate_for_comp);
                                break;
                            }
                        }
                    }
                    break;
                }
            }

            double energy_calc = (1 - v_loss) * initial_energy;
            double random = RandomGenerator::Get().RandomDouble();
            ASSERT_NEAR(random, rnd, random*1e-5);
            ASSERT_EQ((int)type, interaction_type);
            ASSERT_NEAR(energy_calc, energy, std::abs(1e-5 * energy_calc));
            initial_energy = energy;
        }
    }
}

TEST(Sector, EnergyDisplacement)
{
    auto in = getTestFiles("Sector_Energy_Distance.txt");

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(1234);
    std::string particleName;
    std::string mediumName;
    double ecut, vcut;
    double displacement, energy, initial_energy;

    std::shared_ptr<ParticleDef> particle = std::make_shared<MuMinusDef>();
    std::shared_ptr<const Medium> medium = CreateMedium("ice");
    auto cuts = std::make_shared<EnergyCutSettings>(500, 0.05, false);

    auto cross = GetCrossSections(*particle, *medium, cuts, true);
    auto displacement_calc = make_displacement(cross, true);

    std::string str;
    while (getline(in, str)) {
        std::istringstream ss(str);
        ss >> particleName >> mediumName >> ecut >> vcut >> initial_energy;

        //reproduce old behaviour
        if (ecut == -1)
            ecut = INF;
        if (vcut == -1)
            vcut = 1;

        if (particle->name != particleName || cuts->GetEcut() != ecut
            || cuts->GetVcut() != vcut || medium->GetName() != mediumName) {
            cuts = std::make_shared<EnergyCutSettings>(ecut, vcut, false);
            medium = CreateMedium(mediumName);
            particle = std::make_shared<ParticleDef>(getParticleDef(particleName));

            cross = GetCrossSections(*particle, *medium, cuts, true);
            displacement_calc = make_displacement(cross, true);
        }

        while (ss >> displacement >> energy) {
            double energy_calc;
            try {
                energy_calc = displacement_calc->UpperLimitTrackIntegral(initial_energy, displacement * medium->GetMassDensity());
            } catch (std::logic_error& e) {
                // exception thrown if distance can not be reached
                energy_calc = displacement_calc->GetLowerLim();
            }
            if (initial_energy * vcut == ecut) // old PROPOSAL version has inaccuracies for the interpolant here
                EXPECT_NEAR(energy_calc, energy, std::abs(1e-1 * energy_calc));
            else
                EXPECT_NEAR(energy_calc, energy, std::abs(1e-3 * energy_calc));
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
