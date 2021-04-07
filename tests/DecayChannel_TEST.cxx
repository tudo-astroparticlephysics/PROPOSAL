
#include "gtest/gtest.h"
#include <PROPOSAL/particle/Particle.h>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/ManyBodyPhaseSpace.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSALTestUtilities/TestFilesHandling.h"

#include <memory>

using namespace PROPOSAL;

double matrix_element_evaluate(const ParticleState& p_condition,
    const std::vector<ParticleState>& products)
{
    double G_F = 1.1663787 * 1e-2; // MeV

    ParticleState electron = products[0];
    ParticleState numu = products[1];
    ParticleState nuebar = products[2];

    double p1 = p_condition.energy * nuebar.energy
        - (p_condition.GetMomentum() * p_condition.direction)
            * (nuebar.GetMomentum() * nuebar.direction);
    double p2 = electron.energy * numu.energy
        - (electron.GetMomentum() * electron.direction)
            * (numu.GetMomentum() * numu.direction);

    return 64 * G_F * G_F * p1 * p2;
}

inline auto getParticleDef(const std::string& name)
{
    auto ptr = std::make_unique<ParticleDef>();
    if (name == "MuMinus") {
        ptr = std::make_unique<MuMinusDef>();
    } else if (name == "TauMinus") {
        ptr = std::make_unique<TauMinusDef>();
    } else if (name == "EMinus") {
        ptr = std::make_unique<EMinusDef>();
    } else if (name == "NuMu") {
        ptr = std::make_unique<NuMuBarDef>();
    } else if (name == "NuEBar") {
        ptr = std::make_unique<NuEBarDef>();
    } else if (name == "NuTau") {
        ptr = std::make_unique<NuTauDef>();
    } else {
        std::cout << "Unknown return particle";
        ptr = std::make_unique<MuMinusDef>();
    }
    return ptr;
}

const std::string testfile_dir = "tests/TestFiles/";

ParticleDef mu = MuMinusDef();
ParticleDef tau = TauMinusDef();
ParticleDef eminus = EMinusDef();
ParticleDef nue = NuMuDef();
ParticleDef nuebar = NuEBarDef();

TEST(Comparison, Comparison_equal)
{
    LeptonicDecayChannel A(eminus, nue, nuebar);
    LeptonicDecayChannel B(eminus, nue, nuebar);
    EXPECT_TRUE(A == B);

    TwoBodyPhaseSpace C(mu, tau);
    TwoBodyPhaseSpace D(mu, tau);
    EXPECT_TRUE(C == D);

    DecayChannel* E = new LeptonicDecayChannel(eminus, nue, nuebar);
    DecayChannel* F = new LeptonicDecayChannel(eminus, nue, nuebar);

    EXPECT_TRUE(*E == *F);
    delete E;
    delete F;

    E = new TwoBodyPhaseSpace(mu, tau);
    F = new TwoBodyPhaseSpace(mu, tau);

    EXPECT_TRUE(*E == *F);
    delete E;
    delete F;
}

TEST(Comparison, Comparison_not_equal)
{
    LeptonicDecayChannel A(eminus, nue, nuebar);
    TwoBodyPhaseSpace B(mu, tau);
    EXPECT_TRUE(A != B);

    DecayChannel* E = new LeptonicDecayChannel(eminus, nue, nuebar);
    DecayChannel* F = new TwoBodyPhaseSpace(mu, tau);

    EXPECT_TRUE(*E != *F);
    delete E;
    delete F;

    E = new TwoBodyPhaseSpace(mu, mu);
    F = new TwoBodyPhaseSpace(mu, tau);

    EXPECT_TRUE(*E != *F);
    delete E;
    delete F;
}

TEST(Assignment, Copyconstructor)
{
    LeptonicDecayChannel A(eminus, nue, nuebar);
    DecayChannel* B = A.clone();

    EXPECT_TRUE(A == *B);

    delete B;

    DecayChannel* C = new LeptonicDecayChannel(eminus, nue, nuebar);
    DecayChannel* D = C->clone();

    EXPECT_TRUE(*C == *D);

    delete C;
    delete D;

    TwoBodyPhaseSpace E(mu, tau);
    DecayChannel* F = E.clone();

    EXPECT_TRUE(E == *F);

    delete F;

    DecayChannel* G = new TwoBodyPhaseSpace(mu, tau);
    DecayChannel* H = G->clone();

    EXPECT_TRUE(*G == *H);

    delete G;
    delete H;
}

TEST(DecaySpectrum, MuMinus_Rest)
{
    auto in = getTestFiles("Decay_MuMinus_rest.txt");

    int statistic = 1e6;
    int NUM_bins = 50;
    std::string particleName;
    double init_energy;
    std::string particleDecay0;
    std::string particleDecay1;
    std::string particleDecay2;

    in >> statistic >> NUM_bins >> particleName >> init_energy >> particleDecay0
        >> particleDecay1 >> particleDecay2;

    RandomGenerator::Get().SetSeed(1234);

    // LeptonicDecayChannelApproxTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    ParticleDef init_particle_def(*getParticleDef(particleName));
    ParticleState init_particle;
    init_particle.type = init_particle_def.particle_type;

    ParticleDef p0(*getParticleDef(particleDecay0));
    ParticleDef p1(*getParticleDef(particleDecay1));
    ParticleDef p2(*getParticleDef(particleDecay2));

    LeptonicDecayChannelApprox lep_approx(p0, p1, p2);
    std::vector<int> prod_0(NUM_bins, 0);
    std::vector<int> prod_1(NUM_bins, 0);
    std::vector<int> prod_2(NUM_bins, 0);

    double aux_energy, energy_sum;

    init_particle.energy = init_energy;

    double v_max = (std::pow(init_particle_def.mass, 2) + pow(p0.mass, 2))
        / (2 * init_particle_def.mass);
    double gamma = init_particle.energy / init_particle_def.mass;
    double betagamma = init_particle.GetMomentum() / init_particle_def.mass;
    double max_energy = gamma * v_max
        + betagamma * std::sqrt(std::pow(v_max, 2) - std::pow(p0.mass, 2));

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0.;
        auto aux = lep_approx.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    int aux_bin;

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    // LeptonicDecayChannelTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    LeptonicDecayChannel lep(p0, p1, p2);

    std::fill(prod_0.begin(), prod_0.end(), 0); // reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = lep.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    // ManyBody

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<std::shared_ptr<const ParticleDef>> daughters = {
        std::make_shared<ParticleDef>(p0),
        std::make_shared<ParticleDef>(p1),
        std::make_shared<ParticleDef>(p2),
    };

    ManyBodyPhaseSpace many_body(daughters, matrix_element_evaluate);

    std::fill(prod_0.begin(), prod_0.end(), 0); // reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = many_body.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    in.close();
}

TEST(DecaySpectrum, MuMinus_Energy)
{
    auto in = getTestFiles("Decay_MuMinus_energy.txt");

    int statistic = 1e6;
    int NUM_bins = 50;
    std::string particleName;
    double init_energy;
    std::string particleDecay0;
    std::string particleDecay1;
    std::string particleDecay2;

    in >> statistic >> NUM_bins >> particleName >> init_energy >> particleDecay0
        >> particleDecay1 >> particleDecay2;

    RandomGenerator::Get().SetSeed(1234);

    // LeptonicDecayChannelApproxTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    ParticleDef init_particle_def(*getParticleDef(particleName));
    ParticleState init_particle;
    init_particle.type = init_particle_def.particle_type;

    ParticleDef p0(*getParticleDef(particleDecay0));
    ParticleDef p1(*getParticleDef(particleDecay1));
    ParticleDef p2(*getParticleDef(particleDecay2));

    LeptonicDecayChannelApprox lep_approx(p0, p1, p2);
    std::vector<int> prod_0(NUM_bins, 0);
    std::vector<int> prod_1(NUM_bins, 0);
    std::vector<int> prod_2(NUM_bins, 0);

    double aux_energy, energy_sum;

    init_particle.energy = init_energy;

    double v_max = (std::pow(init_particle_def.mass, 2) + pow(p0.mass, 2))
        / (2 * init_particle_def.mass);
    double gamma = init_particle.energy / init_particle_def.mass;
    double betagamma = init_particle.GetMomentum() / init_particle_def.mass;
    double max_energy = gamma * v_max
        + betagamma * std::sqrt(std::pow(v_max, 2) - std::pow(p0.mass, 2));

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = lep_approx.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    int aux_bin;

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    // LeptonicDecayChannelTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    LeptonicDecayChannel lep(p0, p1, p2);

    std::fill(prod_0.begin(), prod_0.end(), 0); // reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = lep.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    // ManyBody

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<std::shared_ptr<const ParticleDef>> daughters = {
        std::make_shared<ParticleDef>(p0),
        std::make_shared<ParticleDef>(p1),
        std::make_shared<ParticleDef>(p2),
    };

    ManyBodyPhaseSpace many_body(daughters, matrix_element_evaluate);

    std::fill(prod_0.begin(), prod_0.end(), 0); // reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = many_body.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    in.close();
}

TEST(DecaySpectrum, TauMinus_Rest)
{
    auto in = getTestFiles("Decay_TauMinus_rest.txt");

    int statistic = 1e6;
    int NUM_bins = 50;
    std::string particleName;
    double init_energy;
    std::string particleDecay0;
    std::string particleDecay1;
    std::string particleDecay2;

    in >> statistic >> NUM_bins >> particleName >> init_energy >> particleDecay0
        >> particleDecay1 >> particleDecay2;

    RandomGenerator::Get().SetSeed(1234);

    // LeptonicDecayChannelApproxTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    ParticleDef init_particle_def(*getParticleDef(particleName));
    ParticleState init_particle;
    init_particle.type = init_particle_def.particle_type;

    ParticleDef p0(*getParticleDef(particleDecay0));
    ParticleDef p1(*getParticleDef(particleDecay1));
    ParticleDef p2(*getParticleDef(particleDecay2));

    LeptonicDecayChannelApprox lep_approx(p0, p1, p2);
    std::vector<int> prod_0(NUM_bins, 0);
    std::vector<int> prod_1(NUM_bins, 0);
    std::vector<int> prod_2(NUM_bins, 0);

    double aux_energy, energy_sum;

    init_particle.energy = init_energy;

    double v_max = (std::pow(init_particle_def.mass, 2) + pow(p0.mass, 2))
        / (2 * init_particle_def.mass);
    double gamma = init_particle.energy / init_particle_def.mass;
    double betagamma = init_particle.GetMomentum() / init_particle_def.mass;
    double max_energy = gamma * v_max
        + betagamma * std::sqrt(std::pow(v_max, 2) - std::pow(p0.mass, 2));

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = lep_approx.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    int aux_bin;

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    // LeptonicDecayChannelTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    LeptonicDecayChannel lep(p0, p1, p2);

    std::fill(prod_0.begin(), prod_0.end(), 0); // reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = lep.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    // ManyBody

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<std::shared_ptr<const ParticleDef>> daughters = {
        std::make_shared<ParticleDef>(p0),
        std::make_shared<ParticleDef>(p1),
        std::make_shared<ParticleDef>(p2),
    };

    ManyBodyPhaseSpace many_body(daughters, matrix_element_evaluate);

    std::fill(prod_0.begin(), prod_0.end(), 0); // reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = many_body.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    in.close();
}

TEST(DecaySpectrum, TauMinus_energy)
{
    auto in = getTestFiles("Decay_TauMinus_energy.txt");

    int statistic = 1e6;
    int NUM_bins = 50;
    std::string particleName;
    double init_energy;
    std::string particleDecay0;
    std::string particleDecay1;
    std::string particleDecay2;

    in >> statistic >> NUM_bins >> particleName >> init_energy >> particleDecay0
        >> particleDecay1 >> particleDecay2;

    RandomGenerator::Get().SetSeed(1234);

    // LeptonicDecayChannelApproxTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    ParticleDef init_particle_def(*getParticleDef(particleName));
    ParticleState init_particle;
    init_particle.type = init_particle_def.particle_type;

    ParticleDef p0(*getParticleDef(particleDecay0));
    ParticleDef p1(*getParticleDef(particleDecay1));
    ParticleDef p2(*getParticleDef(particleDecay2));

    LeptonicDecayChannelApprox lep_approx(p0, p1, p2);
    std::vector<int> prod_0(NUM_bins, 0);
    std::vector<int> prod_1(NUM_bins, 0);
    std::vector<int> prod_2(NUM_bins, 0);

    double aux_energy, energy_sum;

    init_particle.energy = init_energy;

    double v_max = (std::pow(init_particle_def.mass, 2) + pow(p0.mass, 2))
        / (2 * init_particle_def.mass);
    double gamma = init_particle.energy / init_particle_def.mass;
    double betagamma = init_particle.GetMomentum() / init_particle_def.mass;
    double max_energy = gamma * v_max
        + betagamma * std::sqrt(std::pow(v_max, 2) - std::pow(p0.mass, 2));

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = lep_approx.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    int aux_bin;

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    // LeptonicDecayChannelTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    LeptonicDecayChannel lep(p0, p1, p2);

    std::fill(prod_0.begin(), prod_0.end(), 0); // reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = lep.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    // ManyBody

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip
                                                                  // comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<std::shared_ptr<const ParticleDef>> daughters = {
        std::make_shared<ParticleDef>(p0),
        std::make_shared<ParticleDef>(p1),
        std::make_shared<ParticleDef>(p2),
    };

    ManyBodyPhaseSpace many_body(daughters, matrix_element_evaluate);

    std::fill(prod_0.begin(), prod_0.end(), 0); // reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for (int i = 0; i < statistic; i++) {
        init_particle.direction = Cartesian3D(0, 0, -1);
        init_particle.position = Cartesian3D(0, 0, -1);
        init_particle.energy = init_energy;
        init_particle.propagated_distance = 0;
        auto aux = many_body.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for (ParticleState particle : aux) {
            aux_energy = particle.energy;
            if (particle.type == p0.particle_type) {
                prod_0[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p1.particle_type) {
                prod_1[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else if (particle.type == p2.particle_type) {
                prod_2[floor(aux_energy / max_energy * NUM_bins)] += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy,
            1e-7 * energy_sum); // conservation of energy
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1[i], 1e-2 * aux_bin);
    }

    for (int i = 0; i < NUM_bins; i++) {
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2[i], 1e-2 * aux_bin);
    }

    in.close();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
