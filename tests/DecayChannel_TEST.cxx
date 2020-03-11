
#include "gtest/gtest.h"
#include <fstream>
#include <PROPOSAL/particle/Particle.h>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/decay/ManyBodyPhaseSpace.h"
#include "PROPOSAL/math/RandomGenerator.h"


using namespace PROPOSAL;

double matrix_element_evaluate(const DynamicData& p_condition, const std::vector<DynamicData>& products)
{
    double G_F = 1.1663787*1e-2; // MeV

    DynamicData electron = products.at(0);
    DynamicData numu = products.at(1);
    DynamicData nuebar = products.at(2);

    double p1 = p_condition.GetEnergy() * nuebar.GetEnergy() - (p_condition.GetMomentum() * p_condition.GetDirection()) * (nuebar.GetMomentum() * nuebar.GetDirection());
    double p2 = electron.GetEnergy() * numu.GetEnergy() - (electron.GetMomentum() * electron.GetDirection()) * (numu.GetMomentum() * numu.GetDirection());

    return 64 * G_F*G_F * p1 * p2;
}

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus")
    {
        return MuMinusDef::Get();
    } else if (name == "TauMinus")
    {
        return TauMinusDef::Get();
    } else if (name == "EMinus")
    {
        return EMinusDef::Get();
    } else if (name == "NuMu")
    {
        return NuMuBarDef::Get();
    } else if (name == "NuEBar"){
        return NuEBarDef::Get();
    } else if (name == "NuTau"){
        return NuTauDef::Get();
    } else{
        std::cout << "Unknown return particle";
        return MuMinusDef::Get();
    }
}

const std::string testfile_dir = "bin/TestFiles/";

ParticleDef mu  = MuMinusDef::Get();
ParticleDef tau = TauMinusDef::Get();

TEST(Comparison, Comparison_equal)
{
    LeptonicDecayChannel A(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
    LeptonicDecayChannel B(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
    EXPECT_TRUE(A == B);

    TwoBodyPhaseSpace C(mu, tau);
    TwoBodyPhaseSpace D(mu, tau);
    EXPECT_TRUE(C == D);

    DecayChannel* E = new LeptonicDecayChannel(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
    DecayChannel* F = new LeptonicDecayChannel(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());

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
    LeptonicDecayChannel A(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
    TwoBodyPhaseSpace B(mu, tau);
    EXPECT_TRUE(A != B);

    DecayChannel* E = new LeptonicDecayChannel(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
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
    LeptonicDecayChannel A(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
    DecayChannel* B = A.clone();

    EXPECT_TRUE(A == *B);

    delete B;

    DecayChannel* C = new LeptonicDecayChannel(EMinusDef::Get(), NuEDef::Get(), NuEBarDef::Get());
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

TEST(DecaySpectrum, MuMinus_Rest){
    std::string filename = testfile_dir + "Decay_MuMinus_rest.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    int statistic = 1e6;
    int NUM_bins = 50;
    std::string particleName;
    double init_energy;
    std::string particleDecay0;
    std::string particleDecay1;
    std::string particleDecay2;


    in >> statistic >> NUM_bins >> particleName >> init_energy >> particleDecay0 >> particleDecay1 >> particleDecay2;

    RandomGenerator::Get().SetSeed(1234);

    // LeptonicDecayChannelApproxTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    ParticleDef init_particle_def = getParticleDef(particleName);
    DynamicData init_particle(init_particle_def.particle_type);

    ParticleDef p0  = getParticleDef(particleDecay0);
    ParticleDef p1  = getParticleDef(particleDecay1);
    ParticleDef p2  = getParticleDef(particleDecay2);

    LeptonicDecayChannelApprox lep_approx(p0, p1, p2);
    std::vector<int> prod_0(NUM_bins, 0);
    std::vector<int> prod_1(NUM_bins, 0);
    std::vector<int> prod_2(NUM_bins, 0);

    double aux_energy, energy_sum;
    ParticleDef aux_def;

    init_particle.SetEnergy(init_energy);

    double v_max = ( std::pow(init_particle_def.mass,2 ) + pow(p0.mass,2))/(2 * init_particle_def.mass);
    double gamma = init_particle.GetEnergy() / init_particle_def.mass;
    double betagamma = init_particle.GetMomentum() / init_particle_def.mass;
    double max_energy = gamma * v_max + betagamma * std::sqrt( std::pow(v_max, 2) - std::pow(p0.mass,2));

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = lep_approx.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }

    int aux_bin;

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    // LeptonicDecayChannelTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    LeptonicDecayChannel lep(p0, p1, p2);

    std::fill(prod_0.begin(), prod_0.end(), 0); //reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = lep.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }


    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    // ManyBody

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<const ParticleDef*> daughters{&p0, &p1, &p2};

    ManyBodyPhaseSpace many_body(daughters, matrix_element_evaluate);

    std::fill(prod_0.begin(), prod_0.end(), 0); //reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = many_body.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }


    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    in.close();
}

TEST(DecaySpectrum, MuMinus_Energy){
    std::string filename = testfile_dir + "Decay_MuMinus_energy.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    int statistic = 1e6;
    int NUM_bins = 50;
    std::string particleName;
    double init_energy;
    std::string particleDecay0;
    std::string particleDecay1;
    std::string particleDecay2;


    in >> statistic >> NUM_bins >> particleName >> init_energy >> particleDecay0 >> particleDecay1 >> particleDecay2;

    RandomGenerator::Get().SetSeed(1234);

    // LeptonicDecayChannelApproxTest


    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    ParticleDef init_particle_def = getParticleDef(particleName);
    DynamicData init_particle(init_particle_def.particle_type);

    ParticleDef p0  = getParticleDef(particleDecay0);
    ParticleDef p1  = getParticleDef(particleDecay1);
    ParticleDef p2  = getParticleDef(particleDecay2);

    LeptonicDecayChannelApprox lep_approx(p0, p1, p2);
    std::vector<int> prod_0(NUM_bins, 0);
    std::vector<int> prod_1(NUM_bins, 0);
    std::vector<int> prod_2(NUM_bins, 0);

    double aux_energy, energy_sum;
    ParticleDef aux_def;

    init_particle.SetEnergy(init_energy);

    double v_max = ( std::pow(init_particle_def.mass,2 ) + pow(p0.mass,2))/(2 * init_particle_def.mass);
    double gamma = init_particle.GetEnergy() / init_particle_def.mass;
    double betagamma = init_particle.GetMomentum() / init_particle_def.mass;
    double max_energy = gamma * v_max + betagamma * std::sqrt( std::pow(v_max, 2) - std::pow(p0.mass,2));

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = lep_approx.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }

    int aux_bin;

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    // LeptonicDecayChannelTest


    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    LeptonicDecayChannel lep(p0, p1, p2);

    std::fill(prod_0.begin(), prod_0.end(), 0); //reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = lep.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }


    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    // ManyBody

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<const ParticleDef*> daughters{&p0, &p1, &p2};

    ManyBodyPhaseSpace many_body(daughters, matrix_element_evaluate);

    std::fill(prod_0.begin(), prod_0.end(), 0); //reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = many_body.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }


    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    in.close();
}

TEST(DecaySpectrum, TauMinus_Rest){
    std::string filename = testfile_dir + "Decay_TauMinus_rest.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    int statistic = 1e6;
    int NUM_bins = 50;
    std::string particleName;
    double init_energy;
    std::string particleDecay0;
    std::string particleDecay1;
    std::string particleDecay2;


    in >> statistic >> NUM_bins >> particleName >> init_energy >> particleDecay0 >> particleDecay1 >> particleDecay2;

    RandomGenerator::Get().SetSeed(1234);

    // LeptonicDecayChannelApproxTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    ParticleDef init_particle_def = getParticleDef(particleName);
    DynamicData init_particle(init_particle_def.particle_type);

    ParticleDef p0  = getParticleDef(particleDecay0);
    ParticleDef p1  = getParticleDef(particleDecay1);
    ParticleDef p2  = getParticleDef(particleDecay2);

    LeptonicDecayChannelApprox lep_approx(p0, p1, p2);
    std::vector<int> prod_0(NUM_bins, 0);
    std::vector<int> prod_1(NUM_bins, 0);
    std::vector<int> prod_2(NUM_bins, 0);

    double aux_energy, energy_sum;
    ParticleDef aux_def;

    init_particle.SetEnergy(init_energy);

    double v_max = ( std::pow(init_particle_def.mass,2 ) + pow(p0.mass,2))/(2 * init_particle_def.mass);
    double gamma = init_particle.GetEnergy() / init_particle_def.mass;
    double betagamma = init_particle.GetMomentum() / init_particle_def.mass;
    double max_energy = gamma * v_max + betagamma * std::sqrt( std::pow(v_max, 2) - std::pow(p0.mass,2));

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = lep_approx.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }

    int aux_bin;

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    // LeptonicDecayChannelTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    LeptonicDecayChannel lep(p0, p1, p2);

    std::fill(prod_0.begin(), prod_0.end(), 0); //reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = lep.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }


    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    // ManyBody

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<const ParticleDef*> daughters{&p0, &p1, &p2};

    ManyBodyPhaseSpace many_body(daughters, matrix_element_evaluate);

    std::fill(prod_0.begin(), prod_0.end(), 0); //reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = many_body.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }


    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    in.close();
}


TEST(DecaySpectrum, TauMinus_energy){
    std::string filename = testfile_dir + "Decay_TauMinus_energy.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    int statistic = 1e6;
    int NUM_bins = 50;
    std::string particleName;
    double init_energy;
    std::string particleDecay0;
    std::string particleDecay1;
    std::string particleDecay2;


    in >> statistic >> NUM_bins >> particleName >> init_energy >> particleDecay0 >> particleDecay1 >> particleDecay2;

    RandomGenerator::Get().SetSeed(1234);

    // LeptonicDecayChannelApproxTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    ParticleDef init_particle_def = getParticleDef(particleName);
    DynamicData init_particle(init_particle_def.particle_type);

    ParticleDef p0  = getParticleDef(particleDecay0);
    ParticleDef p1  = getParticleDef(particleDecay1);
    ParticleDef p2  = getParticleDef(particleDecay2);

    LeptonicDecayChannelApprox lep_approx(p0, p1, p2);
    std::vector<int> prod_0(NUM_bins, 0);
    std::vector<int> prod_1(NUM_bins, 0);
    std::vector<int> prod_2(NUM_bins, 0);

    double aux_energy, energy_sum;
    ParticleDef aux_def;

    init_particle.SetEnergy(init_energy);

    double v_max = ( std::pow(init_particle_def.mass,2 ) + pow(p0.mass,2))/(2 * init_particle_def.mass);
    double gamma = init_particle.GetEnergy() / init_particle_def.mass;
    double betagamma = init_particle.GetMomentum() / init_particle_def.mass;
    double max_energy = gamma * v_max + betagamma * std::sqrt( std::pow(v_max, 2) - std::pow(p0.mass,2));

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = lep_approx.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }

    int aux_bin;

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    // LeptonicDecayChannelTest

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    LeptonicDecayChannel lep(p0, p1, p2);

    std::fill(prod_0.begin(), prod_0.end(), 0); //reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = lep.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }


    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    // ManyBody

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //skip comment
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<const ParticleDef*> daughters{&p0, &p1, &p2};

    ManyBodyPhaseSpace many_body(daughters, matrix_element_evaluate);

    std::fill(prod_0.begin(), prod_0.end(), 0); //reset histogram
    std::fill(prod_1.begin(), prod_1.end(), 0);
    std::fill(prod_2.begin(), prod_2.end(), 0);

    for(int i=0; i<statistic; i++){
        init_particle.SetDirection(Vector3D(0, 0, -1));
        init_particle.SetPosition(Vector3D(0, 0, -1));
        init_particle.SetEnergy(init_energy);
        init_particle.SetPropagatedDistance(0);
        Secondaries aux = many_body.Decay(init_particle_def, init_particle);
        energy_sum = 0;
        for(DynamicData particle : aux.GetSecondaries()){
            aux_energy = particle.GetEnergy();
            if(particle.GetType() == p0.particle_type) {
                prod_0.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p1.particle_type) {
                prod_1.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            }
            else if(particle.GetType() == p2.particle_type) {
                prod_2.at((unsigned long) floor(aux_energy / max_energy * NUM_bins)) += 1;
            } else {
                FAIL() << "Unknown return particle";
            }
            energy_sum += aux_energy;
        }
        ASSERT_NEAR(energy_sum, init_energy, 1e-7*energy_sum); //conservation of energy
    }


    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_0.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_1.at(i), 1e-2*aux_bin);
    }

    for(int i=0; i<NUM_bins; i++){
        in >> aux_bin;
        ASSERT_NEAR(aux_bin, prod_2.at(i), 1e-2*aux_bin);
    }

    in.close();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
