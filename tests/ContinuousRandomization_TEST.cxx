
// #include <iostream>
// #include <string>
// #include <vector>

#include "gtest/gtest.h"

#include "PROPOSAL/PROPOSAL.h"

using namespace std;
using namespace PROPOSAL;

class RndFromFile
{
private:
    double rnd_;
    string Path_;
    ifstream in_;

public:
    RndFromFile(string Path)
    {
        Path_ = Path;
        in_.open(Path_.c_str());
        in_ >> rnd_;
        if (!in_.good())
            log_warn("less than one rnd_number!");
    }

    double rnd()
    {
        in_ >> rnd_;
        if (!in_.good())
        {
            in_.close();
            in_.clear();
            in_.open(Path_.c_str());
            in_ >> rnd_;
        }
        return rnd_;
    }
};

ParticleDef getParticleDef(const string& name)
{
    if (name == "MuMinus")
    {
        return MuMinusDef::Get();
    }
    else if (name == "TauMinus")
    {
        return TauMinusDef::Get();
    }
    else
    {
        return EMinusDef::Get();
    }
}

class Test_Utilities : public ::testing::Test
{
protected:
    Test_Utilities()
    {
        std::cout << std::endl;
        std::cout << "new created" << std::endl;
        std::cout << std::endl;
    }

    virtual ~Test_Utilities() {}
    // virtual void TearDown() {}

    static Utility a;
    static Utility b;
};

Utility Test_Utilities::a(MuMinusDef::Get(), Ice(), EnergyCutSettings(), Utility::Definition(), InterpolationDef());
Utility Test_Utilities::b(TauMinusDef::Get(), Ice(), EnergyCutSettings(), Utility::Definition(), InterpolationDef());

// Utility utility_a(MuMinusDef::Get(), Ice(), EnergyCutSettings(), Utility::Definition(), InterpolationDef());
// Utility utility_b(TauMinusDef::Get(), Ice(), EnergyCutSettings(), Utility::Definition(), InterpolationDef());

TEST_F(Test_Utilities, Comparison_equal)
{
    ContinuousRandomizer A(a);
    ContinuousRandomizer B(a);

    EXPECT_TRUE(A == B);

    ContinuousRandomizer* C = new ContinuousRandomizer(a);
    ContinuousRandomizer* D = new ContinuousRandomizer(a);

    EXPECT_TRUE(*C == *D);

    delete C;
    delete D;

    A.Randomize(0, 1, 0.5);
    B.Randomize(0, 1, 0.5);

    EXPECT_TRUE(A == B);
}

TEST_F(Test_Utilities, Comparison_not_equal)
{
    ContinuousRandomizer A(a);
    ContinuousRandomizer B(b);

    // Only check different particle. If test does not fail, it should be good, because interanlly the utilitys are
    // compared. So the Utility test must work fine.
    EXPECT_TRUE(A != B);
}

TEST_F(Test_Utilities, Copyconstructor)
{
    ContinuousRandomizer A(a);
    ContinuousRandomizer B = A;

    EXPECT_TRUE(A == B);
}

TEST_F(Test_Utilities, Copyconstructor2)
{
    ContinuousRandomizer A(a);
    ContinuousRandomizer B(A);

    EXPECT_TRUE(A == B);

    Utility utility_c(a);

    EXPECT_TRUE(a == utility_c);

    ContinuousRandomizer C(a);
    ContinuousRandomizer D(utility_c, C);

    EXPECT_TRUE(A == B);
}

// TEST(ContinuousRandomization , Randomize ) {

//    ifstream in;
//    in.open("bin/TestFiles/ContinuousRandomization.txt");

//    char firstLine[256];
//    in.getline(firstLine,256);
//    double initial_energy;
//    double final_energy;
//    double randomized_energy;
//    double randomized_energy_new;
//    double vcut;
//    double ecut;
//    string mediumName;
//    string particleName;
//    double rnd;
//    cout.precision(16);

//    double max_diff = 0;
//    double diff = 0;
//    double sum_diff = 0;
//    int counter = 0;

//    while(in.good())
//    {
//        in>>rnd>>particleName>>mediumName>>ecut>>vcut>>initial_energy>>final_energy>>randomized_energy;
//        counter++;
//        Medium *medium = new Medium(Medium::GetTypeFromName(mediumName),1.);
//        Particle *particle = new Particle(PROPOSALParticle::GetTypeFromName(particleName),1.,1.,1,.20,20,1e5,10);
//        EnergyCutSettings *cut_settings = new EnergyCutSettings(ecut,vcut);

//        vector<CrossSections*> crosssections;

//        crosssections.resize(4);
//        crosssections.at(0) = new Ionization(particle, medium, cut_settings);
//        crosssections.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
//        crosssections.at(2) = new Photonuclear(particle, medium, cut_settings);
//        crosssections.at(3) = new Epairproduction(particle, medium, cut_settings);

//        ContinuousRandomization * cont = new ContinuousRandomization(particle,medium,crosssections);
//        randomized_energy_new = cont->Randomize(initial_energy,final_energy,rnd);

//        diff= abs(1 - randomized_energy_new/randomized_energy);
//        sum_diff += diff;

//        ASSERT_NEAR(randomized_energy_new, randomized_energy, 1e-1*randomized_energy);
//        cout<<randomized_energy_new<<"\t"<<randomized_energy<<"\t"<<diff<<endl;
//        if(diff>max_diff)
//        {
//            max_diff = diff;
//        }
//        for(unsigned int i = 0 ; i < crosssections.size() ; i++){

//            delete crosssections.at(i);
//        }

//        delete cut_settings;
//        delete medium;
//        delete particle;
//        delete cont;
//    }
//    cout<<"Maxi_diff "<<max_diff<<endl;
//    cout<<"Average diff "<<sum_diff/counter<<endl;

//}

TEST(ContinuousRandomization, Randomize_interpol)
{
    ifstream in;
    string filename = "bin/TestFiles/continous_randomization.txt";
    in.open(filename.c_str());

    if (!in.good())
    {
        std::cerr << "File \"" << filename << "\" not found" << std::endl;
    }

    char firstLine[256];
    in.getline(firstLine, 256);
    double initial_energy;
    double final_energy;
    double randomized_energy;
    double randomized_energy_new;
    double vcut;
    double ecut;
    string mediumName;
    string particleName;
    double rnd;
    double energy_old;
    bool first = true;

    cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        if (first)
            in >> rnd >> particleName >> mediumName >> ecut >> vcut >> initial_energy >> final_energy >>
                randomized_energy;

        first      = false;
        energy_old = -1;

        Medium* medium = MediumFactory::Get().CreateMedium(mediumName);
        EnergyCutSettings cut_settings(ecut, vcut);
        ParticleDef particle_def = getParticleDef(particleName);

        Utility utility(particle_def, *medium, cut_settings, Utility::Definition(), InterpolationDef());
        ContinuousRandomizer cont(utility, InterpolationDef());

        while (energy_old < initial_energy)
        {
            energy_old = initial_energy;

            randomized_energy_new = cont.Randomize(initial_energy, final_energy, rnd);

            ASSERT_NEAR(randomized_energy_new, randomized_energy, 1e-1 * randomized_energy);

            in >> rnd >> particleName >> mediumName >> ecut >> vcut >> initial_energy >> final_energy >>
                randomized_energy;
        }

        delete medium;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
