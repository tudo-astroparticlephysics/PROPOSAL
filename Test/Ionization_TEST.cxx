#include "gtest/gtest.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>

using namespace std;

class RndFromFile{
private:
    double rnd_;
    string Path_;
    ifstream in_;

public:
    RndFromFile(string Path){
        Path_ = Path;
        in_.open(Path_.c_str());
        in_>>rnd_;
        if(!in_.good())cout << "less than one rnd_number!" << endl;
    }

    double rnd(){
        in_>>rnd_;
        if(!in_.good()){
            in_.close();
            in_.clear();
            in_.open(Path_.c_str());
            in_>>rnd_;
        }
        return rnd_;
    }
};

TEST(Ionization , Test_of_dEdx ) {

    ifstream in;
    in.open("bin/Ioniz_dEdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dEdx_new;
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;

    cout.precision(16);


    while(in.good())
    {
        in>>ecut>>vcut>>energy>>med>>particleName>>dEdx;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *ioniz = new Ionization(particle, medium, cuts);

        dEdx_new=ioniz->CalculatedEdx();
        ASSERT_NEAR(dEdx_new, dEdx, 1e-12*dEdx);

        delete cuts;
        delete medium;
        delete particle;
        delete ioniz;
    }
}

TEST(Ionization , Test_of_dNdx ) {

    ifstream in;
    in.open("bin/Ioniz_dNdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdx_new;
    double energy;
    double dNdx;
    double ecut;
    double vcut;
    string med;
    string particleName;

    cout.precision(16);


    while(in.good())
    {
        in>>ecut>>vcut>>energy>>med>>particleName>>dNdx;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *ioniz = new Ionization(particle, medium, cuts);

        dNdx_new=ioniz->CalculatedNdx();
        ASSERT_NEAR(dNdx_new, dNdx, 1e-8*dNdx);

        delete cuts;
        delete medium;
        delete particle;
        delete ioniz;

    }
}

TEST(Ionization , Test_of_dNdxrnd ) {

    ifstream in;
    in.open("bin/Ioniz_dNdxrnd.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double dNdxrnd;
    double dNdxrnd_new;
    double energy;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);
    double energy_old=-1;

    RndFromFile* Rand = new RndFromFile("bin/rnd.txt");

    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>energy>>med>>particleName>>dNdxrnd;
        first = false;
        energy_old = -1;
        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *ioniz = new Ionization(particle, medium, cuts);

        while(energy_old < energy){
            energy_old = energy;
            ioniz->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=ioniz->CalculatedNdx(Rand->rnd());

            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, 1E-7*dNdxrnd);

            in>>ecut>>vcut>>energy>>med>>particleName>>dNdxrnd;
        }



        delete cuts;
        delete medium;
        delete particle;
        delete ioniz;
    }
    delete Rand;
}

TEST(Ionization , Test_of_e ) {

    ifstream in;
    in.open("bin/Ioniz_e.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double e;
    double energy;
    double e_new;
    double ecut;
    double vcut;
    string med;
    string particleName;

    cout.precision(16);

    double precision = 1E-8;

    RndFromFile* Rand = new RndFromFile("bin/rnd.txt");
    RndFromFile* Rand2 = new RndFromFile("bin/rnd.txt");
    Rand2->rnd();

    double rnd1,rnd2;

    double energy_old;
    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>energy>>med>>particleName>>e;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
        CrossSections *ioniz = new Ionization(particle, medium, cuts);

        while(energy_old < energy)
        {
            energy_old = energy;
            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            ioniz->GetParticle()->SetEnergy(energy);
            e_new=ioniz->CalculateStochasticLoss(rnd1,rnd2);

            ASSERT_NEAR(e_new, e, precision*e);

            in>>ecut>>vcut>>energy>>med>>particleName>>e;
        }
        delete cuts;
        delete medium;
        delete particle;
        delete ioniz;
    }
}

TEST(Ionization , Test_of_dEdx_Interpolant ) {

    ifstream in;
    in.open("bin/Ioniz_dEdx_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dEdx_new;
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;

    cout.precision(16);

    double precision = 1E-8;

    double energy_old;
    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>energy>>med>>particleName>>dEdx;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
        CrossSections *ioniz = new Ionization(particle, medium, cuts);

        ioniz->EnableDEdxInterpolation();

        while(energy_old < energy)
        {
            energy_old = energy;

            ioniz->GetParticle()->SetEnergy(energy);
            dEdx_new=ioniz->CalculatedEdx();

            ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

            in>>ecut>>vcut>>energy>>med>>particleName>>dEdx;
        }
        delete cuts;
        delete medium;
        delete particle;
        delete ioniz;
    }
}

TEST(Ionization , Test_of_dNdx_Interpolant ) {

    ifstream in;
    in.open("bin/Ioniz_dNdx_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdx_new;
    double energy;
    double dNdx;
    double ecut;
    double vcut;
    string med;
    string particleName;

    cout.precision(16);

    double precision = 1E-5;

    double energy_old;
    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>energy>>med>>particleName>>dNdx;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
        CrossSections *ioniz = new Ionization(particle, medium, cuts);

        ioniz->EnableDNdxInterpolation();

        while(energy_old < energy)
        {
            energy_old = energy;

            ioniz->GetParticle()->SetEnergy(energy);
            dNdx_new=ioniz->CalculatedNdx();

            ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

            in>>ecut>>vcut>>energy>>med>>particleName>>dNdx;
        }
        delete cuts;
        delete medium;
        delete particle;
        delete ioniz;
    }
}

TEST(Ionization , Test_of_dNdxrnd_interpol ) {

    ifstream in;
    in.open("bin/Ioniz_dNdxrnd_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double dNdxrnd;
    double dNdxrnd_new;
    double energy;
    double ecut;
    double vcut;
    string med;
    string particleName;

    cout.precision(16);
    double energy_old=-1;

    RndFromFile* Rand = new RndFromFile("bin/rnd.txt");

    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>energy>>med>>particleName>>dNdxrnd;
        first = false;
        energy_old = -1;
        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *ioniz = new Ionization(particle, medium, cuts);
        ioniz->EnableDNdxInterpolation();

        while(energy_old < energy){
            energy_old = energy;
            ioniz->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=ioniz->CalculatedNdx(Rand->rnd());

            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, 1E-5*dNdxrnd);

            in>>ecut>>vcut>>energy>>med>>particleName>>dNdxrnd;
        }



        delete cuts;
        delete medium;
        delete particle;
        delete ioniz;
    }
    delete Rand;
}


TEST(Ionization , Test_of_e_interpol ) {

    ifstream in;
    in.open("bin/Ioniz_e_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double e;
    double energy;
    double e_new;
    double ecut;
    double vcut;
    string med;
    string particleName;

    cout.precision(16);

    double precision = 1E-2;

    RndFromFile* Rand = new RndFromFile("bin/rnd.txt");
    RndFromFile* Rand2 = new RndFromFile("bin/rnd.txt");
    Rand2->rnd();

    double rnd1,rnd2;

    double energy_old;
    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>energy>>med>>particleName>>e;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
        CrossSections *ioniz = new Ionization(particle, medium, cuts);

        ioniz->EnableDEdxInterpolation();

        while(energy_old < energy)
        {
            energy_old = energy;
            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            ioniz->GetParticle()->SetEnergy(energy);

            e_new=ioniz->CalculateStochasticLoss(rnd1,rnd2);

            if(e!=0)if(log10(fabs(1-e_new/e))>-3){
                cout<< "\t" << ecut<< "\t" << vcut << "\t" << energy<< "\t" << med<< "\t" << particleName<<endl;
                cout << log10(fabs(1-e_new/e)) << endl;
            }
            //cout << "rnd1: " << rnd1 << "\t rnd2: " << rnd2 << endl;
            ASSERT_NEAR(e_new, e, precision*e);

            in>>ecut>>vcut>>energy>>med>>particleName>>e;
        }
        delete cuts;
        delete medium;
        delete particle;
        delete ioniz;
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
