#include "gtest/gtest.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>

using namespace std;

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

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
