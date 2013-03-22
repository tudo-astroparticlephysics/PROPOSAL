#include "gtest/gtest.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>

using namespace std;

TEST(Bremsstrahlung , Test_of_dEdx ) {

    ifstream in;
    in.open("bin/Brems_dEdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dEdx_new;
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);


    while(in.good())
    {
        in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        cout<<para<<"\t"<<ecut<<"\t"<<vcut<<"\t"<<lpm<<"\t"<<energy<<"\t"<<med<<"\t"<<particleName<<"\t"<<dEdx<<"\t";

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);

        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);
        dEdx_new=brems->CalculatedEdx();
        cout<<dEdx_new<<"\t"<<particle->GetMass()<<endl;
        //ASSERT_NEAR(dEdx_new, dEdx, 5e-1*dEdx);

        delete brems;
        delete cuts;
        delete medium;
        delete particle;




    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
