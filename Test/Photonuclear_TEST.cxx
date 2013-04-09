#include "gtest/gtest.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>

using namespace std;


TEST(Photonuclear , Test_of_dEdx ) {

    ifstream in;
    in.open("bin/Photo_dEdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    int shadow;
    int bb;
    double dEdx_new;
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    int para;

    cout.precision(16);


    while(in.good())
    {

        in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dEdx;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *photo = new Photonuclear(particle, medium, cuts);

        if(bb==1&&para<3)photo->SetParametrization(1);
        if(bb==2&&para<3)photo->SetParametrization(2);
        if(bb==3&&para<3)photo->SetParametrization(3);
        if(bb==4&&para<3)photo->SetParametrization(4);
        if(bb==1&&para==3)photo->SetParametrization(5);
        if(bb==2&&para==3)photo->SetParametrization(6);
        if(bb==1&&para==4)photo->SetParametrization(7);



        if(para==2)photo->SetHardComponent(true);
        else photo->SetHardComponent(false);



        photo->SetShadow(shadow);

        dEdx_new=photo->CalculatedEdx();
        ASSERT_NEAR(dEdx_new, dEdx, 1e-6*dEdx);

        delete cuts;
        delete medium;
        delete particle;
        delete photo;



    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
