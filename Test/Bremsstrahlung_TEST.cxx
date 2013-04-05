#include "gtest/gtest.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>

using namespace std;



TEST(Bremsstrahlung , Test_of_dNdx_Interpolant ) {

    ifstream in;
    in.open("bin/Brems_dNdx_interpol.txt");

    ifstream inrnd;
    inrnd.open("bin/rnd.txt");
    double rnd;

    char firstLine[256];
    in.getline(firstLine,256);

    double dNdx;
    double dNdx_new;
    double energy;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);
    double energy_old=-1;
    CrossSections *brems = new Bremsstrahlung();
    bool first = true;
    while(in.good())
    {
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        first=false;
        energy_old = -1;
        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        //inrnd>>rnd;
        //cout << rnd << endl;
        //if(rnd>0.01)continue;
       CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);
        brems->EnableDNdxInterpolation();

        cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << med << "\t" << particleName<< "\t" << dNdx << endl;

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            dNdx_new=brems->CalculatedNdx();

            if(dNdx!=0)cout << fabs(dNdx_new -dNdx)/dNdx << endl;
            ASSERT_NEAR(dNdx_new, dNdx, 1*dNdx);

            in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        }



        delete cuts;
        delete medium;
        delete particle;
        delete brems;
    }
    inrnd.close();
}

TEST(Bremsstrahlung , Test_of_dNdx ) {

    ifstream in;
    in.open("bin/Brems_dNdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdx;
    double dNdx_new;
    double energy;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);


    while(in.good())
    {
        in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);


        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);

        dNdx_new=brems->CalculatedNdx();
        ASSERT_NEAR(dNdx_new, dNdx, 1e-6*dNdx);

        delete cuts;
        delete medium;
        delete particle;
        delete brems;



    }
}

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

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);


        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);



        dEdx_new=brems->CalculatedEdx();

        if(dEdx!=0)cout << fabs(dEdx_new -dEdx)/dEdx << endl;
        ASSERT_NEAR(dEdx_new, dEdx, 1e-8*dEdx);

        delete cuts;
        delete medium;
        delete particle;
        delete brems;



    }
}

TEST(Bremsstrahlung , Test_of_dEdx_Interpolant ) {

    ifstream in;
    //in.open("bin/Brems_dEdx_interpol.txt");
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
    double precision;
    double precisionOld = 1E-2;
    double precisionUran = 5E-2;
    while(in.good())
    {
        in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);

        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);

        brems->EnableDEdxInterpolation();
        dEdx_new=brems->CalculatedEdx();
        if(dEdx!=0){
            if(log10(fabs(dEdx_new -dEdx)/dEdx) > -3){
                cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << med << "\t" << particleName<< "\t" << dEdx << endl;
                cout << fabs(dEdx_new -dEdx)/dEdx << endl;
            }
        }
        if(!med.compare("uranium")){
            precision = precisionUran;
        }
        else{
            precision = precisionOld;
        }

        precision = 1;
        ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

        delete cuts;
        delete medium;
        delete particle;
        delete brems;



    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
