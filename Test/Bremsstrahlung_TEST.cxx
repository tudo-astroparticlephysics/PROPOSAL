#include "gtest/gtest.h"
#include "PROPOSAL/Bremsstrahlung.h"
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
        if(in_.good()){
            in_.close();
            in_.clear();
            in_.open(Path_.c_str());
            in_>>rnd_;
        }
        return rnd_;
    }
};

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
        ASSERT_NEAR(dNdx_new, dNdx, 1e-14*dNdx);

        delete cuts;
        delete medium;
        delete particle;
        delete brems;



    }
}

TEST(Bremsstrahlung , Test_of_dNdxrnd ) {

    ifstream in;
    in.open("bin/Brems_dNdxrnd.txt");

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
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        first=false;
        energy_old = -1;
        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

       CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);

        //cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << med << "\t" << particleName<< "\t" << dNdx << endl;

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=brems->CalculatedNdx(Rand->rnd());

            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, 1E-6*dNdxrnd);

            in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        }



        delete cuts;
        delete medium;
        delete particle;
        delete brems;
    }
    delete Rand;
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

        ASSERT_NEAR(dEdx_new, dEdx, 1e-8*dEdx);

        delete cuts;
        delete medium;
        delete particle;
        delete brems;



    }
}

TEST(Bremsstrahlung , Test_of_dEdx_Interpolant ) {

    ifstream in;
    in.open("bin/Brems_dEdx_interpol.txt");
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
    bool first=true;
    double energy_old;
    while(in.good())
    {
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        first=false;

        precision = precisionOld;
        energy_old =-1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);

        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);
        brems->EnableDEdxInterpolation();

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            dEdx_new=brems->CalculatedEdx();

            if(!particleName.compare("tau") && energy < 10001)precision = 0.5;
            if(!particleName.compare("e") && energy > 1E10)precision = 0.5;

            ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

            in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;

            precision = precisionOld;

        }



        delete cuts;
        delete medium;
        delete particle;
        delete brems;



    }
}

TEST(Bremsstrahlung , Test_of_dNdx_Interpolant ) {
    ifstream in;
    in.open("bin/Brems_dNdx_interpol.txt");

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

       CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);
        brems->EnableDNdxInterpolation();

        //cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << med << "\t" << particleName<< "\t" << dNdx << endl;


        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            dNdx_new=brems->CalculatedNdx();

            //if(dNdx!=0){
            //    if(log10(fabs(dNdx_new -dNdx)/dNdx) > -3){
            //        cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << med << "\t" << particleName<< "\t" << dNdx << endl;
            //        cout << fabs(dNdx_new -dNdx)/dNdx << endl;
            //    }
            //}


            ASSERT_NEAR(dNdx_new, dNdx, 1E-14*dNdx);

            in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        }



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
