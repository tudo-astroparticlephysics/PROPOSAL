#include "gtest/gtest.h"
#include "PROPOSAL/Epairproduction.h"
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

TEST(Epairproduction , Test_of_dEdx ) {
    return;
    ifstream in;
    in.open("bin/Epair_dEdx.txt");

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

    cout.precision(16);


    while(in.good())
    {
        in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *epair = new Epairproduction(particle, medium, cuts);

        epair->EnableLpmEffect(lpm);

        dEdx_new=epair->CalculatedEdx();
        ASSERT_NEAR(dEdx_new, dEdx, 1e-5*dEdx);

        delete cuts;
        delete medium;
        delete particle;
        delete epair;
    }
}

TEST(Epairproduction , Test_of_dEdx_interpol ) {
    return;
    ifstream in;
    in.open("bin/Epair_dEdx_interpol.txt");

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

    double precision = 1E-5;

    double energy_old;
    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *epair = new Epairproduction(particle, medium, cuts);

        epair->EnableLpmEffect(lpm);
        epair->EnableDEdxInterpolation();

        while(energy_old < energy)
        {
            energy_old = energy;

            epair->GetParticle()->SetEnergy(energy);
            dEdx_new=epair->CalculatedEdx();

            ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        }
        delete cuts;
        delete medium;
        delete particle;
        delete epair;
    }

}

TEST(Epairproduction , Test_of_dNdx ) {
    return;
    ifstream in;
    in.open("bin/Epair_dNdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdx_new;
    double energy;
    double dNdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;

    double precision = 1E-6;

    double energy_old;
    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *epair = new Epairproduction(particle, medium, cuts);

        epair->EnableLpmEffect(lpm);
        //epair->EnableDNdxInterpolation();


        while(energy_old < energy)
        {
            energy_old = energy;

            epair->GetParticle()->SetEnergy(energy);
            dNdx_new=epair->CalculatedNdx();

            if(dNdx!=0)
                if(log10(fabs(1-dNdx_new/dNdx))>-13)
                {
                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
                    cout << energy << "\t" << log10(fabs(1-dNdx_new/dNdx)) << endl;
                }
            ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        }

        delete cuts;
        delete medium;
        delete particle;
        delete epair;
    }
}

TEST(Epairproduction , Test_of_dNdx_interpol ) {
    return;
    ifstream in;
    in.open("bin/Epair_dNdx_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdx_new;
    double energy;
    double dNdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;

    double precision = 1E-2;

    double energy_old;
    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *epair = new Epairproduction(particle, medium, cuts);

        epair->EnableLpmEffect(lpm);
        epair->EnableDNdxInterpolation();

        while(energy_old < energy)
        {
            energy_old = energy;

            epair->GetParticle()->SetEnergy(energy);
            dNdx_new=epair->CalculatedNdx();

            if(dNdx!=0)
                if(log10(fabs(1-dNdx_new/dNdx))>-8)
                {
                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
                    cout << energy << "\t" << log10(fabs(1-dNdx_new/dNdx)) << endl;
                }

            ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        }
        delete cuts;
        delete medium;
        delete particle;
        delete epair;
    }
}

TEST(Epairproduction , Test_of_dNdxrnd ) {
    return;
    ifstream in;
    in.open("bin/Epair_dNdxrnd.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdxrnd_new;
    double energy;
    double dNdxrnd;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;

    double precision = 1E-6;

    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/rnd.txt");

    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *epair = new Epairproduction(particle, medium, cuts);

        epair->EnableLpmEffect(lpm);
        //epair->EnableDNdxInterpolation();


        while(energy_old < energy)
        {
            energy_old = energy;

            epair->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=epair->CalculatedNdx(Rand->rnd());

            /*
            if(dNdxrnd!=0){
                if(log10(fabs(1-dNdxrnd_new/dNdxrnd))>-14)
                {
                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
                    cout << energy << "\t" << log10(fabs(1-dNdxrnd_new/dNdxrnd)) << endl;
                }
            }
            */
            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, precision*dNdxrnd);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        }

        delete cuts;
        delete medium;
        delete particle;
        delete epair;
    }
    delete Rand;
}

TEST(Epairproduction , Test_of_dNdxrnd_interpol ) {
    return;
    ifstream in;
    in.open("bin/Epair_dNdxrnd_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdxrnd_new;
    double energy;
    double dNdxrnd;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;

    double precision = 1E-2;

    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/rnd.txt");

    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *epair = new Epairproduction(particle, medium, cuts);

        epair->EnableLpmEffect(lpm);
        epair->EnableDNdxInterpolation();


        while(energy_old < energy)
        {
            energy_old = energy;

            epair->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=epair->CalculatedNdx(Rand->rnd());


            if(dNdxrnd!=0){
                if(log10(fabs(1-dNdxrnd_new/dNdxrnd))>-6)
                {
                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
                    cout << energy << "\t" << log10(fabs(1-dNdxrnd_new/dNdxrnd)) << endl;
                }
            }

            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, precision*dNdxrnd);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        }

        delete cuts;
        delete medium;
        delete particle;
        delete epair;
    }
    delete Rand;
}

TEST(Epairproduction , Test_of_e ) {
    return;
    ifstream in;
    in.open("bin/Epair_e.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double e_new;
    double energy;
    double e;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;

    double precision = 1E-6;

    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/rnd.txt");
    RndFromFile* Rand2 = new RndFromFile("bin/rnd.txt");
    Rand2->rnd();
    double rnd1,rnd2;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *epair = new Epairproduction(particle, medium, cuts);

        epair->EnableLpmEffect(lpm);


        while(energy_old < energy)
        {
            energy_old = energy;

            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            epair->GetParticle()->SetEnergy(energy);
            e_new=epair->CalculateStochasticLoss(rnd1,rnd2);

/*
            if(e!=0){
                if(log10(fabs(1-e_new/e))>-8)
                {
                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
                    cout << energy << "\t" << log10(fabs(1-e_new/e)) << endl;
                }
            }
*/
            ASSERT_NEAR(e_new, e, precision*e);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
        }

        delete cuts;
        delete medium;
        delete particle;
        delete epair;
    }
    delete Rand;
    delete Rand2;
}

TEST(Epairproduction , Test_of_e_interpol ) {
    //return;
    ifstream in;
    in.open("bin/Epair_e_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double e_new;
    double energy;
    double e;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;

    double precision = 1E-2;

    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/rnd.txt");
    RndFromFile* Rand2 = new RndFromFile("bin/rnd.txt");
    Rand2->rnd();
    double rnd1,rnd2;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *epair = new Epairproduction(particle, medium, cuts);

        epair->EnableLpmEffect(lpm);
        cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
        epair->EnableDNdxInterpolation();
        cout << "Interpolation table build!" << endl;


        while(energy_old < energy)
        {
            energy_old = energy;

            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            epair->GetParticle()->SetEnergy(energy);
            e_new=epair->CalculateStochasticLoss(rnd1,rnd2);


            if(e!=0){
                if(log10(fabs(1-e_new/e))>-14)
                {
                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
                    cout << energy << "\t" << log10(fabs(1-e_new/e)) << endl;
                }
            }

            ASSERT_NEAR(e_new, e, precision*e);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
        }

        delete cuts;
        delete medium;
        delete particle;
        delete epair;
    }
    delete Rand;
    delete Rand2;
}



int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
