#include "gtest/gtest.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>
#include "PROPOSAL/Output.h"


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
        if(!in_.good())log_warn("less than one rnd_number!");
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

TEST(Comparison , Comparison_equal ) {

    double dEdx;
    Medium *medium = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    Epairproduction *A = new Epairproduction(particle, medium, cuts);
    Epairproduction *B = new Epairproduction(particle, medium, cuts);
    EXPECT_TRUE(*A==*B);
    Epairproduction *C = new Epairproduction();
    Epairproduction *D = new Epairproduction();
    EXPECT_TRUE(*C==*D);

    A->GetParticle()->SetEnergy(1e6);
    B->GetParticle()->SetEnergy(1e6);
    EXPECT_TRUE(*A==*B);

    dEdx = A->CalculatedNdx();
    dEdx = B->CalculatedNdx();
    EXPECT_TRUE(*A==*B);

    A->EnableDEdxInterpolation();
    A->EnableDNdxInterpolation();

    B->EnableDEdxInterpolation();
    B->EnableDNdxInterpolation();

    EXPECT_TRUE(*A==*B);
}

TEST(Comparison , Comparison_not_equal ) {
    double dEdx;
    Medium *medium = new Medium("air",1.);
    Medium *medium2 = new Medium("water",1.);
    Particle *particle = new Particle("mu",1.,1.,1,20,20,1e5,10);
    Particle *particle2 = new Particle("tau",1.,1.,1,20,20,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    Epairproduction *A = new Epairproduction(particle, medium, cuts);
    Epairproduction *B = new Epairproduction(particle, medium2, cuts);
    Epairproduction *C = new Epairproduction(particle2, medium, cuts);
    Epairproduction *D = new Epairproduction(particle2, medium2, cuts);
    Epairproduction *E = new Epairproduction(particle2, medium2, cuts);

    EXPECT_TRUE(*A!=*B);
    EXPECT_TRUE(*C!=*D);
    EXPECT_TRUE(*B!=*D);
    EXPECT_TRUE(*D==*E);

    E->SetParticle(particle);
    EXPECT_TRUE(*D!=*E);
    D->SetParticle(particle);
    EXPECT_TRUE(*D==*E);
    D->SetParametrization(6);
    EXPECT_TRUE(*D!=*E);


}

TEST(Assignment , Copyconstructor ) {
    Epairproduction A;
    Epairproduction B =A;

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Copyconstructor2 ) {
    Medium *medium = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);

    Epairproduction A(particle, medium, cuts);
    Epairproduction B(A);

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Operator ) {
    Medium *medium = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    Epairproduction A(particle, medium, cuts);
    Epairproduction B(particle, medium, cuts);
    A.SetParametrization(6);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);

    Medium *medium2 = new Medium("water",1.);
    Particle *particle2 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts2 = new EnergyCutSettings(200,-1);
    Epairproduction *C = new Epairproduction(particle2, medium2, cuts2);
    EXPECT_TRUE(A!=*C);

    A=*C;

    EXPECT_TRUE(A==*C);

}

TEST(Assignment , Swap ) {
    Medium *medium = new Medium("air",1.);
    Medium *medium2 = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,20,20,1e5,10);
    Particle *particle2 = new Particle("mu",1.,1.,1,20,20,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    EnergyCutSettings *cuts2 = new EnergyCutSettings(500,-1);
    Epairproduction A(particle, medium, cuts);
    Epairproduction B(particle2, medium2, cuts2);
    A.EnableDEdxInterpolation();
    B.EnableDEdxInterpolation();
    EXPECT_TRUE(A==B);

    Medium *medium3 = new Medium("water",1.);
    Medium *medium4 = new Medium("water",1.);
    Particle *particle3 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
    Particle *particle4 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts3 = new EnergyCutSettings(200,-1);
    EnergyCutSettings *cuts4 = new EnergyCutSettings(200,-1);
    Epairproduction *C = new Epairproduction(particle3, medium3, cuts3);
    Epairproduction *D = new Epairproduction(particle4, medium4, cuts4);
    EXPECT_TRUE(*C==*D);

    A.swap(*C);

    EXPECT_TRUE(A==*D);
    EXPECT_TRUE(*C==B);


}


std::vector<Medium*>                CombOfMedium;
std::vector<Particle*>              CombOfParticle;
std::vector<EnergyCutSettings*>     CombOfEnergyCutSettings;
std::vector<CrossSections*>         CombOfEpair;

TEST(Epairproduction , Set_Up ) {
    ifstream in;
    in.open("bin/TestFiles/Epair_dEdx.txt");

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

    int i=-1;
    double energy_old;
    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        energy_old  = -1;

        i++;
        CombOfMedium.push_back(new Medium(med,1.));
        CombOfParticle.push_back(new Particle(particleName,1.,1.,1,.20,20,1e5,10));
        CombOfParticle.at(i)->SetEnergy(energy);
        CombOfEnergyCutSettings.push_back(new EnergyCutSettings(ecut,vcut));
        CombOfEpair.push_back(new Epairproduction(CombOfParticle.at(i), CombOfMedium.at(i), CombOfEnergyCutSettings.at(i)));

        CombOfEpair.at(i)->EnableLpmEffect(lpm);

        while(energy_old < energy)
        {
            energy_old = energy;
            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        }
    }
}

TEST(Epairproduction , Test_of_dEdx ) {
    ifstream in;
    in.open("bin/TestFiles/Epair_dEdx.txt");

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

    int i = -1;
    double energy_old;
    while(in.good())
    {
        in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        energy_old = -1;
        i++;
        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
            dEdx_new=CombOfEpair.at(i)->CalculatedEdx();

            ASSERT_NEAR(dEdx_new, dEdx, 1E-5*dEdx);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        }
    }
}

TEST(Epairproduction , Test_of_dNdx ) {
    ifstream in;
    in.open("bin/TestFiles/Epair_dNdx.txt");

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

    double precision = 1E-5;

    double energy_old;
    bool first = true;
    int i=-1;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        first = false;
        energy_old = -1;

        i++;
        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
            dNdx_new=CombOfEpair.at(i)->CalculatedNdx();

            if(dNdx!=0)
            {
                if(log10(fabs(1-dNdx_new/dNdx))>-13)
                {
//                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
//                    cout << energy << "\t" << log10(fabs(1-dNdx_new/dNdx)) << endl;
                }
            }
            ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        }
    }
}


TEST(Epairproduction , Test_of_dNdxrnd ) {
    ifstream in;
    in.open("bin/TestFiles/Epair_dNdxrnd.txt");

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

    double precision = 1E-5;

    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

    int i=-1;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        first = false;
        energy_old = -1;

        i++;
        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=CombOfEpair.at(i)->CalculatedNdx(Rand->rnd());


            if(dNdxrnd!=0){
                if(log10(fabs(1-dNdxrnd_new/dNdxrnd))>-14)
                {
//                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
//                    cout << energy << "\t" << log10(fabs(1-dNdxrnd_new/dNdxrnd)) << endl;
                }
            }

            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, precision*dNdxrnd);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        }
    }
    delete Rand;
}


TEST(Epairproduction , Test_of_e ) {
    ifstream in;
    in.open("bin/TestFiles/Epair_e.txt");

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

    double precision = 1E-5;

    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
    RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
    Rand2->rnd();
    double rnd1,rnd2;
    int i = -1;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
        first = false;
        energy_old = -1;


        i++;
        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);

            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
            e_new=CombOfEpair.at(i)->CalculateStochasticLoss(rnd1,rnd2);


            if(e!=0){
                if(log10(fabs(1-e_new/e))>-14)
                {
                    //cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
                    //cout << energy << "\t" << log10(fabs(1-e_new/e)) << endl;
                }
            }

            ASSERT_NEAR(e_new, e, precision*e);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
        }
    }
    delete Rand;
    delete Rand2;
}

TEST(Epairproduction , Test_of_dEdx_interpol ) {
    ifstream in;
    in.open("bin/TestFiles/Epair_dEdx_interpol.txt");

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
    int i = -1;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        first = false;
        energy_old = -1;


        i++;
        CombOfEpair.at(i)->EnableDEdxInterpolation();
        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
            dEdx_new=CombOfEpair.at(i)->CalculatedEdx();

            ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        }
    }
}

TEST(Epairproduction , Test_of_dNdx_interpol ) {
    ifstream in;
    in.open("bin/TestFiles/Epair_dNdx_interpol.txt");

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

    int i=-1;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        first = false;
        energy_old = -1;

        i++;
        CombOfEpair.at(i)->EnableDNdxInterpolation();

        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
            dNdx_new=CombOfEpair.at(i)->CalculatedNdx();

            if(dNdx!=0)
                if(log10(fabs(1-dNdx_new/dNdx))>-8)
                {
                    //cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
                    //cout << energy << "\t" << log10(fabs(1-dNdx_new/dNdx)) << endl;
                }

            ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        }
    }
}

TEST(Epairproduction , Test_of_dNdxrnd_interpol ) {
    ifstream in;
    in.open("bin/TestFiles/Epair_dNdxrnd_interpol.txt");

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
    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

    int i=-1;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        first = false;
        energy_old = -1;

        i++;
        CombOfEpair.at(i)->EnableDNdxInterpolation();
        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=CombOfEpair.at(i)->CalculatedNdx(Rand->rnd());


            if(dNdxrnd!=0){
                if(log10(fabs(1-dNdxrnd_new/dNdxrnd))>-5)
                {
                    log_warn("%f \t %f \t %i \t %f \t %s \t %s \t",ecut,vcut,lpm,energy,med.c_str(),particleName.c_str());
                    log_warn("%f \t %f",energy,log10(fabs(1-dNdxrnd_new/dNdxrnd)));
                }
            }

            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, precision*dNdxrnd);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        }
    }
    delete Rand;
}

TEST(Epairproduction , Test_of_e_interpol ) {
    ifstream in;
    in.open("bin/TestFiles/Epair_e_interpol.txt");

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
    double precision_old=precision;
    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
    RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
    Rand2->rnd();
    double rnd1,rnd2;
    int i=-1;
    int ctr=0;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
        first = false;
        energy_old = -1;

        i++;
        CombOfEpair.at(i)->EnableDEdxInterpolation();
        CombOfEpair.at(i)->EnableDNdxInterpolation();

        while(energy_old < energy)
        {
            precision=precision_old;
            energy_old = energy;

            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            //The Cross section for such high energy is just extrapolated and
            //therefore pretty unceartain.
            if(!particleName.compare("e") && energy > 1E12)precision = 1E-1;

            CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
            e_new=CombOfEpair.at(i)->CalculateStochasticLoss(rnd1,rnd2);


            if(e!=0){
                if(log10(fabs(1-e_new/e))>-3)
                {
                    //cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
                    cout << energy << "\t" << log10(fabs(1-e_new/e)) << endl;
                }
            }
            //cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
            //cout << "ctr: " << ctr++ << endl;
            EXPECT_NEAR(e_new, e, precision*e);

            in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
        }
    }
    delete Rand;
    delete Rand2;
}



int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
