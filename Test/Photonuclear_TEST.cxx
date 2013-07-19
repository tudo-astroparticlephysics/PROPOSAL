#include "gtest/gtest.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>
#include "PROPOSAL/Output.h"

using namespace std;

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
    Photonuclear *A = new Photonuclear(particle, medium, cuts);
    Photonuclear *B = new Photonuclear(particle, medium, cuts);
    EXPECT_TRUE(*A==*B);

    Photonuclear *C = new Photonuclear();
    Photonuclear *D = new Photonuclear();
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
    Photonuclear *A = new Photonuclear(particle, medium, cuts);
    Photonuclear *B = new Photonuclear(particle, medium2, cuts);
    Photonuclear *C = new Photonuclear(particle2, medium, cuts);
    Photonuclear *D = new Photonuclear(particle2, medium2, cuts);
    Photonuclear *E = new Photonuclear(particle2, medium2, cuts);

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
    Photonuclear A;
    Photonuclear B =A;

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Copyconstructor2 ) {
    Medium *medium = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);

    Photonuclear A(particle, medium, cuts);
    Photonuclear B(A);

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Operator ) {
    Medium *medium = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    Photonuclear A(particle, medium, cuts);
    Photonuclear B(particle, medium, cuts);
    A.SetParametrization(6);
    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);

    Medium *medium2 = new Medium("water",1.);
    Particle *particle2 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts2 = new EnergyCutSettings(200,-1);
    Photonuclear *C = new Photonuclear(particle2, medium2, cuts2);
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
    Photonuclear A(particle, medium, cuts);
    Photonuclear B(particle2, medium2, cuts2);
    A.SetParametrization(2);
    B.SetParametrization(2);
    A.EnableDEdxInterpolation();
    B.EnableDEdxInterpolation();
    A.EnableDNdxInterpolation();
    B.EnableDNdxInterpolation();
    EXPECT_TRUE(A==B);

    Medium *medium3 = new Medium("water",1.);
    Medium *medium4 = new Medium("water",1.);
    Particle *particle3 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
    Particle *particle4 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts3 = new EnergyCutSettings(200,-1);
    EnergyCutSettings *cuts4 = new EnergyCutSettings(200,-1);
    Photonuclear *C = new Photonuclear(particle3, medium3, cuts3);
    Photonuclear *D = new Photonuclear(particle4, medium4, cuts4);
    EXPECT_TRUE(*C==*D);

    A.swap(*C);

    EXPECT_TRUE(A==*D);
    EXPECT_TRUE(*C==B);


}
//----------------------------------------------------------------------------//

std::vector<Medium*>                CombOfMedium;
std::vector<Particle*>              CombOfParticle;
std::vector<EnergyCutSettings*>     CombOfEnergyCutSettings;
std::vector<Photonuclear*>          CombOfPhoto;

//----------------------------------------------------------------------------//

TEST(Photonuclear , Set_Up ) {
    ifstream in;
    in.open("bin/TestFiles/Photo_dEdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    int shadow;
    int bb;
    int para;

    cout.precision(16);

    int i=-1;
    double energy_old;
    bool first = true;
    while(in.good())
    {
        if(first)in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dEdx;
        energy_old  = -1;

        i++;
        CombOfMedium.push_back(new Medium(med,1.));
        CombOfParticle.push_back(new Particle(particleName,1.,1.,1,.20,20,1e5,10));
        CombOfParticle.at(i)->SetEnergy(energy);
        CombOfEnergyCutSettings.push_back(new EnergyCutSettings(ecut,vcut));
        CombOfPhoto.push_back(new Photonuclear(CombOfParticle.at(i), CombOfMedium.at(i), CombOfEnergyCutSettings.at(i)));

        // Now: parametrization_ = 1,  Former: form=1 and bb=1 Kokoulin
        // Now: parametrization_ = 2,  Former: form=2 and bb=1 Kokoulin + hard component
        // Now: parametrization_ = 3,  Former: form=1 and bb=2 Rhode
        // Now: parametrization_ = 4,  Former: form=2 and bb=2 Rhode + hard component
        // Now: parametrization_ = 5,  Former: form=1 and bb=3 Bezrukov/Bugaev
        // Now: parametrization_ = 6,  Former: form=2 and bb=3 Bezrukov/Bugaev + hard component
        // Now: parametrization_ = 7,  Former: form=1 and bb=4 Zeus
        // Now: parametrization_ = 8,  Former: form=2 and bb=4 Zeus + hard component
        // Now: parametrization_ = 9,  Former: form=3 and bb=1 shadow=1 ALLM 91
        // Now: parametrization_ = 10, Former: form=3 and bb=1 shadow=2 ALLM 91
        // Now: parametrization_ = 11, Former: form=3 and bb=2 shadow=1 ALLM 97
        // Now: parametrization_ = 12, Former: form=3 and bb=2 shadow=2 ALLM 97
        // Now: parametrization_ = 13, Former: form=4 and bb=1 shadow=1 Butkevich/Mikhailov
        // Now: parametrization_ = 14, Former: form=4 and bb=1 shadow=2 Butkevich/Mikhailov

        if(bb==1&&para==1)CombOfPhoto.at(i)->SetParametrization(1);
        if(bb==1&&para==2)CombOfPhoto.at(i)->SetParametrization(2);
        if(bb==2&&para==1)CombOfPhoto.at(i)->SetParametrization(3);
        if(bb==2&&para==2)CombOfPhoto.at(i)->SetParametrization(4);
        if(bb==3&&para==1)CombOfPhoto.at(i)->SetParametrization(5);
        if(bb==3&&para==2)CombOfPhoto.at(i)->SetParametrization(6);
        if(bb==4&&para==1)CombOfPhoto.at(i)->SetParametrization(7);
        if(bb==4&&para==2)CombOfPhoto.at(i)->SetParametrization(8);
        if(bb==1&&para==3&&shadow==1)CombOfPhoto.at(i)->SetParametrization(9);
        if(bb==1&&para==3&&shadow==2)CombOfPhoto.at(i)->SetParametrization(10);
        if(bb==2&&para==3&&shadow==1)CombOfPhoto.at(i)->SetParametrization(11);
        if(bb==2&&para==3&&shadow==2)CombOfPhoto.at(i)->SetParametrization(12);
        if(bb==1&&para==4&&shadow==1)CombOfPhoto.at(i)->SetParametrization(13);
        if(bb==1&&para==4&&shadow==2)CombOfPhoto.at(i)->SetParametrization(14);

        while(energy_old < energy)
        {
            energy_old = energy;
            in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dEdx;
        }
    }
}

//----------------------------------------------------------------------------//

TEST(Photonuclear , Test_of_e ) {
    ifstream in;
    in.open("bin/TestFiles/Photo_e.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double e_new;
    double energy;
    double e;
    double ecut;
    double vcut;
    string med;
    string particleName;
    int shadow;
    int bb;
    int para;
    double precision = 1E-6;

    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
    RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
    Rand2->rnd();
    double rnd1,rnd2;
    int i = -1;
    while(in.good())
    {
        if(first)in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>e;
        first = false;
        energy_old = -1;


        i++;

        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);

            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
            e_new=CombOfPhoto.at(i)->CalculateStochasticLoss(rnd1,rnd2);


            if(e!=0){
                if(log10(fabs(1-e_new/e))>-8)
                {
//                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
//                    cout << energy << "\t" << log10(fabs(1-e_new/e)) << endl;
                }
            }

            ASSERT_NEAR(e_new, e, precision*e);

            in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>e;
        }
    }
    delete Rand;
    delete Rand2;
}

//----------------------------------------------------------------------------//

TEST(Photonuclear , Test_of_dNdxrnd ) {
    ifstream in;
    in.open("bin/TestFiles/Photo_dNdxrnd.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdxrnd_new;
    double energy;
    double dNdxrnd;
    double ecut;
    double vcut;
    string med;
    string particleName;
    int shadow;
    int bb;
    int para;
    double precision = 1E-6;

    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

    int i=-1;
    while(in.good())
    {
        if(first)in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dNdxrnd;
        first = false;
        energy_old = -1;

        i++;

        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=CombOfPhoto.at(i)->CalculatedNdx(Rand->rnd());


            if(dNdxrnd!=0){
                if(log10(fabs(1-dNdxrnd_new/dNdxrnd))>-14)
                {
//                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
//                    cout << energy << "\t" << log10(fabs(1-dNdxrnd_new/dNdxrnd)) << endl;
                }
            }

            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, precision*dNdxrnd);

            in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dNdxrnd;
        }
    }
    delete Rand;
}

//----------------------------------------------------------------------------//


TEST(Photonuclear , Test_of_dNdx ) {
    ifstream in;
    in.open("bin/TestFiles/Photo_dNdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdx_new;
    double energy;
    double dNdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    int shadow;
    int bb;
    int para;

    double precision = 1E-6;

    double energy_old;
    bool first = true;
    int i=-1;
    while(in.good())
    {
        if(first)in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dNdx;
        first = false;
        energy_old = -1;

        i++;

        while(energy_old < energy)
        {
            energy_old = energy;
            CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
            dNdx_new=CombOfPhoto.at(i)->CalculatedNdx();

            if(dNdx!=0)
            {
                if(log10(fabs(1-dNdx_new/dNdx))>-13)
                {
//                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
//                    cout << energy << "\t" << log10(fabs(1-dNdx_new/dNdx)) << endl;
                }
            }
            ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

            in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dNdx;
        }
    }
}

//----------------------------------------------------------------------------//


TEST(Photonuclear , Test_of_dEdx ) {
    ifstream in;
    in.open("bin/TestFiles/Photo_dEdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    int shadow;
    int bb;
    int para;
    double dEdx_new;

    cout.precision(16);

    int i = -1;
    double energy_old;
    while(in.good())
    {
        in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dEdx;
        energy_old = -1;
        i++;
        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
            dEdx_new=CombOfPhoto.at(i)->CalculatedEdx();

            ASSERT_NEAR(dEdx_new, dEdx, 1E-6*dEdx);

            in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dEdx;
        }
    }
}

//----------------------------------------------------------------------------//


TEST(Photonuclear , Test_of_dEdx_interpol ) {
    ifstream in;
    in.open("bin/TestFiles/Photo_dEdx_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dEdx_new;
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    int shadow;
    int bb;
    int para;
    double precision = 1E-5;

    double energy_old;
    bool first = true;
    int i = -1;
    while(in.good())
    {
        if(first)in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dEdx;
        first = false;
        energy_old = -1;


        i++;

        CombOfPhoto.at(i)->EnableDEdxInterpolation();
        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
            dEdx_new=CombOfPhoto.at(i)->CalculatedEdx();

            ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

            in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dEdx;
        }
    }
}

//----------------------------------------------------------------------------//

TEST(Photonuclear , Test_of_dNdx_interpol ) {
    //return;
    ifstream in;
    in.open("bin/TestFiles/Photo_dNdx_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdx_new;
    double energy;
    double dNdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    int shadow;
    int bb;
    int para;
    double precision = 1E-2;

    double energy_old;
    bool first = true;

    int i=-1;
    while(in.good())
    {
        if(first)in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dNdx;
        first = false;
        energy_old = -1;

        i++;

        CombOfPhoto.at(i)->EnableDNdxInterpolation();

        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
            dNdx_new=CombOfPhoto.at(i)->CalculatedNdx();

            if(dNdx!=0)
                if(log10(fabs(1-dNdx_new/dNdx))>-8)
                {
//                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
//                    cout << energy << "\t" << log10(fabs(1-dNdx_new/dNdx)) << endl;
                }

            ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

            in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dNdx;
        }
    }
}

//----------------------------------------------------------------------------//

TEST(Photonuclear , Test_of_dNdxrnd_interpol ) {
    ifstream in;
    in.open("bin/TestFiles/Photo_dNdxrnd_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdxrnd_new;
    double energy;
    double dNdxrnd;
    double ecut;
    double vcut;
    string med;
    string particleName;
    int shadow;
    int bb;
    int para;
    double precision = 1E-2;

    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

    int i=-1;
    while(in.good())
    {
        if(first)in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dNdxrnd;
        first = false;
        energy_old = -1;

        i++;

        CombOfPhoto.at(i)->EnableDNdxInterpolation();
        while(energy_old < energy)
        {
            energy_old = energy;

            CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=CombOfPhoto.at(i)->CalculatedNdx(Rand->rnd());


            if(dNdxrnd!=0){
                if(log10(fabs(1-dNdxrnd_new/dNdxrnd))>-6)
                {
//                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
//                    cout << energy << "\t" << log10(fabs(1-dNdxrnd_new/dNdxrnd)) << endl;
                }
            }

            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, precision*dNdxrnd);

            in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>dNdxrnd;
        }
    }
    delete Rand;
}

//----------------------------------------------------------------------------//

TEST(Photonuclear , Test_of_e_interpol ) {
    ifstream in;
    in.open("bin/TestFiles/Photo_e_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double e_new;
    double energy;
    double e;
    double ecut;
    double vcut;
    string med;
    string particleName;
    int shadow;
    int bb;
    int para;
    double precision = 1E-2;

    double energy_old;
    bool first = true;
    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
    RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
    Rand2->rnd();
    double rnd1,rnd2;
    int i=-1;
    while(in.good())
    {
        if(first)in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>e;
        first = false;
        energy_old = -1;

        i++;

        CombOfPhoto.at(i)->EnableDNdxInterpolation();

        while(energy_old < energy)
        {
            energy_old = energy;

            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
            e_new=CombOfPhoto.at(i)->CalculateStochasticLoss(rnd1,rnd2);


            if(e!=0){
                if(log10(fabs(1-e_new/e))>-14)
                {
//                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << med<< "\t" << particleName << endl;
//                    cout << energy << "\t" << log10(fabs(1-e_new/e)) << endl;
                }
            }

            ASSERT_NEAR(e_new, e, precision*e);

            in>>para>>bb>>shadow>>ecut>>vcut>>energy>>med>>particleName>>e;
        }
    }
    delete Rand;
    delete Rand2;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
