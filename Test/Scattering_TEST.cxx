#include "gtest/gtest.h"
#include "PROPOSAL/Scattering.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Photonuclear.h"
#include <iostream>

using namespace std;

TEST(Comparison , Comparison_equal ) {
    Scattering A;
    Scattering B;

    EXPECT_TRUE(A==B);


    EnergyCutSettings* cut = new EnergyCutSettings(100,0.01);
    Medium* med = new Medium("uranium",1.);
    Particle* particle2 = new Particle("tau",-1,0,3,0,2,0,0);

    vector<CrossSections*> crosssections2;
    crosssections2.push_back(new Ionization(particle2,med,cut));
    crosssections2.push_back(new Bremsstrahlung(particle2,med,cut));
    crosssections2.push_back(new Photonuclear(particle2,med,cut));
    crosssections2.push_back(new Epairproduction(particle2,med,cut));
    Scattering* C = new Scattering(crosssections2);
    Scattering* D = new Scattering(crosssections2);


    EXPECT_TRUE(*C == *D);

    vector<CrossSections*> crosssections;
    crosssections.push_back(new Ionization());
    crosssections.push_back(new Bremsstrahlung());
    crosssections.push_back(new Photonuclear());
    crosssections.push_back(new Epairproduction());

    Ionization* asdsad = new Ionization();
    Ionization* gg = new Ionization();

    Particle* particle = new Particle("mu");

    Scattering* E = new Scattering(crosssections);
    E->SetParticle(particle);
    EXPECT_TRUE(A==*E);
}

TEST(Comparison , Comparison_not_equal ) {
    EnergyCutSettings* cut = new EnergyCutSettings(100,0.01);
    Medium* med = new Medium("uranium",1.);
    Particle* particle2 = new Particle("tau",-1,0,3,0,2,0,0);

    vector<CrossSections*> crosssections2;
    crosssections2.push_back(new Ionization(particle2,med,cut));
    crosssections2.push_back(new Bremsstrahlung(particle2,med,cut));
    crosssections2.push_back(new Photonuclear(particle2,med,cut));
    crosssections2.push_back(new Epairproduction(particle2,med,cut));

    Scattering A;
    Scattering B(crosssections2);

    EXPECT_TRUE(A!=B);

    vector<CrossSections*> crosssections;
    crosssections.push_back(new Ionization());
    crosssections.push_back(new Bremsstrahlung());
    crosssections.push_back(new Photonuclear());
    crosssections.push_back(new Epairproduction());

    Particle* particle = new Particle("mu",0,0,0,0,0,0,0);
    for(unsigned int i =0;i<crosssections.size();i++)
    {
        crosssections.at(i)->SetParticle(particle);
    }

    Scattering* C = new Scattering(crosssections2);
    Scattering* D = new Scattering(crosssections);
    Scattering* E = new Scattering(crosssections);

    EXPECT_TRUE(*C != *D);
    EXPECT_TRUE(*D == *E);
    E->SetParticle(particle2);
    EXPECT_TRUE(*D!=*E);
    D->SetParticle(particle2);
    EXPECT_TRUE(*D==*E);

}

TEST(Assignment , Swap ) {
    EnergyCutSettings* cut = new EnergyCutSettings(100,0.01);
    Medium* med = new Medium("uranium",1.);
    Particle* particle2 = new Particle("tau",-1,0,3,0,2,0,0);

    vector<CrossSections*> crosssections2;
    crosssections2.push_back(new Ionization(particle2,med,cut));
    crosssections2.push_back(new Bremsstrahlung(particle2,med,cut));
    crosssections2.push_back(new Photonuclear(particle2,med,cut));
    crosssections2.push_back(new Epairproduction(particle2,med,cut));

    Scattering A;
    Scattering B;
    Scattering *C = new Scattering(crosssections2);
    Scattering *D = new Scattering(crosssections2);

    EXPECT_TRUE(A==B);
    EXPECT_TRUE(*C == *D);

    A.swap(*C);
    EXPECT_TRUE(A == *D);
    EXPECT_TRUE(B == *C);
}

TEST(Assignment , Copyconstructor ) {
    Scattering A;
    Scattering B;

    EXPECT_TRUE(A==B);
}

TEST(Assignment , Copyconstructor2 ) {
    EnergyCutSettings* cut = new EnergyCutSettings(100,0.01);
    Medium* med = new Medium("uranium",1.);
    Particle* particle2 = new Particle("tau",-1,0,3,0,2,0,0);

    vector<CrossSections*> crosssections2;
    crosssections2.push_back(new Ionization(particle2,med,cut));
    crosssections2.push_back(new Bremsstrahlung(particle2,med,cut));
    crosssections2.push_back(new Photonuclear(particle2,med,cut));
    crosssections2.push_back(new Epairproduction(particle2,med,cut));

    Scattering A(crosssections2);
    Scattering B(A);

    EXPECT_TRUE(A==B);


}

TEST(Assignment , Operator ) {
    EnergyCutSettings* cut = new EnergyCutSettings(100,0.01);
    Medium* med = new Medium("uranium",1.);
    Particle* particle2 = new Particle("tau",-1,0,3,0,2,0,0);

    vector<CrossSections*> crosssections2;
    crosssections2.push_back(new Ionization(particle2,med,cut));
    crosssections2.push_back(new Bremsstrahlung(particle2,med,cut));
    crosssections2.push_back(new Photonuclear(particle2,med,cut));
    crosssections2.push_back(new Epairproduction(particle2,med,cut));

    Scattering A;
    Scattering B(crosssections2);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);

    A.SetCrossSections(crosssections2);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);
}


TEST(Scattering , Theta0 ) {

    ifstream in;
    in.open("bin/TestFiles/Scattering_Theta0.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    if(!in.good())
    {
        cerr << "File Scattering_Theta0.txt not found!!" << endl;
        EXPECT_TRUE(false);
    }

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

    double dr,ef,Theta0,Theta0_new;
    double energy_old=-1;
    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>med>>particleName>>dr>>energy>>ef>>Theta0;
        first=false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        std::vector<CrossSections*> vecOfProcColl;
        CrossSections* ion = new Ionization(particle,medium,cuts);

        CrossSections* brems = new Bremsstrahlung(particle,medium,cuts);
        brems->SetParametrization(1);

        CrossSections* epair = new Epairproduction(particle,medium,cuts);

        CrossSections* photo = new Photonuclear(particle,medium,cuts);
        photo->SetParametrization(12);


        vecOfProcColl.push_back(ion);
        vecOfProcColl.push_back(brems);
        vecOfProcColl.push_back(epair);
        vecOfProcColl.push_back(photo);

        Scattering* scat = new Scattering(vecOfProcColl);
        scat->EnableInterpolation("/data/LocalApps/LocalFiles/tables");

        while(energy_old<=energy)
        {
            energy_old = energy;
            scat->GetParticle()->SetEnergy(energy);
            Theta0_new = scat->CalculateTheta0(dr,energy,ef);


            if(fabs(Theta0 -  Theta0_new)>1e-4*Theta0)cout << med << "\t" << particleName << "\t" << ecut << "\t" << vcut << endl;

            EXPECT_NEAR(Theta0, Theta0_new, 1e-2*Theta0);


            in>>ecut>>vcut>>lpm>>med>>particleName>>dr>>energy>>ef>>Theta0;
            if(in.good() == false)break;
        }


        vecOfProcColl.clear();
        delete scat;
        delete medium;
        delete particle;
        delete cuts;
    }
}

TEST(Scattering , Advance ) {

    ifstream in;
    in.open("bin/TestFiles/Scattering_Advance.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    if(!in.good())
    {
        cerr << "File Scattering_Advance.txt not found!!" << endl;
        EXPECT_TRUE(false);
    }

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

    double dr,ef;
    double x,y,z,costheta,cosphi;
    double x_new,y_new,z_new,theta_new,phi_new,costheta_new,cosphi_new;
    double energy_old=-1;
    bool first = true;
    while(in.good())
    {
        if(first)in>>ecut>>vcut>>lpm>>med>>particleName>>dr>>energy>>ef>>x>>y>>z>>cosphi>>costheta;
        first=false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        std::vector<CrossSections*> vecOfProcColl;
        CrossSections* ion = new Ionization(particle,medium,cuts);

        CrossSections* brems = new Bremsstrahlung(particle,medium,cuts);
        brems->SetParametrization(1);

        CrossSections* epair = new Epairproduction(particle,medium,cuts);

        CrossSections* photo = new Photonuclear(particle,medium,cuts);
        photo->SetParametrization(12);


        vecOfProcColl.push_back(ion);
        vecOfProcColl.push_back(brems);
        vecOfProcColl.push_back(epair);
        vecOfProcColl.push_back(photo);

        Scattering* scat = new Scattering(vecOfProcColl);
        scat->EnableInterpolation("/data/LocalApps/LocalFiles/tables");


        while(energy_old<=energy)
        {
            energy_old = energy;
            scat->GetParticle()->SetEnergy(energy);

            scat->GetParticle()->SetPhi(0);
            scat->GetParticle()->SetTheta(0);
            scat->GetParticle()->SetX(0);
            scat->GetParticle()->SetY(0);
            scat->GetParticle()->SetZ(0);


            scat->Scatter(dr,energy,ef);

            phi_new = scat->GetParticle()->GetPhi();
            theta_new = scat->GetParticle()->GetTheta();
            x_new = scat->GetParticle()->GetX();
            y_new = scat->GetParticle()->GetY();
            z_new = scat->GetParticle()->GetZ();

            cosphi_new = cos(phi_new);
            costheta_new = cos(theta_new);

            //if(fabs(x -  Theta0_new)>1e-4*Theta0)cout << med << "\t" << particleName << "\t" << ecut << "\t" << vcut << endl;

            ASSERT_NEAR(x, x_new, fabs(1e-4*x));
            ASSERT_NEAR(y, y_new, fabs(1e-4*y));
            ASSERT_NEAR(z, z_new, fabs(1e-8*z));
            ASSERT_NEAR(costheta, costheta_new, fabs(1e-8*costheta));
            ASSERT_NEAR(cosphi, cosphi_new, fabs(1e-4*cosphi));

            in>>ecut>>vcut>>lpm>>med>>particleName>>dr>>energy>>ef>>x>>y>>z>>cosphi>>costheta;
            if(in.good() == false)break;
        }


        vecOfProcColl.clear();
        delete scat;
        delete medium;
        delete particle;
        delete cuts;
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
