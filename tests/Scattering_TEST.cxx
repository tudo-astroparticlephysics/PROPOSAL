
// #include <iostream>

#include "gtest/gtest.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/scattering/ScatteringNoScattering.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"
#include "PROPOSAL/scattering/ScatteringDefault.h"
// #include "PROPOSAL/Scattering.h"
// #include "PROPOSAL/Ionization.h"
// #include "PROPOSAL/Bremsstrahlung.h"
// #include "PROPOSAL/Epairproduction.h"
// #include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"

using namespace std;
using namespace PROPOSAL;

TEST(Comparison , Comparison_equal )
{
    Particle mu = Particle(MuMinusDef::Get());
    Water water(1.0);

    EnergyCutSettings ecuts;
    Utility::Definition utility_defs;
    Utility utils(MuMinusDef::Get(), water, ecuts, utility_defs);

    ScatteringNoScattering noScat1(mu, water);
    ScatteringNoScattering noScat2(mu, water);

    ScatteringMoliere moliere1(mu, water);
    ScatteringMoliere moliere2(mu, water);

    ScatteringHighland high1(mu, water);
    ScatteringHighland high2(mu, water);

    ScatteringDefault default1(mu, utils);
    ScatteringDefault default2(mu, utils);
    // Scattering A = ScatteringFactory::Get().CreateScattering(ScatteringFactory::Moliere, particle_, utility_)
    // Scattering B(mu);

    EXPECT_TRUE(noScat1==noScat2);
    EXPECT_TRUE(moliere1==moliere2);
    EXPECT_TRUE(high1==high2);
    EXPECT_TRUE(default1==default2);

}

TEST(Comparison , Comparison_not_equal ) {

    Particle mu = Particle(MuMinusDef::Get());
    Particle tau = Particle(TauMinusDef::Get());
    Water water(1.0);

    EnergyCutSettings ecuts;
    Utility::Definition utility_defs;
    Utility utils(MuMinusDef::Get(), water, ecuts, utility_defs);

    ScatteringNoScattering noScat1(mu, water);
    ScatteringNoScattering noScat2(tau, water);

    ScatteringMoliere moliere1(mu, water);
    ScatteringMoliere moliere2(tau, water);

    ScatteringHighland high1(mu, water);
    ScatteringHighland high2(tau, water);

    ScatteringDefault default1(mu, utils);
    ScatteringDefault default2(tau, utils);

    EXPECT_TRUE(noScat1!=noScat2);
    EXPECT_TRUE(moliere1!=moliere2);
    EXPECT_TRUE(high1!=high2);
    EXPECT_TRUE(default1!=default2);


}

// TEST(Assignment , Swap ) {
//     Scattering A;
//     Scattering B;
//     Scattering *C = new Scattering();
//     Scattering *D = new Scattering(*C);

//     EXPECT_TRUE(A==B);
//     EXPECT_TRUE(*C == *D);

//     A.swap(*C);
//     EXPECT_TRUE(A == *D);
//     EXPECT_TRUE(B == *C);
// }

// TEST(Assignment , Copyconstructor ) {
//     Scattering A;
//     Scattering B;

//     A=B;
//     EXPECT_TRUE(A==B);
// }

// TEST(Assignment , Copyconstructor2 ) {
//     Scattering A;
//     Scattering B(A);

//     EXPECT_TRUE(A==B);


// }


TEST(Scattering , Theta0 )
{
    cout << "!!! NOT COMPARABLE: BUG IN OLD VERSION !!!" << endl;
    return;
//    ifstream in;
//    in.open("bin/TestFiles/Scattering_Theta0.txt");

//    char firstLine[256];
//    in.getline(firstLine,256);
//    if(!in.good())
//    {
//        cerr << "File Scattering_Theta0.txt not found!!" << endl;
//        EXPECT_TRUE(false);
//    }

//    double dEdx_new;
//    double energy;
//    double dEdx;
//    double ecut;
//    double vcut;
//    string mediumName;
//    string particleName;
//    bool lpm;
//    int para;

//    cout.precision(16);

//    double dr,ef,Theta0,Theta0_new;
//    double energy_old=-1;
//    bool first = true;
//    while(in.good())
//    {
//        if(first)in>>ecut>>vcut>>lpm>>mediumName>>particleName>>dr>>energy>>ef>>Theta0;
//        first=false;
//        energy_old = -1;

//        Medium *medium = new Medium(Medium::GetTypeFromName(mediumName),1.);
//        Particle *particle = new Particle(PROPOSALParticle::GetTypeFromName(particleName),1.,1.,1,.20,20,1e5,10);
//        particle->SetEnergy(energy);
//        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//        std::vector<CrossSections*> vecOfProcColl;
//        CrossSections* ion = new Ionization(particle,medium,cuts);

//        CrossSections* brems = new Bremsstrahlung(particle,medium,cuts);
//        brems->SetParametrization(1);

//        CrossSections* epair = new Epairproduction(particle,medium,cuts);

//        CrossSections* photo = new Photonuclear(particle,medium,cuts);
//        photo->SetParametrization(12);


//        vecOfProcColl.push_back(ion);
//        vecOfProcColl.push_back(brems);
//        vecOfProcColl.push_back(epair);
//        vecOfProcColl.push_back(photo);

//        Scattering* scat = new Scattering(vecOfProcColl);
//        scat->EnableInterpolation("/data/LocalApps/LocalFiles/tables");

//        while(energy_old<=energy)
//        {
//            energy_old = energy;
//            scat->GetParticle()->SetEnergy(energy);
//            Theta0_new = scat->CalculateTheta0(dr,energy,ef);


//            //if(fabs(Theta0 -  Theta0_new)>1e-4*Theta0)cout << mediumName << "\t" << particleName << "\t" << ecut << "\t" << vcut << endl;

//            EXPECT_NEAR(Theta0, Theta0_new, 1e-2*Theta0);


//            in>>ecut>>vcut>>lpm>>mediumName>>particleName>>dr>>energy>>ef>>Theta0;
//            if(in.good() == false)break;
//        }


//        vecOfProcColl.clear();
//        delete scat;
//        delete medium;
//        delete particle;
//        delete cuts;
//    }
}

TEST(Scattering , Advance ) {
    cout << "!!! NOT COMPARABLE: BUG IN OLD VERSION !!!" << endl;
    return;
//    ifstream in;
//    in.open("bin/TestFiles/Scattering_Advance.txt");

//    char firstLine[256];
//    in.getline(firstLine,256);
//    if(!in.good())
//    {
//        cerr << "File Scattering_Advance.txt not found!!" << endl;
//        EXPECT_TRUE(false);
//    }

//    double dEdx_new;
//    double energy;
//    double dEdx;
//    double ecut;
//    double vcut;
//    string mediumName;
//    string particleName;
//    bool lpm;
//    int para;

//    cout.precision(16);

//    double dr,ef;
//    double x,y,z,costheta,cosphi;
//    double x_new,y_new,z_new,theta_new,phi_new,costheta_new,cosphi_new;
//    double energy_old=-1;
//    bool first = true;
//    while(in.good())
//    {
//        if(first)in>>ecut>>vcut>>lpm>>mediumName>>particleName>>dr>>energy>>ef>>x>>y>>z>>cosphi>>costheta;
//        first=false;
//        energy_old = -1;

//        Medium *medium = new Medium(Medium::GetTypeFromName(mediumName),1.);
//        Particle *particle = new Particle(PROPOSALParticle::GetTypeFromName(particleName),1.,1.,1,.20,20,1e5,10);
//        particle->SetEnergy(energy);
//        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//        std::vector<CrossSections*> vecOfProcColl;
//        CrossSections* ion = new Ionization(particle,medium,cuts);

//        CrossSections* brems = new Bremsstrahlung(particle,medium,cuts);
//        brems->SetParametrization(1);

//        CrossSections* epair = new Epairproduction(particle,medium,cuts);

//        CrossSections* photo = new Photonuclear(particle,medium,cuts);
//        photo->SetParametrization(12);


//        vecOfProcColl.push_back(ion);
//        vecOfProcColl.push_back(brems);
//        vecOfProcColl.push_back(epair);
//        vecOfProcColl.push_back(photo);

//        Scattering* scat = new Scattering(vecOfProcColl);
//        scat->EnableInterpolation("/data/LocalApps/LocalFiles/tables");


//        while(energy_old<=energy)
//        {
//            energy_old = energy;
//            scat->GetParticle()->SetEnergy(energy);

//            scat->GetParticle()->SetPhi(0);
//            scat->GetParticle()->SetTheta(0);
//            scat->GetParticle()->SetX(0);
//            scat->GetParticle()->SetY(0);
//            scat->GetParticle()->SetZ(0);


//            scat->Scatter(dr,energy,ef);

//            phi_new = scat->GetParticle()->GetPhi();
//            theta_new = scat->GetParticle()->GetTheta();
//            x_new = scat->GetParticle()->GetX();
//            y_new = scat->GetParticle()->GetY();
//            z_new = scat->GetParticle()->GetZ();

//            cosphi_new = cos(phi_new);
//            costheta_new = cos(theta_new);

//            //if(fabs(x -  Theta0_new)>1e-4*Theta0)cout << mediumName << "\t" << particleName << "\t" << ecut << "\t" << vcut << endl;

//            ASSERT_NEAR(x, x_new, fabs(1e-4*x));
//            ASSERT_NEAR(y, y_new, fabs(1e-4*y));
//            ASSERT_NEAR(z, z_new, fabs(1e-8*z));
//            ASSERT_NEAR(costheta, costheta_new, fabs(1e-8*costheta));
//            ASSERT_NEAR(cosphi, cosphi_new, fabs(1e-4*cosphi));

//            in>>ecut>>vcut>>lpm>>mediumName>>particleName>>dr>>energy>>ef>>x>>y>>z>>cosphi>>costheta;
//            if(in.good() == false)break;
//        }


//        vecOfProcColl.clear();
//        delete scat;
//        delete medium;
//        delete particle;
//        delete cuts;
//    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
