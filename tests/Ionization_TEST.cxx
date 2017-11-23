
// #include <iostream>
// #include <string>
// #include <math.h>

#include "gtest/gtest.h"

#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;
using namespace std;

// class RndFromFile{
// private:
//     double rnd_;
//     string Path_;
//     ifstream in_;

// public:
//     RndFromFile(string Path){
//         Path_ = Path;
//         in_.open(Path_.c_str());
//         in_>>rnd_;
//         if(!in_.good())log_warn("less than one rnd_number!");
//     }

//     double rnd(){
//         in_>>rnd_;
//         if(!in_.good()){
//             in_.close();
//             in_.clear();
//             in_.open(Path_.c_str());
//             in_>>rnd_;
//         }
//         return rnd_;
//     }
// };

TEST(Comparison, Comparison_equal)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;

    Ionization* Ioniz_A = new Ionization(particle_def, medium, ecuts, multiplier);
    Parametrization* Ioniz_B = new Ionization(particle_def, medium, ecuts, multiplier);
    EXPECT_TRUE(*Ioniz_A == *Ioniz_B);

    Ionization param(particle_def, medium, ecuts, multiplier);
    EXPECT_TRUE(param == *Ioniz_A);

    IonizIntegral* Int_A = new IonizIntegral(param);
    CrossSectionIntegral* Int_B = new IonizIntegral(param);
    EXPECT_TRUE(*Int_A == *Int_B);

    InterpolationDef InterpolDef;
    IonizInterpolant* Interpol_A = new IonizInterpolant(param, InterpolDef);
    CrossSectionInterpolant* Interpol_B = new IonizInterpolant(param, InterpolDef);
    EXPECT_TRUE(*Interpol_A == *Interpol_B);

    delete Ioniz_A;
    delete Ioniz_B;
    delete Int_A;
    delete Int_B;
    delete Interpol_A;
    delete Interpol_B;
}

TEST(Comparison, Comparison_not_equal)
{
    ParticleDef mu_def = MuMinusDef::Get();
    ParticleDef tau_def = TauMinusDef::Get();
    Water medium_1;
    Ice medium_2;
    EnergyCutSettings ecuts_1(500, -1);
    EnergyCutSettings ecuts_2(-1, 0.05);
    double multiplier_1 = 1.;
    double multiplier_2 = 2.;

    Ionization Ioniz_A(mu_def, medium_1, ecuts_1, multiplier_1);
    Ionization Ioniz_B(tau_def, medium_1, ecuts_1, multiplier_1);
    Ionization Ioniz_C(mu_def, medium_2, ecuts_1, multiplier_1);
    Ionization Ioniz_D(mu_def, medium_1, ecuts_2, multiplier_1);
    Ionization Ioniz_E(mu_def, medium_1, ecuts_1, multiplier_2);
    EXPECT_TRUE(Ioniz_A != Ioniz_B);
    EXPECT_TRUE(Ioniz_A != Ioniz_C);
    EXPECT_TRUE(Ioniz_A != Ioniz_D);
    EXPECT_TRUE(Ioniz_A != Ioniz_E);

    IonizIntegral Int_A(Ioniz_A);
    IonizIntegral Int_B(Ioniz_B);
    EXPECT_TRUE(Int_A != Int_B);

    InterpolationDef InterpolDef;
    IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
    IonizInterpolant Interpol_B(Ioniz_B, InterpolDef);
    EXPECT_TRUE(Interpol_A != Interpol_B);
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef mu_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts(500, -1);
    double multiplier = 1.;

    Ionization Ioniz_A(mu_def, medium, ecuts, multiplier);
    Ionization Ioniz_B = Ioniz_A;

    IonizIntegral Int_A(Ioniz_A);
    IonizIntegral Int_B = Int_A;
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
    IonizInterpolant Interpol_B = Interpol_A;
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

TEST(Assignment, Copyconstructor2)
{

    ParticleDef mu_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts(500, -1);
    double multiplier = 1.;

    Ionization Ioniz_A(mu_def, medium, ecuts, multiplier);
    Ionization Ioniz_B(Ioniz_A);

    IonizIntegral Int_A(Ioniz_A);
    IonizIntegral Int_B(Int_A);
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
    IonizInterpolant Interpol_B(Interpol_A);
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

// in polymorphism an assignmant and swap operator doesn't make sense


// TEST(Ionization , Test_of_dEdx ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Ioniz_dEdx.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dEdx_new;
//     double energy;
//     double dEdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;

//     cout.precision(16);


//     while(in.good())
//     {
//         in>>ecut>>vcut>>energy>>mediumName>>particleName>>dEdx;

//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//         CrossSections *ioniz = new Ionization(particle, medium, cuts);

//         dEdx_new=ioniz->CalculatedEdx();
//         ASSERT_NEAR(dEdx_new, dEdx, 1e-12*dEdx);

//         delete cuts;
//         delete medium;
//         delete particle;
//         delete ioniz;
//     }
// }

// TEST(Ionization , Test_of_dNdx ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Ioniz_dNdx.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdx_new;
//     double energy;
//     double dNdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;

//     cout.precision(16);


//     while(in.good())
//     {
//         in>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdx;

//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//         CrossSections *ioniz = new Ionization(particle, medium, cuts);

//         dNdx_new=ioniz->CalculatedNdx();
//         ASSERT_NEAR(dNdx_new, dNdx, 1e-8*dNdx);

//         delete cuts;
//         delete medium;
//         delete particle;
//         delete ioniz;

//     }
// }

// TEST(Ionization , Test_of_dNdxrnd ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Ioniz_dNdxrnd.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);

//     double dNdxrnd;
//     double dNdxrnd_new;
//     double energy;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;
//     int para;

//     cout.precision(16);
//     double energy_old=-1;

//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdxrnd;
//         first = false;
//         energy_old = -1;
//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//         CrossSections *ioniz = new Ionization(particle, medium, cuts);

//         while(energy_old < energy){
//             energy_old = energy;
//             ioniz->GetParticle()->SetEnergy(energy);
//             dNdxrnd_new=ioniz->CalculatedNdx(Rand->rnd());

//             ASSERT_NEAR(dNdxrnd_new, dNdxrnd, 1E-7*dNdxrnd);

//             in>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdxrnd;
//         }



//         delete cuts;
//         delete medium;
//         delete particle;
//         delete ioniz;
//     }
//     delete Rand;
// }

// TEST(Ionization , Test_of_e ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Ioniz_e.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double e;
//     double energy;
//     double e_new;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;

//     cout.precision(16);

//     double precision = 1E-8;

//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
//     RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
//     Rand2->rnd();

//     double rnd1,rnd2;

//     double energy_old;
//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>energy>>mediumName>>particleName>>e;
//         first = false;
//         energy_old = -1;

//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
//         CrossSections *ioniz = new Ionization(particle, medium, cuts);

//         while(energy_old < energy)
//         {
//             energy_old = energy;
//             rnd1 = Rand->rnd();
//             rnd2 = Rand2->rnd();

//             ioniz->GetParticle()->SetEnergy(energy);
//             e_new=ioniz->CalculateStochasticLoss(rnd1,rnd2);

//             ASSERT_NEAR(e_new, e, precision*e);

//             in>>ecut>>vcut>>energy>>mediumName>>particleName>>e;
//         }
//         delete cuts;
//         delete medium;
//         delete particle;
//         delete ioniz;
//     }
// }

// TEST(Ionization , Test_of_dEdx_Interpolant ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Ioniz_dEdx_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dEdx_new;
//     double energy;
//     double dEdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;

//     cout.precision(16);

//     double precision = 1E-8;

//     double energy_old;
//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>energy>>mediumName>>particleName>>dEdx;
//         first = false;
//         energy_old = -1;

//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
//         CrossSections *ioniz = new Ionization(particle, medium, cuts);

//         ioniz->EnableDEdxInterpolation();

//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             ioniz->GetParticle()->SetEnergy(energy);
//             dEdx_new=ioniz->CalculatedEdx();

//             ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

//             in>>ecut>>vcut>>energy>>mediumName>>particleName>>dEdx;
//         }
//         delete cuts;
//         delete medium;
//         delete particle;
//         delete ioniz;
//     }
// }

// TEST(Ionization , Test_of_dNdx_Interpolant ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Ioniz_dNdx_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdx_new;
//     double energy;
//     double dNdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;

//     cout.precision(16);

//     double precision = 1E-5;

//     double energy_old;
//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdx;
//         first = false;
//         energy_old = -1;

//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
//         CrossSections *ioniz = new Ionization(particle, medium, cuts);

//         ioniz->EnableDNdxInterpolation();

//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             ioniz->GetParticle()->SetEnergy(energy);
//             dNdx_new=ioniz->CalculatedNdx();

//             ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

//             in>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdx;
//         }
//         delete cuts;
//         delete medium;
//         delete particle;
//         delete ioniz;
//     }
// }

// TEST(Ionization , Test_of_dNdxrnd_interpol ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Ioniz_dNdxrnd_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);

//     double dNdxrnd;
//     double dNdxrnd_new;
//     double energy;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;

//     cout.precision(16);
//     double energy_old=-1;

//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdxrnd;
//         first = false;
//         energy_old = -1;
//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//         CrossSections *ioniz = new Ionization(particle, medium, cuts);
//         ioniz->EnableDNdxInterpolation();

//         while(energy_old < energy){
//             energy_old = energy;
//             ioniz->GetParticle()->SetEnergy(energy);
//             dNdxrnd_new=ioniz->CalculatedNdx(Rand->rnd());

//             ASSERT_NEAR(dNdxrnd_new, dNdxrnd, 1E-5*dNdxrnd);

//             in>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdxrnd;
//         }



//         delete cuts;
//         delete medium;
//         delete particle;
//         delete ioniz;
//     }
//     delete Rand;
// }


// TEST(Ionization , Test_of_e_interpol ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Ioniz_e_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double e;
//     double energy;
//     double e_new;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;

//     cout.precision(16);

//     double precision = 1E-2;

//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
//     RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
//     Rand2->rnd();

//     double rnd1,rnd2;

//     double energy_old;
//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>energy>>mediumName>>particleName>>e;
//         first = false;
//         energy_old = -1;

//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
//         CrossSections *ioniz = new Ionization(particle, medium, cuts);

//         ioniz->EnableDEdxInterpolation();

//         while(energy_old < energy)
//         {
//             energy_old = energy;
//             rnd1 = Rand->rnd();
//             rnd2 = Rand2->rnd();

//             ioniz->GetParticle()->SetEnergy(energy);

//             e_new=ioniz->CalculateStochasticLoss(rnd1,rnd2);

//             if(e!=0)if(log10(fabs(1-e_new/e))>-3){
// //                cout<< "\t" << ecut<< "\t" << vcut << "\t" << energy<< "\t" << mediumName<< "\t" << particleName<<endl;
// //                cout << log10(fabs(1-e_new/e)) << endl;
//             }

//             ASSERT_NEAR(e_new, e, precision*e);

//             in>>ecut>>vcut>>energy>>mediumName>>particleName>>e;
//         }
//         delete cuts;
//         delete medium;
//         delete particle;
//         delete ioniz;
//     }
// }

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
