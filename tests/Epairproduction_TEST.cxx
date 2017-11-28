
// #include <iostream>
// #include <string>
// #include <math.h>

#include "gtest/gtest.h"

#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/crossection/EpairIntegral.h"
#include "PROPOSAL/crossection/EpairInterpolant.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"


using namespace std;
using namespace PROPOSAL;

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
    bool lpm = true;

    EpairProductionRhoIntegral* EpairInt_A = new EpairProductionRhoIntegral(particle_def, medium, ecuts, multiplier, lpm);
    Parametrization* EpairInt_B = new EpairProductionRhoIntegral(particle_def, medium, ecuts, multiplier, lpm);
    EXPECT_TRUE(*EpairInt_A == *EpairInt_B);

    EpairProductionRhoIntegral param_int(particle_def, medium, ecuts, multiplier, lpm);
    EXPECT_TRUE(param_int == *EpairInt_A);

    EpairIntegral* Int_A = new EpairIntegral(param_int);
    CrossSectionIntegral* Int_B = new EpairIntegral(param_int);
    EXPECT_TRUE(*Int_A == *Int_B);

    InterpolationDef InterpolDef;
    EpairProductionRhoInterpolant* EpairInterpol_A = new EpairProductionRhoInterpolant(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
    Parametrization* EpairInterpol_B = new EpairProductionRhoInterpolant(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
    EXPECT_TRUE(*EpairInterpol_A == *EpairInterpol_B);

    EpairProductionRhoInterpolant param_interpol(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
    EXPECT_TRUE(param_interpol == *EpairInterpol_A);

    EpairInterpolant* Interpol_A = new EpairInterpolant(param_interpol, InterpolDef);
    CrossSectionInterpolant* Interpol_B = new EpairInterpolant(param_interpol, InterpolDef);
    EXPECT_TRUE(*Interpol_A == *Interpol_B);

    delete EpairInt_A;
    delete EpairInt_B;
    delete Int_A;
    delete Int_B;
    delete EpairInterpol_A;
    delete EpairInterpol_B;
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
    bool lpm_1 = true;
    bool lpm_2 = false;

    EpairProductionRhoIntegral EpairInt_A(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    EpairProductionRhoIntegral EpairInt_B(tau_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    EpairProductionRhoIntegral EpairInt_C(mu_def, medium_2, ecuts_1, multiplier_1, lpm_1);
    EpairProductionRhoIntegral EpairInt_D(mu_def, medium_1, ecuts_2, multiplier_1, lpm_1);
    EpairProductionRhoIntegral EpairInt_E(mu_def, medium_1, ecuts_1, multiplier_2, lpm_1);
    EpairProductionRhoIntegral EpairInt_F(mu_def, medium_1, ecuts_1, multiplier_1, lpm_2);
    EXPECT_TRUE(EpairInt_A != EpairInt_B);
    EXPECT_TRUE(EpairInt_A != EpairInt_C);
    EXPECT_TRUE(EpairInt_A != EpairInt_D);
    EXPECT_TRUE(EpairInt_A != EpairInt_E);
    EXPECT_TRUE(EpairInt_A != EpairInt_F);

    EpairIntegral Int_A(EpairInt_A);
    EpairIntegral Int_B(EpairInt_B);
    EXPECT_TRUE(Int_A != Int_B);

    InterpolationDef InterpolDef;
    EpairProductionRhoInterpolant EpairInterpol_A(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1, InterpolDef);
    EpairProductionRhoInterpolant EpairInterpol_B(tau_def, medium_1, ecuts_1, multiplier_1, lpm_1, InterpolDef);
    EpairProductionRhoInterpolant EpairInterpol_C(mu_def, medium_2, ecuts_1, multiplier_1, lpm_1, InterpolDef);
    EpairProductionRhoInterpolant EpairInterpol_D(mu_def, medium_1, ecuts_2, multiplier_1, lpm_1, InterpolDef);
    EpairProductionRhoInterpolant EpairInterpol_E(mu_def, medium_1, ecuts_1, multiplier_2, lpm_1, InterpolDef);
    EpairProductionRhoInterpolant EpairInterpol_F(mu_def, medium_1, ecuts_1, multiplier_1, lpm_2, InterpolDef);
    EXPECT_TRUE(EpairInterpol_A != EpairInterpol_B);
    EXPECT_TRUE(EpairInterpol_A != EpairInterpol_C);
    EXPECT_TRUE(EpairInterpol_A != EpairInterpol_D);
    EXPECT_TRUE(EpairInterpol_A != EpairInterpol_E);
    EXPECT_TRUE(EpairInterpol_A != EpairInterpol_F);

    EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
    EpairInterpolant Interpol_B(EpairInterpol_B, InterpolDef);
    EXPECT_TRUE(Interpol_A != Interpol_B);
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool lpm = true;

    EpairProductionRhoIntegral EpairInt_A(particle_def, medium, ecuts, multiplier, lpm);
    EpairProductionRhoIntegral EpairInt_B = EpairInt_A;
    EXPECT_TRUE(EpairInt_A == EpairInt_B);

    EpairIntegral Int_A(EpairInt_A);
    EpairIntegral Int_B = Int_A;
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    EpairProductionRhoInterpolant EpairInterpol_A(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
    EpairProductionRhoInterpolant EpairInterpol_B = EpairInterpol_A;
    EXPECT_TRUE(EpairInterpol_A == EpairInterpol_B);

    EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
    EpairInterpolant Interpol_B = Interpol_A;
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool lpm = true;

    EpairProductionRhoIntegral EpairInt_A(particle_def, medium, ecuts, multiplier, lpm);
    EpairProductionRhoIntegral EpairInt_B(EpairInt_A);
    EXPECT_TRUE(EpairInt_A == EpairInt_B);

    EpairIntegral Int_A(EpairInt_A);
    EpairIntegral Int_B(Int_A);
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    EpairProductionRhoInterpolant EpairInterpol_A(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
    EpairProductionRhoInterpolant EpairInterpol_B(EpairInterpol_A);
    EXPECT_TRUE(EpairInterpol_A == EpairInterpol_B);

    EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
    EpairInterpolant Interpol_B(Interpol_A);
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

// in polymorphism an assignmant and swap operator doesn't make sense


// TEST(Assignment , Operator ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     Medium *medium = new Air();
//     PROPOSALParticle *particle = new PROPOSALParticle(ParticleType::MuMinus,position,direction,1e5,10);
//     EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
//     Epairproduction A(particle, medium, cuts);
//     Epairproduction B(particle, medium, cuts);
//     A.SetParametrization(ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard);

//     EXPECT_TRUE(A!=B);

//     B=A;

//     EXPECT_TRUE(A==B);

//     Medium *medium2 = new Water();
//     PROPOSALParticle *particle2 = new PROPOSALParticle(ParticleType::TauMinus,position,direction,1e5,10);
//     EnergyCutSettings *cuts2 = new EnergyCutSettings(200,-1);
//     Epairproduction *C = new Epairproduction(particle2, medium2, cuts2);
//     EXPECT_TRUE(A!=*C);

//     A=*C;

//     EXPECT_TRUE(A==*C);

// }

// TEST(Assignment , Swap ) {
//     direction.SetSphericalCoordinates(1,20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     Medium *medium = new Air();
//     Medium *medium2 = new Air();
//     PROPOSALParticle *particle = new PROPOSALParticle(ParticleType::MuMinus,position,direction,1e5,10);
//     PROPOSALParticle *particle2 = new PROPOSALParticle(ParticleType::MuMinus,position,direction,1e5,10);
//     EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
//     EnergyCutSettings *cuts2 = new EnergyCutSettings(500,-1);
//     Epairproduction A(particle, medium, cuts);
//     Epairproduction B(particle2, medium2, cuts2);
//     A.EnableDEdxInterpolation();
//     B.EnableDEdxInterpolation();
//     EXPECT_TRUE(A==B);

//     Medium *medium3 = new Water();
//     Medium *medium4 = new Water();
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     PROPOSALParticle *particle3 = new PROPOSALParticle(ParticleType::TauMinus,position,direction,1e5,10);
//     PROPOSALParticle *particle4 = new PROPOSALParticle(ParticleType::TauMinus,position,direction,1e5,10);
//     EnergyCutSettings *cuts3 = new EnergyCutSettings(200,-1);
//     EnergyCutSettings *cuts4 = new EnergyCutSettings(200,-1);
//     Epairproduction *C = new Epairproduction(particle3, medium3, cuts3);
//     Epairproduction *D = new Epairproduction(particle4, medium4, cuts4);
//     EXPECT_TRUE(*C==*D);

//     A.swap(*C);

//     EXPECT_TRUE(A==*D);
//     EXPECT_TRUE(*C==B);


// }


// std::vector<Medium*>                CombOfMedium;
// std::vector<PROPOSALParticle*>      CombOfParticle;
// std::vector<EnergyCutSettings*>     CombOfEnergyCutSettings;
// std::vector<CrossSections*>         CombOfEpair;

// TEST(Epairproduction , Set_Up ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     ifstream in;
//     in.open("bin/TestFiles/Epair_dEdx.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dEdx_new;
//     double energy;
//     double dEdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;

//     cout.precision(16);

//     int i=-1;
//     double energy_old;
//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;
//         energy_old  = -1;

//         i++;
//         CombOfMedium.push_back(MediumFactory::Get()->CreateMedium(mediumName));
//         CombOfParticle.push_back(new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10));
//         CombOfParticle.at(i)->SetEnergy(energy);
//         CombOfEnergyCutSettings.push_back(new EnergyCutSettings(ecut,vcut));
//         CombOfEpair.push_back(new Epairproduction(CombOfParticle.at(i), CombOfMedium.at(i), CombOfEnergyCutSettings.at(i)));

//         CombOfEpair.at(i)->EnableLpmEffect(lpm);

//         while(energy_old < energy)
//         {
//             energy_old = energy;
//             in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;
//         }
//     }
// }

// TEST(Epairproduction , Test_of_dEdx ) {
//     ifstream in;
//     in.open("bin/TestFiles/Epair_dEdx.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dEdx_new;
//     double energy;
//     double dEdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;

//     cout.precision(16);

//     int i = -1;
//     double energy_old;
//     while(in.good())
//     {
//         in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;
//         energy_old = -1;
//         i++;
//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
//             dEdx_new=CombOfEpair.at(i)->CalculatedEdx();

//             ASSERT_NEAR(dEdx_new, dEdx, 1E-5*dEdx);

//             in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;
//         }
//     }
// }

// TEST(Epairproduction , Test_of_dNdx ) {
//     ifstream in;
//     in.open("bin/TestFiles/Epair_dNdx.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdx_new;
//     double energy;
//     double dNdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;

//     double precision = 1E-5;

//     double energy_old;
//     bool first = true;
//     int i=-1;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdx;
//         first = false;
//         energy_old = -1;

//         i++;
//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
//             dNdx_new=CombOfEpair.at(i)->CalculatedNdx();

//             if(dNdx!=0)
//             {
//                 if(log10(fabs(1-dNdx_new/dNdx))>-13)
//                 {
// //                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
// //                    cout << energy << "\t" << log10(fabs(1-dNdx_new/dNdx)) << endl;
//                 }
//             }
//             ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

//             in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdx;
//         }
//     }
// }


// TEST(Epairproduction , Test_of_dNdxrnd ) {
//     ifstream in;
//     in.open("bin/TestFiles/Epair_dNdxrnd.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdxrnd_new;
//     double energy;
//     double dNdxrnd;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;

//     double precision = 1E-5;

//     double energy_old;
//     bool first = true;
//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

//     int i=-1;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdxrnd;
//         first = false;
//         energy_old = -1;

//         i++;
//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
//             dNdxrnd_new=CombOfEpair.at(i)->CalculatedNdx(Rand->rnd());


//             if(dNdxrnd!=0){
//                 if(log10(fabs(1-dNdxrnd_new/dNdxrnd))>-14)
//                 {
// //                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
// //                    cout << energy << "\t" << log10(fabs(1-dNdxrnd_new/dNdxrnd)) << endl;
//                 }
//             }

//             ASSERT_NEAR(dNdxrnd_new, dNdxrnd, precision*dNdxrnd);

//             in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdxrnd;
//         }
//     }
//     delete Rand;
// }


// TEST(Epairproduction , Test_of_e ) {
//     ifstream in;
//     in.open("bin/TestFiles/Epair_e.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double e_new;
//     double energy;
//     double e;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;

//     double precision = 1E-5;

//     double energy_old;
//     bool first = true;
//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
//     RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
//     Rand2->rnd();
//     double rnd1,rnd2;
//     int i = -1;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
//         first = false;
//         energy_old = -1;


//         i++;
//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);

//             rnd1 = Rand->rnd();
//             rnd2 = Rand2->rnd();

//             CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
//             e_new=CombOfEpair.at(i)->CalculateStochasticLoss(rnd1,rnd2);


//             if(e!=0){
//                 if(log10(fabs(1-e_new/e))>-14)
//                 {
//                     //cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
//                     //cout << energy << "\t" << log10(fabs(1-e_new/e)) << endl;
//                 }
//             }

//             ASSERT_NEAR(e_new, e, precision*e);

//             in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
//         }
//     }
//     delete Rand;
//     delete Rand2;
// }

// TEST(Epairproduction , Test_of_dEdx_interpol ) {
//     ifstream in;
//     in.open("bin/TestFiles/Epair_dEdx_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dEdx_new;
//     double energy;
//     double dEdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;

//     double precision = 1E-5;

//     double energy_old;
//     bool first = true;
//     int i = -1;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;
//         first = false;
//         energy_old = -1;


//         i++;
//         CombOfEpair.at(i)->EnableDEdxInterpolation();
//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
//             dEdx_new=CombOfEpair.at(i)->CalculatedEdx();

//             ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

//             in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;
//         }
//     }
// }

// TEST(Epairproduction , Test_of_dNdx_interpol ) {
//     ifstream in;
//     in.open("bin/TestFiles/Epair_dNdx_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdx_new;
//     double energy;
//     double dNdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;

//     double precision = 1E-2;

//     double energy_old;
//     bool first = true;

//     int i=-1;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdx;
//         first = false;
//         energy_old = -1;

//         i++;
//         CombOfEpair.at(i)->EnableDNdxInterpolation();

//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
//             dNdx_new=CombOfEpair.at(i)->CalculatedNdx();

//             if(dNdx!=0)
//                 if(log10(fabs(1-dNdx_new/dNdx))>-8)
//                 {
//                     //cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
//                     //cout << energy << "\t" << log10(fabs(1-dNdx_new/dNdx)) << endl;
//                 }

//             ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

//             in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdx;
//         }
//     }
// }

// TEST(Epairproduction , Test_of_dNdxrnd_interpol ) {
//     ifstream in;
//     in.open("bin/TestFiles/Epair_dNdxrnd_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdxrnd_new;
//     double energy;
//     double dNdxrnd;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;

//     double precision = 1E-2;

//     double energy_old;
//     bool first = true;
//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

//     int i=-1;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdxrnd;
//         first = false;
//         energy_old = -1;

//         i++;
//         CombOfEpair.at(i)->EnableDNdxInterpolation();
//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
//             dNdxrnd_new=CombOfEpair.at(i)->CalculatedNdx(Rand->rnd());


//             if(dNdxrnd!=0){
//                 if(log10(fabs(1-dNdxrnd_new/dNdxrnd))>-5)
//                 {
//                     log_warn("%f \t %f \t %i \t %f \t %s \t %s \t",ecut,vcut,lpm,energy,mediumName.c_str(),particleName.c_str());
//                     log_warn("%f \t %f",energy,log10(fabs(1-dNdxrnd_new/dNdxrnd)));
//                 }
//             }

//             ASSERT_NEAR(dNdxrnd_new, dNdxrnd, precision*dNdxrnd);

//             in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdxrnd;
//         }
//     }
//     delete Rand;
// }

// TEST(Epairproduction , Test_of_e_interpol ) {
//     ifstream in;
//     in.open("bin/TestFiles/Epair_e_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double e_new;
//     double energy;
//     double e;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;

//     double precision = 1E-2;
//     double precision_old=precision;
//     double energy_old;
//     bool first = true;
//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
//     RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
//     Rand2->rnd();
//     double rnd1,rnd2;
//     int i=-1;
//     int ctr=0;
//     while(in.good())
//     {
//         if(first)in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
//         first = false;
//         energy_old = -1;

//         i++;
//         CombOfEpair.at(i)->EnableDEdxInterpolation();
//         CombOfEpair.at(i)->EnableDNdxInterpolation();

//         while(energy_old < energy)
//         {
//             precision=precision_old;
//             energy_old = energy;

//             rnd1 = Rand->rnd();
//             rnd2 = Rand2->rnd();

//             //The Cross section for such high energy is just extrapolated and
//             //therefore pretty unceartain.
//             if(!particleName.compare("e") && energy > 1E12)precision = 1E-1;

//             CombOfEpair.at(i)->GetParticle()->SetEnergy(energy);
//             e_new=CombOfEpair.at(i)->CalculateStochasticLoss(rnd1,rnd2);


//             if(e!=0){
//                 if(log10(fabs(1-e_new/e))>-3)
//                 {
//                     //cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
//                     cout << energy << "\t" << log10(fabs(1-e_new/e)) << endl;
//                 }
//             }
//             //cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
//             //cout << "ctr: " << ctr++ << endl;
//             EXPECT_NEAR(e_new, e, precision*e);

//             in>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
//         }
//     }
//     delete Rand;
//     delete Rand2;
// }



int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
