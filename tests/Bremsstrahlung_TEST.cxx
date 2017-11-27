
// #include <iostream>
// #include <string>

#include "gtest/gtest.h"

#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/BremsIntegral.h"
#include "PROPOSAL/crossection/BremsInterpolant.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"

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

    BremsKelnerKokoulinPetrukhin* Brems_A = new BremsKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, lpm);
    Parametrization* Brems_B = new BremsKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, lpm);
    EXPECT_TRUE(*Brems_A == *Brems_B);

    BremsKelnerKokoulinPetrukhin param(particle_def, medium, ecuts, multiplier, lpm);
    EXPECT_TRUE(param == *Brems_A);

    BremsIntegral* Int_A = new BremsIntegral(param);
    CrossSectionIntegral* Int_B = new BremsIntegral(param);
    EXPECT_TRUE(*Int_A == *Int_B);

    InterpolationDef InterpolDef;
    BremsInterpolant* Interpol_A = new BremsInterpolant(param, InterpolDef);
    CrossSectionInterpolant* Interpol_B = new BremsInterpolant(param, InterpolDef);
    EXPECT_TRUE(*Interpol_A == *Interpol_B);

    delete Brems_A;
    delete Brems_B;
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
    bool lpm_1 = true;
    bool lpm_2 = false;

    BremsKelnerKokoulinPetrukhin Brems_A(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    BremsKelnerKokoulinPetrukhin Brems_B(tau_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    BremsKelnerKokoulinPetrukhin Brems_C(mu_def, medium_2, ecuts_1, multiplier_1, lpm_1);
    BremsKelnerKokoulinPetrukhin Brems_D(mu_def, medium_1, ecuts_2, multiplier_1, lpm_1);
    BremsKelnerKokoulinPetrukhin Brems_E(mu_def, medium_1, ecuts_1, multiplier_2, lpm_1);
    BremsKelnerKokoulinPetrukhin Brems_F(mu_def, medium_1, ecuts_1, multiplier_1, lpm_2);
    EXPECT_TRUE(Brems_A != Brems_B);
    EXPECT_TRUE(Brems_A != Brems_C);
    EXPECT_TRUE(Brems_A != Brems_C);
    EXPECT_TRUE(Brems_A != Brems_E);
    EXPECT_TRUE(Brems_A != Brems_F);

    BremsAndreevBezrukovBugaev param_2(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    BremsPetrukhinShestakov param_3(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    BremsCompleteScreening param_4(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
    EXPECT_TRUE(Brems_A != param_2);
    EXPECT_TRUE(Brems_A != param_3);
    EXPECT_TRUE(Brems_A != param_4);
    EXPECT_TRUE(param_2 != param_3);
    EXPECT_TRUE(param_2 != param_4);
    EXPECT_TRUE(param_3 != param_4);

    BremsIntegral Int_A(param_2);
    BremsIntegral Int_B(param_3);
    EXPECT_TRUE(Int_A != Int_B);

    InterpolationDef InterpolDef;
    BremsInterpolant Interpol_A(param_2, InterpolDef);
    BremsInterpolant Interpol_B(param_3, InterpolDef);
    EXPECT_TRUE(Interpol_A != Interpol_B);
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool lpm = true;

    BremsKelnerKokoulinPetrukhin Brems_A(particle_def, medium, ecuts, multiplier, lpm);
    BremsKelnerKokoulinPetrukhin Brems_B = Brems_A;
    EXPECT_TRUE(Brems_A == Brems_B);

    BremsIntegral Int_A(Brems_A);
    BremsIntegral Int_B = Int_A;
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    BremsInterpolant Interpol_A(Brems_A, InterpolDef);
    BremsInterpolant Interpol_B = Interpol_A;
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool lpm = true;

    BremsKelnerKokoulinPetrukhin Brems_A(particle_def, medium, ecuts, multiplier, lpm);
    BremsKelnerKokoulinPetrukhin Brems_B(Brems_A);
    EXPECT_TRUE(Brems_A == Brems_B);

    BremsIntegral Int_A(Brems_A);
    BremsIntegral Int_B(Int_A);
    EXPECT_TRUE(Int_A == Int_B);

    InterpolationDef InterpolDef;
    BremsInterpolant Interpol_A(Brems_A, InterpolDef);
    BremsInterpolant Interpol_B(Interpol_A);
    EXPECT_TRUE(Interpol_A == Interpol_B);
}

// in polymorphism an assignmant and swap operator doesn't make sense


// TEST(Bremsstrahlung , Test_of_dEdx ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Brems_dEdx.txt");

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
//     int para;

//     cout.precision(16);


//     while(in.good())
//     {
//         in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;

//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//         CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);


//         brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
//         brems->EnableLpmEffect(lpm);



//         dEdx_new=brems->CalculatedEdx();

//         ASSERT_NEAR(dEdx_new, dEdx, 1e-7*dEdx);

//         delete cuts;
//         delete medium;
//         delete particle;
//         delete brems;
//     }
// }

// TEST(Bremsstrahlung , Test_of_dNdx ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Brems_dNdx.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdx;
//     double dNdx_new;
//     double energy;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;
//     int para;

//     cout.precision(16);


//     while(in.good())
//     {
//         in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdx;

//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//         CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);


//         brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
//         brems->EnableLpmEffect(lpm);

//         dNdx_new=brems->CalculatedNdx();
//         ASSERT_NEAR(dNdx_new, dNdx, 1e-7*dNdx);

//         delete cuts;
//         delete medium;
//         delete particle;
//         delete brems;



//     }
// }

// TEST(Bremsstrahlung , Test_of_dNdxrnd ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Brems_dNdxrnd.txt");

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
//         if(first)in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdxrnd;
//         first=false;
//         energy_old = -1;
//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//         CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
//         brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
//         brems->EnableLpmEffect(lpm);

//         //cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << mediumName << "\t" << particleName<< "\t" << dNdxrnd << endl;

//         while(energy_old < energy){
//             energy_old = energy;
//             brems->GetParticle()->SetEnergy(energy);
//             dNdxrnd_new=brems->CalculatedNdx(Rand->rnd());

//             ASSERT_NEAR(dNdxrnd_new, dNdxrnd, 1E-7*dNdxrnd);

//             in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdxrnd;
//         }



//         delete cuts;
//         delete medium;
//         delete particle;
//         delete brems;
//     }
//     delete Rand;
// }



// TEST(Bremsstrahlung , Test_of_e ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Brems_e.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);

//     double e;
//     double e_new;
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
//     RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
//     Rand2->rnd();

//     double rnd1, rnd2;
//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
//         first=false;
//         energy_old = -1;
//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//         CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
//         brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
//         brems->EnableLpmEffect(lpm);

//         //cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << mediumName << "\t" << particleName<< "\t" << e << endl;

//         while(energy_old < energy){
//             energy_old = energy;
//             brems->GetParticle()->SetEnergy(energy);

//             rnd1 = Rand->rnd();
//             rnd2 = Rand2->rnd();

//             e_new=brems->CalculateStochasticLoss(rnd1,rnd2);
//             ASSERT_NEAR(e_new, e, 1E-7*e);

//             in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
//         }

//         delete cuts;
//         delete medium;
//         delete particle;
//         delete brems;
//     }
//     delete Rand2;
//     delete Rand;
// }


// TEST(Bremsstrahlung , Test_of_dEdx_Interpolant ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();

//     ifstream in;
//     in.open("bin/TestFiles/Brems_dEdx_interpol.txt");
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
//     int para;

//     cout.precision(16);
//     double precision;
//     double precisionOld = 1E-2;
//     bool first=true;
//     double energy_old;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;
//         first=false;

//         precision = precisionOld;
//         energy_old =-1;

//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
//         CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);

//         brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
//         brems->EnableLpmEffect(lpm);
//         brems->EnableDEdxInterpolation();

//         while(energy_old < energy){
//             energy_old = energy;
//             brems->GetParticle()->SetEnergy(energy);
//             dEdx_new=brems->CalculatedEdx();

//             if(!particleName.compare("tau") && energy < 10001)precision = 0.5;
//             if(!particleName.compare("e") && energy > 1E10)precision = 0.5;

//             ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

//             in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;

//             precision = precisionOld;

//         }



//         delete cuts;
//         delete medium;
//         delete particle;
//         delete brems;



//     }
// }

// TEST(Bremsstrahlung , Test_of_dNdx_Interpolant ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     ifstream in;
//     in.open("bin/TestFiles/Brems_dNdx_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);

//     double dNdx;
//     double dNdx_new;
//     double energy;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;
//     int para;

//     cout.precision(16);
//     double energy_old=-1;

//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdx;
//         first=false;
//         energy_old = -1;
//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//         CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
//         brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
//         brems->EnableLpmEffect(lpm);
//         brems->EnableDNdxInterpolation();

//         while(energy_old < energy){
//             energy_old = energy;
//             brems->GetParticle()->SetEnergy(energy);
//             dNdx_new=brems->CalculatedNdx();

//             ASSERT_NEAR(dNdx_new, dNdx, 1E-6*dNdx);

//             in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdx;
//         }



//         delete cuts;
//         delete medium;
//         delete particle;
//         delete brems;
//     }
// }

// TEST(Bremsstrahlung , Test_of_e_interpol ) {
// return;
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     ifstream in;
//     in.open("bin/TestFiles/Brems_e_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);

//     double e;
//     double e_new;
//     double energy;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     bool lpm;
//     int para;

//     cout.precision(16);
//     double energy_old=-1;
//     double precision = 1E-5;
//     double precision_old = precision;
//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
//     RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
//     Rand2->rnd();

//     double rnd1,rnd2;
//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
//         first=false;
//         energy_old = -1;
//         Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
//         PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
//         particle->SetEnergy(energy);
//         EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

//         CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
//         brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
//         brems->EnableLpmEffect(lpm);
//         brems->EnableDNdxInterpolation();


//         while(energy_old < energy){
//             energy_old = energy;
//             brems->GetParticle()->SetEnergy(energy);
//             rnd1 = Rand->rnd();
//             rnd2 = Rand2->rnd();

//             e_new = brems->CalculateStochasticLoss(rnd1,rnd2);

//             ASSERT_NEAR(e_new, e, 1*e);

//             in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
//             precision = precision_old;
//         }



//         delete cuts;
//         delete medium;
//         delete particle;
//         delete brems;
//     }
//     delete Rand2;
//     delete Rand;
// }

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
