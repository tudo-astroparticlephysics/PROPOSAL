
// #include <iostream>
// #include <string>

// #include <math.h>

#include "gtest/gtest.h"

#include "PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crossection/PhotoIntegral.h"
#include "PROPOSAL/crossection/PhotoInterpolant.h"
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
    bool hardbb = true;
    ShadowButkevichMikhailov shadow;
    InterpolationDef InterpolDef;

    PhotoKokoulin* PhotoReal_A = new PhotoKokoulin(particle_def, medium, ecuts, hardbb, multiplier);
    Parametrization* PhotoReal_B = new PhotoKokoulin(particle_def, medium, ecuts, hardbb, multiplier);
    EXPECT_TRUE(*PhotoReal_A == *PhotoReal_B);

    PhotoKokoulin param_PhotoReal(particle_def, medium, ecuts, hardbb, multiplier);
    EXPECT_TRUE(param_PhotoReal == *PhotoReal_A);

    PhotoIntegral* Int_PhotoReal_A = new PhotoIntegral(param_PhotoReal);
    CrossSectionIntegral* Int_PhotoReal_B = new PhotoIntegral(param_PhotoReal);
    EXPECT_TRUE(*Int_PhotoReal_A == *Int_PhotoReal_B);

    PhotoInterpolant* Interpol_PhotoReal_A = new PhotoInterpolant(param_PhotoReal, InterpolDef);
    CrossSectionInterpolant* Interpol_PhotoReal_B = new PhotoInterpolant(param_PhotoReal, InterpolDef);
    EXPECT_TRUE(*Interpol_PhotoReal_A == *Interpol_PhotoReal_B);

    delete PhotoReal_A;
    delete PhotoReal_B;
    delete Int_PhotoReal_A;
    delete Int_PhotoReal_B;
    delete Interpol_PhotoReal_A;
    delete Interpol_PhotoReal_B;

    PhotoAbramowiczLevinLevyMaor97* PhotoQ2_A = new PhotoAbramowiczLevinLevyMaor97(particle_def, medium, ecuts, shadow, multiplier);
    Parametrization* PhotoQ2_B = new PhotoAbramowiczLevinLevyMaor97(particle_def, medium, ecuts, shadow, multiplier);
    EXPECT_TRUE(*PhotoQ2_A == *PhotoQ2_B);

    PhotoAbramowiczLevinLevyMaor97 param_Q2(particle_def, medium, ecuts, shadow, multiplier);
    EXPECT_TRUE(param_Q2 == *PhotoQ2_A);

    PhotoIntegral* Int_PhotoQ2_A = new PhotoIntegral(param_Q2);
    CrossSectionIntegral* Int_PhotoQ2_B = new PhotoIntegral(param_Q2);
    EXPECT_TRUE(*Int_PhotoQ2_A == *Int_PhotoQ2_B);

    PhotoInterpolant* Interpol_PhotoQ2_A = new PhotoInterpolant(param_Q2, InterpolDef);
    CrossSectionInterpolant* Interpol_PhotoQ2_B = new PhotoInterpolant(param_Q2, InterpolDef);
    EXPECT_TRUE(*Interpol_PhotoQ2_A == *Interpol_PhotoQ2_B);

    delete PhotoQ2_A;
    delete PhotoQ2_B;
    delete Int_PhotoQ2_A;
    delete Int_PhotoQ2_B;
    delete Interpol_PhotoQ2_A;
    delete Interpol_PhotoQ2_B;
}

TEST(Comparison, Comparison_not_equal)
{
    ParticleDef mu_def = MuMinusDef::Get();
    ParticleDef tau_def = TauMinusDef::Get();
    Water medium_1;
    Ice medium_2;
    EnergyCutSettings ecuts_1(500, 0.05);
    EnergyCutSettings ecuts_2(-1, 0.05);
    double multiplier_1 = 1.;
    double multiplier_2 = 2.;
    bool hardbb = true;
    bool softbb = false;
    ShadowButkevichMikhailov shadow_1;
    ShadowDuttaRenoSarcevicSeckel shadow_2;
    InterpolationDef InterpolDef;

    PhotoKokoulin PhotoReal_A(mu_def, medium_1, ecuts_1, hardbb, multiplier_1);
    // PhotoKokoulin PhotoReal_B(tau_def, medium_1, ecuts_1, hardbb, multiplier_1);
    // PhotoKokoulin PhotoReal_C(mu_def, medium_2, ecuts_1, hardbb, multiplier_1);
    // PhotoKokoulin PhotoReal_D(mu_def, medium_1, ecuts_2, hardbb, multiplier_1);
    PhotoKokoulin PhotoReal_E(mu_def, medium_1, ecuts_1, softbb, multiplier_1);
    // PhotoKokoulin PhotoReal_F(mu_def, medium_1, ecuts_1, hardbb, multiplier_2);
    // EXPECT_TRUE(PhotoReal_A != PhotoReal_B);
    // EXPECT_TRUE(PhotoReal_A != PhotoReal_C);
    // EXPECT_TRUE(PhotoReal_A != PhotoReal_D);
    EXPECT_TRUE(PhotoReal_A != PhotoReal_E);
    // EXPECT_TRUE(PhotoReal_A != PhotoReal_F);

    // PhotoZeus param_Real_2(mu_def, medium_1, ecuts_1, hardbb, multiplier_1);
    // PhotoBezrukovBugaev param_Real_3(mu_def, medium_1, ecuts_1, hardbb, multiplier_1);
    // PhotoRhode param_Real_4(mu_def, medium_1, ecuts_1, hardbb, multiplier_1);
    // EXPECT_TRUE(PhotoReal_A != param_Real_2);
    // EXPECT_TRUE(PhotoReal_A != param_Real_3);
    // EXPECT_TRUE(PhotoReal_A != param_Real_4);
    // EXPECT_TRUE(param_Real_2 != param_Real_3);
    // EXPECT_TRUE(param_Real_2 != param_Real_4);
    // EXPECT_TRUE(param_Real_3 != param_Real_4);
    //
    // PhotoIntegral Int_PhotoReal_A(PhotoReal_A);
    // PhotoIntegral Int_PhotoReal_B(PhotoReal_B);
    // EXPECT_TRUE(Int_PhotoReal_A != Int_PhotoReal_B);
    //
    // PhotoInterpolant Interpol_PhotoReal_A(PhotoReal_A, InterpolDef);
    // PhotoInterpolant Interpol_PhotoReal_B(PhotoReal_B, InterpolDef);
    // EXPECT_TRUE(Interpol_PhotoReal_A != Interpol_PhotoReal_B);
    //
    PhotoAbramowiczLevinLevyMaor97 PhotoQ2_A(mu_def, medium_1, ecuts_1, shadow_1, multiplier_1);
    // PhotoAbramowiczLevinLevyMaor97 PhotoQ2_B(tau_def, medium_1, ecuts_1, shadow_1, multiplier_1);
    // PhotoAbramowiczLevinLevyMaor97 PhotoQ2_C(mu_def, medium_2, ecuts_1, shadow_1, multiplier_1);
    // PhotoAbramowiczLevinLevyMaor97 PhotoQ2_D(mu_def, medium_1, ecuts_2, shadow_1, multiplier_1);
    PhotoAbramowiczLevinLevyMaor97 PhotoQ2_E(mu_def, medium_1, ecuts_1, shadow_2, multiplier_1);
    // PhotoAbramowiczLevinLevyMaor97 PhotoQ2_F(mu_def, medium_1, ecuts_1, shadow_1, multiplier_2);
    // EXPECT_TRUE(PhotoQ2_A != PhotoQ2_B);
    // EXPECT_TRUE(PhotoQ2_A != PhotoQ2_C);
    // EXPECT_TRUE(PhotoQ2_A != PhotoQ2_D);
    EXPECT_TRUE(PhotoQ2_A != PhotoQ2_E);
    // EXPECT_TRUE(PhotoQ2_A != PhotoQ2_F);
    //
    // EXPECT_TRUE(PhotoReal_A != PhotoQ2_A);
    //
    // PhotoAbramowiczLevinLevyMaor91 param_Q2_2(mu_def, medium_1, ecuts_1, shadow_1, multiplier_1);
    // PhotoButkevichMikhailov param_Q2_3(mu_def, medium_1, ecuts_1, shadow_1, multiplier_1);
    // PhotoRenoSarcevicSu param_Q2_4(mu_def, medium_1, ecuts_1, shadow_1, multiplier_1);
    // EXPECT_TRUE(PhotoQ2_A != param_Q2_2);
    // EXPECT_TRUE(PhotoQ2_A != param_Q2_3);
    // EXPECT_TRUE(PhotoQ2_A != param_Q2_4);
    // EXPECT_TRUE(param_Q2_2 != param_Q2_3);
    // EXPECT_TRUE(param_Q2_2 != param_Q2_4);
    // EXPECT_TRUE(param_Q2_3 != param_Q2_4);
    //
    // PhotoIntegral Int_PhotoQ2_A(PhotoQ2_A);
    // PhotoIntegral Int_PhotoQ2_B(PhotoQ2_B);
    // EXPECT_TRUE(Int_PhotoQ2_A != Int_PhotoQ2_B);
    //
    // PhotoInterpolant Interpol_PhotoQ2_A(PhotoQ2_A, InterpolDef);
    // PhotoInterpolant Interpol_PhotoQ2_B(PhotoQ2_B, InterpolDef);
    // EXPECT_TRUE(Interpol_PhotoQ2_A != Interpol_PhotoQ2_B);
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool hardbb = true;
    ShadowButkevichMikhailov shadow;
    InterpolationDef InterpolDef;

    PhotoKokoulin PhotoReal_A(particle_def, medium, ecuts, hardbb, multiplier);
    PhotoKokoulin PhotoReal_B = PhotoReal_A;
    EXPECT_TRUE(PhotoReal_A == PhotoReal_B);

    PhotoIntegral Int_PhotoReal_A(PhotoReal_A);
    PhotoIntegral Int_PhotoReal_B = Int_PhotoReal_A;
    EXPECT_TRUE(Int_PhotoReal_A == Int_PhotoReal_B);

    PhotoInterpolant Interpol_PhotoReal_A(PhotoReal_A, InterpolDef);
    PhotoInterpolant Interpol_PhotoReal_B = Interpol_PhotoReal_A;
    EXPECT_TRUE(Interpol_PhotoReal_A == Interpol_PhotoReal_B);

    PhotoAbramowiczLevinLevyMaor97 PhotoQ2_A(particle_def, medium, ecuts, shadow, multiplier);
    PhotoAbramowiczLevinLevyMaor97 PhotoQ2_B = PhotoQ2_A;
    EXPECT_TRUE(PhotoQ2_A == PhotoQ2_B);

    PhotoIntegral Int_PhotoQ2_A(PhotoQ2_A);
    PhotoIntegral Int_PhotoQ2_B = Int_PhotoQ2_A;
    EXPECT_TRUE(Int_PhotoQ2_A == Int_PhotoQ2_B);

    PhotoInterpolant Interpol_PhotoQ2_A(PhotoQ2_A, InterpolDef);
    PhotoInterpolant Interpol_PhotoQ2_B = Interpol_PhotoQ2_A;
    EXPECT_TRUE(Interpol_PhotoQ2_A == Interpol_PhotoQ2_B);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef particle_def = MuMinusDef::Get();
    Water medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    bool hardbb = true;
    ShadowButkevichMikhailov shadow;
    InterpolationDef InterpolDef;

    PhotoKokoulin PhotoReal_A(particle_def, medium, ecuts, hardbb, multiplier);
    PhotoKokoulin PhotoReal_B(PhotoReal_A);
    EXPECT_TRUE(PhotoReal_A == PhotoReal_B);

    PhotoIntegral Int_PhotoReal_A(PhotoReal_A);
    PhotoIntegral Int_PhotoReal_B(Int_PhotoReal_A);
    EXPECT_TRUE(Int_PhotoReal_A == Int_PhotoReal_B);

    PhotoInterpolant Interpol_PhotoReal_A(PhotoReal_A, InterpolDef);
    PhotoInterpolant Interpol_PhotoReal_B(Interpol_PhotoReal_A);
    EXPECT_TRUE(Interpol_PhotoReal_A == Interpol_PhotoReal_B);

    PhotoAbramowiczLevinLevyMaor97 PhotoQ2_A(particle_def, medium, ecuts, shadow, multiplier);
    PhotoAbramowiczLevinLevyMaor97 PhotoQ2_B(PhotoQ2_A);
    EXPECT_TRUE(PhotoQ2_A == PhotoQ2_B);

    PhotoIntegral Int_PhotoQ2_A(PhotoQ2_A);
    PhotoIntegral Int_PhotoQ2_B(Int_PhotoQ2_A);
    EXPECT_TRUE(Int_PhotoQ2_A == Int_PhotoQ2_B);

    PhotoInterpolant Interpol_PhotoQ2_A(PhotoQ2_A, InterpolDef);
    PhotoInterpolant Interpol_PhotoQ2_B(Interpol_PhotoQ2_A);
    EXPECT_TRUE(Interpol_PhotoQ2_A == Interpol_PhotoQ2_B);
}

// in polymorphism an assignmant and swap operator doesn't make sense


// TEST(Assignment , Operator ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     Medium *medium = new Air();
//     PROPOSALParticle *particle = new PROPOSALParticle(ParticleType::MuMinus,position,direction,1e5,10);
//     EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
//     Photonuclear A(particle, medium, cuts);
//     Photonuclear B(particle, medium, cuts);
//     A.SetParametrization(ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard);
//     EXPECT_TRUE(A!=B);

//     B=A;

//     EXPECT_TRUE(A==B);

//     Medium *medium2 = new Water();
//     PROPOSALParticle *particle2 = new PROPOSALParticle(ParticleType::TauMinus,position,direction,1e5,10);
//     EnergyCutSettings *cuts2 = new EnergyCutSettings(200,-1);
//     Photonuclear *C = new Photonuclear(particle2, medium2, cuts2);
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
//     Photonuclear A(particle, medium, cuts);
//     Photonuclear B(particle2, medium2, cuts2);
//     A.SetParametrization(ParametrizationType::PhotoKokoulinShadowBezrukovHard);
//     B.SetParametrization(ParametrizationType::PhotoKokoulinShadowBezrukovHard);
//     A.EnableDEdxInterpolation();
//     B.EnableDEdxInterpolation();
//     A.EnableDNdxInterpolation();
//     B.EnableDNdxInterpolation();
//     EXPECT_TRUE(A==B);

//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     Medium *medium3 = new Water();
//     Medium *medium4 = new Water();
//     PROPOSALParticle *particle3 = new PROPOSALParticle(ParticleType::TauMinus,position,direction,1e5,10);
//     PROPOSALParticle *particle4 = new PROPOSALParticle(ParticleType::TauMinus,position,direction,1e5,10);
//     EnergyCutSettings *cuts3 = new EnergyCutSettings(200,-1);
//     EnergyCutSettings *cuts4 = new EnergyCutSettings(200,-1);
//     Photonuclear *C = new Photonuclear(particle3, medium3, cuts3);
//     Photonuclear *D = new Photonuclear(particle4, medium4, cuts4);
//     EXPECT_TRUE(*C==*D);

//     A.swap(*C);

//     EXPECT_TRUE(A==*D);
//     EXPECT_TRUE(*C==B);


// }
// //----------------------------------------------------------------------------//

// std::vector<Medium*>                CombOfMedium;
// std::vector<PROPOSALParticle*>      CombOfParticle;
// std::vector<EnergyCutSettings*>     CombOfEnergyCutSettings;
// std::vector<Photonuclear*>          CombOfPhoto;

// //----------------------------------------------------------------------------//

// TEST(Photonuclear , Set_Up ) {
//     direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
//     direction.CalculateCartesianFromSpherical();
//     ifstream in;
//     in.open("bin/TestFiles/Photo_dEdx.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double energy;
//     double dEdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     int para;

//     cout.precision(16);

//     int i=-1;
//     double energy_old;
//     bool first = true;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dEdx;
//         energy_old  = -1;

//         i++;
//         CombOfMedium.push_back(MediumFactory::Get()->CreateMedium(mediumName));
//         CombOfParticle.push_back(new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10));
//         CombOfParticle.at(i)->SetEnergy(energy);
//         CombOfEnergyCutSettings.push_back(new EnergyCutSettings(ecut,vcut));
//         CombOfPhoto.push_back(new Photonuclear(CombOfParticle.at(i), CombOfMedium.at(i), CombOfEnergyCutSettings.at(i)));

//         // Now: parametrization_ =  31,  Former: form=1 and bb=1 Kokoulin
//         // Now: parametrization_ = -31,  Former: form=2 and bb=1 Kokoulin + hard component
//         // Now: parametrization_ =  32,  Former: form=1 and bb=2 Rhode
//         // Now: parametrization_ = -32,  Former: form=2 and bb=2 Rhode + hard component
//         // Now: parametrization_ =  33,  Former: form=1 and bb=3 Bezrukov/Bugaev
//         // Now: parametrization_ = -33,  Former: form=2 and bb=3 Bezrukov/Bugaev + hard component
//         // Now: parametrization_ =  34,  Former: form=1 and bb=4 Zeus
//         // Now: parametrization_ = -34,  Former: form=2 and bb=4 Zeus + hard component
//         // Now: parametrization_ =  35,  Former: form=3 and bb=1 shadow=1 ALLM 91
//         // Now: parametrization_ = -35, Former: form=3 and bb=1 shadow=2 ALLM 91
//         // Now: parametrization_ =  36, Former: form=3 and bb=2 shadow=1 ALLM 97
//         // Now: parametrization_ = -36, Former: form=3 and bb=2 shadow=2 ALLM 97
//         // Now: parametrization_ =  37, Former: form=4 and bb=1 shadow=1 Butkevich/Mikhailov
//         // Now: parametrization_ = -37, Former: form=4 and bb=1 shadow=2 Butkevich/Mikhailov
//         // Now: parametrization_ =  38, New parametrization for supersymmetric particles
//         // Now: parametrization_ = -38, New parametrization for supersymmetric particles

//         if(para ==  31)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoKokoulinShadowBezrukovSoft);
//         if(para == -31)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoKokoulinShadowBezrukovHard);
//         if(para ==  32)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoRhodeShadowBezrukovSoft);
//         if(para == -32)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoRhodeShadowBezrukovHard);
//         if(para ==  33)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft);
//         if(para == -33)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard);
//         if(para ==  34)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoZeusShadowBezrukovSoft);
//         if(para == -34)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoZeusShadowBezrukovHard);
//         if(para ==  35)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta);
//         if(para == -35)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich);
//         if(para ==  36)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta);
//         if(para == -36)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich);
//         if(para ==  37)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoButkevichMikhailovShadowDutta);
//         if(para == -37)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoButkevichMikhailovShadowButkevich);
//         if(para ==  38)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoRenoSarcevicSuShadowDutta);
//         if(para == -38)CombOfPhoto.at(i)->SetParametrization(ParametrizationType::PhotoRenoSarcevicSuShadowButkevich);

//         while(energy_old < energy)
//         {
//             energy_old = energy;
//             in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dEdx;
//         }
//     }
// }

// //----------------------------------------------------------------------------//

// TEST(Photonuclear , Test_of_e ) {
//     ifstream in;
//     in.open("bin/TestFiles/Photo_e.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double e_new;
//     double energy;
//     double e;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     int para;
//     double precision = 1E-6;

//     double energy_old;
//     bool first = true;
//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
//     RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
//     Rand2->rnd();
//     double rnd1,rnd2;
//     int i = -1;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>e;
//         first = false;
//         energy_old = -1;


//         i++;

//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);

//             rnd1 = Rand->rnd();
//             rnd2 = Rand2->rnd();

//             CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
//             e_new=CombOfPhoto.at(i)->CalculateStochasticLoss(rnd1,rnd2);


//             if(e!=0){
//                 if(log10(fabs(1-e_new/e))>-8)
//                 {
// //                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
// //                    cout << energy << "\t" << log10(fabs(1-e_new/e)) << endl;
//                 }
//             }

//             ASSERT_NEAR(e_new, e, precision*e);

//             in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>e;
//         }
//     }
//     delete Rand;
//     delete Rand2;
// }

// //----------------------------------------------------------------------------//

// TEST(Photonuclear , Test_of_dNdxrnd ) {
//     ifstream in;
//     in.open("bin/TestFiles/Photo_dNdxrnd.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdxrnd_new;
//     double energy;
//     double dNdxrnd;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     int para;
//     double precision = 1E-6;

//     double energy_old;
//     bool first = true;
//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

//     int i=-1;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdxrnd;
//         first = false;
//         energy_old = -1;

//         i++;

//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
//             dNdxrnd_new=CombOfPhoto.at(i)->CalculatedNdx(Rand->rnd());


//             if(dNdxrnd!=0){
//                 if(log10(fabs(1-dNdxrnd_new/dNdxrnd))>-14)
//                 {
// //                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
// //                    cout << energy << "\t" << log10(fabs(1-dNdxrnd_new/dNdxrnd)) << endl;
//                 }
//             }

//             ASSERT_NEAR(dNdxrnd_new, dNdxrnd, precision*dNdxrnd);

//             in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdxrnd;
//         }
//     }
//     delete Rand;
// }

// //----------------------------------------------------------------------------//


// TEST(Photonuclear , Test_of_dNdx ) {
//     ifstream in;
//     in.open("bin/TestFiles/Photo_dNdx.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdx_new;
//     double energy;
//     double dNdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     int para;

//     double precision = 1E-6;

//     double energy_old;
//     bool first = true;
//     int i=-1;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdx;
//         first = false;
//         energy_old = -1;

//         i++;

//         while(energy_old < energy)
//         {
//             energy_old = energy;
//             CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
//             dNdx_new=CombOfPhoto.at(i)->CalculatedNdx();

//             if(dNdx!=0)
//             {
//                 if(log10(fabs(1-dNdx_new/dNdx))>-13)
//                 {
// //                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
// //                    cout << energy << "\t" << log10(fabs(1-dNdx_new/dNdx)) << endl;
//                 }
//             }
//             ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

//             in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdx;
//         }
//     }
// }

// //----------------------------------------------------------------------------//


// TEST(Photonuclear , Test_of_dEdx ) {
//     ifstream in;
//     in.open("bin/TestFiles/Photo_dEdx.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double energy;
//     double dEdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     int para;
//     double dEdx_new;

//     cout.precision(16);

//     int i = -1;
//     double energy_old;
//     while(in.good())
//     {
//         in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dEdx;
//         energy_old = -1;
//         i++;
//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
//             dEdx_new=CombOfPhoto.at(i)->CalculatedEdx();

//             ASSERT_NEAR(dEdx_new, dEdx, 1E-6*dEdx);

//             in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dEdx;
//         }
//     }
// }

// //----------------------------------------------------------------------------//


// TEST(Photonuclear , Test_of_dEdx_interpol ) {
//     ifstream in;
//     in.open("bin/TestFiles/Photo_dEdx_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dEdx_new;
//     double energy;
//     double dEdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     int para;
//     double precision = 1E-5;

//     double energy_old;
//     bool first = true;
//     int i = -1;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dEdx;
//         first = false;
//         energy_old = -1;


//         i++;

//         CombOfPhoto.at(i)->EnableDEdxInterpolation();
//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
//             dEdx_new=CombOfPhoto.at(i)->CalculatedEdx();

//             ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

//             in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dEdx;
//         }
//     }
// }

// //----------------------------------------------------------------------------//

// TEST(Photonuclear , Test_of_dNdx_interpol ) {
//     //return;
//     ifstream in;
//     in.open("bin/TestFiles/Photo_dNdx_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdx_new;
//     double energy;
//     double dNdx;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     int para;
//     double precision = 1E-2;

//     double energy_old;
//     bool first = true;

//     int i=-1;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdx;
//         first = false;
//         energy_old = -1;

//         i++;

//         CombOfPhoto.at(i)->EnableDNdxInterpolation();

//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
//             dNdx_new=CombOfPhoto.at(i)->CalculatedNdx();

//             if(dNdx!=0)
//                 if(log10(fabs(1-dNdx_new/dNdx))>-8)
//                 {
// //                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
// //                    cout << energy << "\t" << log10(fabs(1-dNdx_new/dNdx)) << endl;
//                 }

//             ASSERT_NEAR(dNdx_new, dNdx, precision*dNdx);

//             in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdx;
//         }
//     }
// }

// //----------------------------------------------------------------------------//

// TEST(Photonuclear , Test_of_dNdxrnd_interpol ) {
//     ifstream in;
//     in.open("bin/TestFiles/Photo_dNdxrnd_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double dNdxrnd_new;
//     double energy;
//     double dNdxrnd;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     int para;
//     double precision = 1E-2;

//     double energy_old;
//     bool first = true;
//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

//     int i=-1;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdxrnd;
//         first = false;
//         energy_old = -1;

//         i++;

//         CombOfPhoto.at(i)->EnableDNdxInterpolation();
//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
//             dNdxrnd_new=CombOfPhoto.at(i)->CalculatedNdx(Rand->rnd());


//             if(dNdxrnd!=0){
//                 if(log10(fabs(1-dNdxrnd_new/dNdxrnd))>-6)
//                 {
// //                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
// //                    cout << energy << "\t" << log10(fabs(1-dNdxrnd_new/dNdxrnd)) << endl;
//                 }
//             }

//             ASSERT_NEAR(dNdxrnd_new, dNdxrnd, precision*dNdxrnd);

//             in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>dNdxrnd;
//         }
//     }
//     delete Rand;
// }

// //----------------------------------------------------------------------------//

// TEST(Photonuclear , Test_of_e_interpol ) {
//     ifstream in;
//     in.open("bin/TestFiles/Photo_e_interpol.txt");

//     char firstLine[256];
//     in.getline(firstLine,256);
//     double e_new;
//     double energy;
//     double e;
//     double ecut;
//     double vcut;
//     string mediumName;
//     string particleName;
//     int para;
//     double precision = 1E-2;

//     double energy_old;
//     bool first = true;
//     RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
//     RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
//     Rand2->rnd();
//     double rnd1,rnd2;
//     int i=-1;
//     while(in.good())
//     {
//         if(first)in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>e;
//         first = false;
//         energy_old = -1;

//         i++;

//         CombOfPhoto.at(i)->EnableDNdxInterpolation();

//         while(energy_old < energy)
//         {
//             energy_old = energy;

//             rnd1 = Rand->rnd();
//             rnd2 = Rand2->rnd();

//             CombOfPhoto.at(i)->GetParticle()->SetEnergy(energy);
//             e_new=CombOfPhoto.at(i)->CalculateStochasticLoss(rnd1,rnd2);


//             if(e!=0){
//                 if(log10(fabs(1-e_new/e))>-14)
//                 {
// //                    cout <<ecut << "\t" << vcut<< "\t" << lpm<< "\t" << energy<< "\t" << mediumName<< "\t" << particleName << endl;
// //                    cout << energy << "\t" << log10(fabs(1-e_new/e)) << endl;
//                 }
//             }

//             ASSERT_NEAR(e_new, e, precision*e);

//             in>>para>>ecut>>vcut>>energy>>mediumName>>particleName>>e;
//         }
//     }
//     delete Rand;
//     delete Rand2;
// }

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
