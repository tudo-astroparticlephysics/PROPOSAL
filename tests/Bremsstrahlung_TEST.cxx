
// #include <iostream>
// #include <string>

#include "gtest/gtest.h"

#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Constants.h"
// #include "PROPOSAL/CrossSections.h"

using namespace std;
using namespace PROPOSAL;

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

Vector3D position(1.,1.,1.);
Vector3D direction(0.,0.,0.);

TEST(Comparison , Comparison_equal ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    double dEdx;
    Medium* medium = new Air();
    PROPOSALParticle *particle = new PROPOSALParticle(ParticleType::MuMinus,position,direction,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    Bremsstrahlung *A = new Bremsstrahlung(particle, medium, cuts);
    Bremsstrahlung *B = new Bremsstrahlung(particle, medium, cuts);
    EXPECT_TRUE(*A==*B);

    Bremsstrahlung *C = new Bremsstrahlung();
    Bremsstrahlung *D = new Bremsstrahlung();

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
    direction.SetSphericalCoordinates(1,20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    Medium *medium = new Air();
    Medium *medium2 = new Water();
    PROPOSALParticle *particle = new PROPOSALParticle(ParticleType::MuMinus,position,direction,1e5,10);
    PROPOSALParticle *particle2 = new PROPOSALParticle(ParticleType::TauMinus,position,direction,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    Bremsstrahlung *A = new Bremsstrahlung(particle, medium, cuts);
    Bremsstrahlung *B = new Bremsstrahlung(particle, medium2, cuts);
    Bremsstrahlung *C = new Bremsstrahlung(particle2, medium, cuts);
    Bremsstrahlung *D = new Bremsstrahlung(particle2, medium2, cuts);
    Bremsstrahlung *E = new Bremsstrahlung(particle2, medium2, cuts);

    EXPECT_TRUE(*A!=*B);
    EXPECT_TRUE(*C!=*D);
    EXPECT_TRUE(*B!=*D);
    EXPECT_TRUE(*D==*E);

    E->SetParticle(particle);
    EXPECT_TRUE(*D!=*E);
    D->SetParticle(particle);
    EXPECT_TRUE(*D==*E);


}

TEST(Assignment , Copyconstructor ) {
    Bremsstrahlung A;
    Bremsstrahlung B =A;

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Copyconstructor2 ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();
    Medium *medium = new Air();
    PROPOSALParticle *particle = new PROPOSALParticle(ParticleType::MuMinus,position,direction,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);

    Bremsstrahlung A(particle, medium, cuts);
    Bremsstrahlung B(A);

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Operator ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();
    Medium *medium = new Air();
    PROPOSALParticle *particle = new PROPOSALParticle(ParticleType::MuMinus,position,direction,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    Bremsstrahlung A(particle, medium, cuts);
    Bremsstrahlung B(particle, medium, cuts);
    A.SetParametrization(ParametrizationType::BremsPetrukhinShestakov);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);

    Medium *medium2 = new Water();
    PROPOSALParticle *particle2 = new PROPOSALParticle(ParticleType::TauMinus,position,direction,1e5,10);
    EnergyCutSettings *cuts2 = new EnergyCutSettings(200,-1);
    Bremsstrahlung *C = new Bremsstrahlung(particle2, medium2, cuts2);
    EXPECT_TRUE(A!=*C);

    A=*C;

    EXPECT_TRUE(A==*C);

}

TEST(Assignment , Swap ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();
    Medium *medium = new Air();
    Medium *medium2 = new Air();
    PROPOSALParticle *particle = new PROPOSALParticle(ParticleType::MuMinus,position,direction,1e5,10);
    PROPOSALParticle *particle2 = new PROPOSALParticle(ParticleType::MuMinus,position,direction,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    EnergyCutSettings *cuts2 = new EnergyCutSettings(500,-1);
    Bremsstrahlung A(particle, medium, cuts);
    Bremsstrahlung B(particle2, medium2, cuts2);
    A.EnableDEdxInterpolation();
    B.EnableDEdxInterpolation();
    EXPECT_TRUE(A==B);

    Medium *medium3 = new Water();
    Medium *medium4 = new Water();
    PROPOSALParticle *particle3 = new PROPOSALParticle(ParticleType::TauMinus,position,direction,1e5,10);
    PROPOSALParticle *particle4 = new PROPOSALParticle(ParticleType::TauMinus,position,direction,1e5,10);
    EnergyCutSettings *cuts3 = new EnergyCutSettings(200,-1);
    EnergyCutSettings *cuts4 = new EnergyCutSettings(200,-1);
    Bremsstrahlung *C = new Bremsstrahlung(particle3, medium3, cuts3);
    Bremsstrahlung *D = new Bremsstrahlung(particle4, medium4, cuts4);
    EXPECT_TRUE(*C==*D);

    A.swap(*C);

    EXPECT_TRUE(A==*D);
    EXPECT_TRUE(*C==B);


}

TEST(Bremsstrahlung , Test_of_dEdx ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    ifstream in;
    in.open("bin/TestFiles/Brems_dEdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dEdx_new;
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string mediumName;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);


    while(in.good())
    {
        in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;

        Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
        PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);


        brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
        brems->EnableLpmEffect(lpm);



        dEdx_new=brems->CalculatedEdx();

        ASSERT_NEAR(dEdx_new, dEdx, 1e-7*dEdx);

        delete cuts;
        delete medium;
        delete particle;
        delete brems;
    }
}

TEST(Bremsstrahlung , Test_of_dNdx ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    ifstream in;
    in.open("bin/TestFiles/Brems_dNdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdx;
    double dNdx_new;
    double energy;
    double ecut;
    double vcut;
    string mediumName;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);


    while(in.good())
    {
        in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdx;

        Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
        PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);


        brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
        brems->EnableLpmEffect(lpm);

        dNdx_new=brems->CalculatedNdx();
        ASSERT_NEAR(dNdx_new, dNdx, 1e-7*dNdx);

        delete cuts;
        delete medium;
        delete particle;
        delete brems;



    }
}

TEST(Bremsstrahlung , Test_of_dNdxrnd ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    ifstream in;
    in.open("bin/TestFiles/Brems_dNdxrnd.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double dNdxrnd;
    double dNdxrnd_new;
    double energy;
    double ecut;
    double vcut;
    string mediumName;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);
    double energy_old=-1;

    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

    bool first = true;
    while(in.good())
    {
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdxrnd;
        first=false;
        energy_old = -1;
        Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
        PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
        brems->EnableLpmEffect(lpm);

        //cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << mediumName << "\t" << particleName<< "\t" << dNdxrnd << endl;

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=brems->CalculatedNdx(Rand->rnd());

            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, 1E-7*dNdxrnd);

            in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdxrnd;
        }



        delete cuts;
        delete medium;
        delete particle;
        delete brems;
    }
    delete Rand;
}



TEST(Bremsstrahlung , Test_of_e ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    ifstream in;
    in.open("bin/TestFiles/Brems_e.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double e;
    double e_new;
    double energy;
    double ecut;
    double vcut;
    string mediumName;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);
    double energy_old=-1;

    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
    RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
    Rand2->rnd();

    double rnd1, rnd2;
    bool first = true;
    while(in.good())
    {
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
        first=false;
        energy_old = -1;
        Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
        PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
        brems->EnableLpmEffect(lpm);

        //cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << mediumName << "\t" << particleName<< "\t" << e << endl;

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);

            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            e_new=brems->CalculateStochasticLoss(rnd1,rnd2);
            ASSERT_NEAR(e_new, e, 1E-7*e);

            in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
        }

        delete cuts;
        delete medium;
        delete particle;
        delete brems;
    }
    delete Rand2;
    delete Rand;
}


TEST(Bremsstrahlung , Test_of_dEdx_Interpolant ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    ifstream in;
    in.open("bin/TestFiles/Brems_dEdx_interpol.txt");
    char firstLine[256];
    in.getline(firstLine,256);
    double dEdx_new;
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string mediumName;
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
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;
        first=false;

        precision = precisionOld;
        energy_old =-1;

        Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
        PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);

        brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
        brems->EnableLpmEffect(lpm);
        brems->EnableDEdxInterpolation();

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            dEdx_new=brems->CalculatedEdx();

            if(!particleName.compare("tau") && energy < 10001)precision = 0.5;
            if(!particleName.compare("e") && energy > 1E10)precision = 0.5;

            ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

            in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dEdx;

            precision = precisionOld;

        }



        delete cuts;
        delete medium;
        delete particle;
        delete brems;



    }
}

TEST(Bremsstrahlung , Test_of_dNdx_Interpolant ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();
    ifstream in;
    in.open("bin/TestFiles/Brems_dNdx_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double dNdx;
    double dNdx_new;
    double energy;
    double ecut;
    double vcut;
    string mediumName;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);
    double energy_old=-1;

    bool first = true;
    while(in.good())
    {
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdx;
        first=false;
        energy_old = -1;
        Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
        PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
        brems->EnableLpmEffect(lpm);
        brems->EnableDNdxInterpolation();

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            dNdx_new=brems->CalculatedNdx();

            ASSERT_NEAR(dNdx_new, dNdx, 1E-6*dNdx);

            in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>dNdx;
        }



        delete cuts;
        delete medium;
        delete particle;
        delete brems;
    }
}

TEST(Bremsstrahlung , Test_of_e_interpol ) {
return;
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();
    ifstream in;
    in.open("bin/TestFiles/Brems_e_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double e;
    double e_new;
    double energy;
    double ecut;
    double vcut;
    string mediumName;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);
    double energy_old=-1;
    double precision = 1E-5;
    double precision_old = precision;
    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
    RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
    Rand2->rnd();

    double rnd1,rnd2;
    bool first = true;
    while(in.good())
    {
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
        first=false;
        energy_old = -1;
        Medium *medium = MediumFactory::Get()->CreateMedium(mediumName);
        PROPOSALParticle *particle = new PROPOSALParticle(PROPOSALParticle::GetTypeFromName(particleName),position,direction,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(static_cast<ParametrizationType::Enum>(para));
        brems->EnableLpmEffect(lpm);
        brems->EnableDNdxInterpolation();


        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            e_new = brems->CalculateStochasticLoss(rnd1,rnd2);

            ASSERT_NEAR(e_new, e, 1*e);

            in>>para>>ecut>>vcut>>lpm>>energy>>mediumName>>particleName>>e;
            precision = precision_old;
        }



        delete cuts;
        delete medium;
        delete particle;
        delete brems;
    }
    delete Rand2;
    delete Rand;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
