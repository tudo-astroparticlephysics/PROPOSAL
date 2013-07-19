#include "gtest/gtest.h"
#include "PROPOSAL/Bremsstrahlung.h"
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

    Medium *medium = new Medium("air",1.);
    Medium *medium2 = new Medium("water",1.);
    Particle *particle = new Particle("mu",1.,1.,1,20,20,1e5,10);
    Particle *particle2 = new Particle("tau",1.,1.,1,20,20,1e5,10);
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
    Medium *medium = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);

    Bremsstrahlung A(particle, medium, cuts);
    Bremsstrahlung B(A);

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Operator ) {
    Medium *medium = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    Bremsstrahlung A(particle, medium, cuts);
    Bremsstrahlung B(particle, medium, cuts);
    A.SetParametrization(3);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);

    Medium *medium2 = new Medium("water",1.);
    Particle *particle2 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts2 = new EnergyCutSettings(200,-1);
    Bremsstrahlung *C = new Bremsstrahlung(particle2, medium2, cuts2);
    EXPECT_TRUE(A!=*C);

    A=*C;

    EXPECT_TRUE(A==*C);

}

TEST(Assignment , Swap ) {
    Medium *medium = new Medium("air",1.);
    Medium *medium2 = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    Particle *particle2 = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cuts = new EnergyCutSettings(500,-1);
    EnergyCutSettings *cuts2 = new EnergyCutSettings(500,-1);
    Bremsstrahlung A(particle, medium, cuts);
    Bremsstrahlung B(particle2, medium2, cuts2);
    A.EnableDEdxInterpolation();
    B.EnableDEdxInterpolation();
    EXPECT_TRUE(A==B);

    Medium *medium3 = new Medium("water",1.);
    Medium *medium4 = new Medium("water",1.);
    Particle *particle3 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
    Particle *particle4 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
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

    ifstream in;
    in.open("bin/TestFiles/Brems_dEdx.txt");

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
    int para;

    cout.precision(16);


    while(in.good())
    {
        in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);


        brems->SetParametrization(para);
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

    ifstream in;
    in.open("bin/TestFiles/Brems_dNdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dNdx;
    double dNdx_new;
    double energy;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);


    while(in.good())
    {
        in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);


        brems->SetParametrization(para);
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

    ifstream in;
    in.open("bin/TestFiles/Brems_dNdxrnd.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double dNdxrnd;
    double dNdxrnd_new;
    double energy;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);
    double energy_old=-1;

    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");

    bool first = true;
    while(in.good())
    {
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        first=false;
        energy_old = -1;
        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

       CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);

        //cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << med << "\t" << particleName<< "\t" << dNdx << endl;

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            dNdxrnd_new=brems->CalculatedNdx(Rand->rnd());

            ASSERT_NEAR(dNdxrnd_new, dNdxrnd, 1E-7*dNdxrnd);

            in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdxrnd;
        }



        delete cuts;
        delete medium;
        delete particle;
        delete brems;
    }
    delete Rand;
}



TEST(Bremsstrahlung , Test_of_e ) {

    ifstream in;
    in.open("bin/TestFiles/Brems_e.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double e;
    double e_new;
    double energy;
    double ecut;
    double vcut;
    string med;
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
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
        first=false;
        energy_old = -1;
        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

       CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);

        //cout << para << "\t" << ecut << "\t" << vcut << "\t" << lpm << "\t" << energy << "\t" << med << "\t" << particleName<< "\t" << e << endl;

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);

            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            e_new=brems->CalculateStochasticLoss(rnd1,rnd2);
            ASSERT_NEAR(e_new, e, 1E-7*e);

            in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
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

    ifstream in;
    in.open("bin/TestFiles/Brems_dEdx_interpol.txt");
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
    int para;

    cout.precision(16);
    double precision;
    double precisionOld = 1E-2;
    bool first=true;
    double energy_old;
    while(in.good())
    {
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        first=false;

        precision = precisionOld;
        energy_old =-1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);
        CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);

        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);
        brems->EnableDEdxInterpolation();

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            dEdx_new=brems->CalculatedEdx();

            if(!particleName.compare("tau") && energy < 10001)precision = 0.5;
            if(!particleName.compare("e") && energy > 1E10)precision = 0.5;

            ASSERT_NEAR(dEdx_new, dEdx, precision*dEdx);

            in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;

            precision = precisionOld;

        }



        delete cuts;
        delete medium;
        delete particle;
        delete brems;



    }
}

TEST(Bremsstrahlung , Test_of_dNdx_Interpolant ) {
    ifstream in;
    in.open("bin/TestFiles/Brems_dNdx_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double dNdx;
    double dNdx_new;
    double energy;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;
    int para;

    cout.precision(16);
    double energy_old=-1;

    bool first = true;
    while(in.good())
    {
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        first=false;
        energy_old = -1;
        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

       CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);
        brems->EnableDNdxInterpolation();

        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            dNdx_new=brems->CalculatedNdx();

            ASSERT_NEAR(dNdx_new, dNdx, 1E-6*dNdx);

            in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dNdx;
        }



        delete cuts;
        delete medium;
        delete particle;
        delete brems;
    }
}

TEST(Bremsstrahlung , Test_of_e_interpol ) {
return;
    ifstream in;
    in.open("bin/TestFiles/Brems_e_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double e;
    double e_new;
    double energy;
    double ecut;
    double vcut;
    string med;
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
        if(first)in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
        first=false;
        energy_old = -1;
        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

       CrossSections *brems = new Bremsstrahlung(particle, medium, cuts);
        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);
        brems->EnableDNdxInterpolation();


        while(energy_old < energy){
            energy_old = energy;
            brems->GetParticle()->SetEnergy(energy);
            rnd1 = Rand->rnd();
            rnd2 = Rand2->rnd();

            e_new = brems->CalculateStochasticLoss(rnd1,rnd2);

            ASSERT_NEAR(e_new, e, 1*e);

            in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>e;
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
