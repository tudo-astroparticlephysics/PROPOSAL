#include "gtest/gtest.h"
#include "PROPOSAL/ContinuousRandomization.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>
#include <vector>
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/EnergyCutSettings.h"
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

    double dNdx;

    Medium *medium = new Medium("hydrogen",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cut_settings = new EnergyCutSettings(500,-1);

    vector<CrossSections*> crosssections;

    crosssections.resize(4);
    crosssections.at(0) = new Ionization(particle, medium, cut_settings);
    crosssections.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
    crosssections.at(2) = new Photonuclear(particle, medium, cut_settings);
    crosssections.at(3) = new Epairproduction(particle, medium, cut_settings);


    ContinuousRandomization *A = new ContinuousRandomization(particle, medium, crosssections);
    ContinuousRandomization *B = new ContinuousRandomization(particle, medium, crosssections);
    EXPECT_TRUE(*A==*B);

    A->EnableDE2dxInterpolation();
    B->EnableDE2dxInterpolation();

    EXPECT_TRUE(*A==*B);

    ContinuousRandomization *C = new ContinuousRandomization();
    ContinuousRandomization *D = new ContinuousRandomization();

    EXPECT_TRUE(*C==*D);


    A->GetParticle()->SetEnergy(1e6);
    B->GetParticle()->SetEnergy(1e6);
    EXPECT_TRUE(*A==*B);

    dNdx = A->GetCrosssections().at(0)->CalculatedNdx();
    dNdx = B->GetCrosssections().at(0)->CalculatedNdx();
    EXPECT_TRUE(*A==*B);

}

TEST(Comparison , Comparison_not_equal ) {
    Medium *medium = new Medium("air",1.);
    Medium *medium2 = new Medium("water",1.);
    Particle *particle = new Particle("mu",1.,1.,1,20,20,1e5,10);
    Particle *particle2 = new Particle("tau",1.,1.,1,20,20,1e5,10);
    EnergyCutSettings *cut_settings = new EnergyCutSettings(500,-1);

    vector<CrossSections*> crosssections;

    crosssections.resize(4);
    crosssections.at(0) = new Ionization(particle, medium, cut_settings);
    crosssections.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
    crosssections.at(2) = new Photonuclear(particle, medium, cut_settings);
    crosssections.at(3) = new Epairproduction(particle, medium, cut_settings);

    vector<CrossSections*> crosssections2;

    crosssections2.resize(4);
    crosssections2.at(0) = new Ionization(particle, medium, cut_settings);
    crosssections2.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
    crosssections2.at(2) = new Photonuclear(particle, medium, cut_settings);
    crosssections2.at(3) = new Ionization(particle, medium, cut_settings);

    ContinuousRandomization *A = new ContinuousRandomization(particle, medium, crosssections);
    ContinuousRandomization *B = new ContinuousRandomization(particle2, medium, crosssections);
    ContinuousRandomization *C = new ContinuousRandomization(particle, medium2, crosssections);
    ContinuousRandomization *D = new ContinuousRandomization(particle, medium2, crosssections2);
    ContinuousRandomization *E = new ContinuousRandomization(particle, medium2, crosssections);

    EXPECT_TRUE(*A!=*B);
    EXPECT_TRUE(*C!=*D);
    EXPECT_TRUE(*B!=*D);
    EXPECT_TRUE(*E==*C);

    E->SetParticle(particle2);
    EXPECT_TRUE(*C!=*E);
    C->SetParticle(particle2);
    EXPECT_TRUE(*C==*E);
    C->GetCrosssections().at(2)->SetParametrization(6);

    EXPECT_TRUE(*D!=*E);


}

TEST(Assignment , Copyconstructor ) {
    ContinuousRandomization A;
    ContinuousRandomization B =A;

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Copyconstructor2 ) {
    Medium *medium = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cut_settings = new EnergyCutSettings(500,-1);

    vector<CrossSections*> crosssections;

    crosssections.resize(4);
    crosssections.at(0) = new Ionization(particle, medium, cut_settings);
    crosssections.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
    crosssections.at(2) = new Photonuclear(particle, medium, cut_settings);
    crosssections.at(3) = new Epairproduction(particle, medium, cut_settings);

    ContinuousRandomization A(particle, medium, crosssections);
    ContinuousRandomization B(A);

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Operator ) {
    Medium *medium = new Medium("air",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cut_settings = new EnergyCutSettings(500,-1);

    vector<CrossSections*> crosssections;

    crosssections.resize(4);
    crosssections.at(0) = new Ionization(particle, medium, cut_settings);
    crosssections.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
    crosssections.at(2) = new Photonuclear(particle, medium, cut_settings);
    crosssections.at(3) = new Epairproduction(particle, medium, cut_settings);

    vector<CrossSections*> crosssections2;

    crosssections2.resize(4);
    crosssections2.at(0) = new Ionization(particle, medium, cut_settings);
    crosssections2.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
    crosssections2.at(2) = new Photonuclear(particle, medium, cut_settings);
    crosssections2.at(3) = new Epairproduction(particle, medium, cut_settings);

    ContinuousRandomization A(particle, medium, crosssections);
    ContinuousRandomization B(particle, medium, crosssections2);
    B.GetCrosssections().at(2)->SetParametrization(6);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);

    vector<CrossSections*> crosssections3;

    crosssections3.resize(4);
    crosssections3.at(0) = new Ionization(particle, medium, cut_settings);
    crosssections3.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
    crosssections3.at(2) = new Photonuclear(particle, medium, cut_settings);
    crosssections3.at(3) = new Ionization(particle, medium, cut_settings);

    ContinuousRandomization *C = new ContinuousRandomization(particle, medium, crosssections3);
    EXPECT_TRUE(A!=*C);

    A=*C;

    EXPECT_TRUE(A==*C);

}

TEST(Assignment , Swap ) {
    Medium *medium = new Medium("hydrogen",1.);
    Medium *medium2 = new Medium("hydrogen",1.);
    Particle *particle = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    Particle *particle2 = new Particle("mu",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cut_settings = new EnergyCutSettings(500,0.05);
    EnergyCutSettings *cut_settings2 = new EnergyCutSettings(500,0.05);

    vector<CrossSections*> crosssections;

    crosssections.resize(4);
    crosssections.at(0) = new Ionization(particle, medium, cut_settings);
    crosssections.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
    crosssections.at(2) = new Photonuclear(particle, medium, cut_settings);
    crosssections.at(3) = new Epairproduction(particle, medium, cut_settings);


    vector<CrossSections*> crosssections2;

    crosssections2.resize(4);
    crosssections2.at(0) = new Ionization(particle2, medium2, cut_settings2);
    crosssections2.at(1) = new Bremsstrahlung(particle2, medium2, cut_settings2);
    crosssections2.at(2) = new Photonuclear(particle2, medium2, cut_settings2);
    crosssections2.at(3) = new Epairproduction(particle2, medium2, cut_settings2);

    ContinuousRandomization A(particle, medium, crosssections);
    ContinuousRandomization B(particle2, medium2, crosssections2);

    for(unsigned int i=0 ; i<crosssections.size();i++)
    {
        crosssections.at(i)->EnableDEdxInterpolation();
        crosssections.at(i)->EnableDNdxInterpolation();
    }
    for(unsigned int i=0 ; i<crosssections2.size();i++)
    {
        crosssections2.at(i)->EnableDEdxInterpolation();
        crosssections2.at(i)->EnableDNdxInterpolation();
    }

    A.EnableDE2dxInterpolation();
    B.EnableDE2dxInterpolation();
    A.EnableDE2deInterpolation();
    B.EnableDE2deInterpolation();
    EXPECT_TRUE(A==B);


    Medium *medium3 = new Medium("water",1.);
    Medium *medium4 = new Medium("water",1.);
    Particle *particle3 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
    Particle *particle4 = new Particle("tau",1.,1.,1,.20,20,1e5,10);
    EnergyCutSettings *cut_settings3 = new EnergyCutSettings(200,-1);
    EnergyCutSettings *cut_settings4 = new EnergyCutSettings(200,-1);

    vector<CrossSections*> crosssections3;

    crosssections3.resize(4);
    crosssections3.at(0) = new Ionization(particle3, medium3, cut_settings3);
    crosssections3.at(1) = new Bremsstrahlung(particle3, medium3, cut_settings3);
    crosssections3.at(2) = new Photonuclear(particle3, medium3, cut_settings3);
    crosssections3.at(3) = new Ionization(particle3, medium3, cut_settings3);


    vector<CrossSections*> crosssections4;

    crosssections4.resize(4);
    crosssections4.at(0) = new Ionization(particle4, medium4, cut_settings4);
    crosssections4.at(1) = new Bremsstrahlung(particle4, medium4, cut_settings4);
    crosssections4.at(2) = new Photonuclear(particle4, medium4, cut_settings4);
    crosssections4.at(3) = new Ionization(particle4, medium4, cut_settings4);

    ContinuousRandomization *C = new ContinuousRandomization(particle3, medium3, crosssections3);
    ContinuousRandomization *D = new ContinuousRandomization(particle4, medium4, crosssections4);
    EXPECT_TRUE(*C==*D);

    A.swap(*C);

    EXPECT_TRUE(A==*D);
    EXPECT_TRUE(*C==B);


}

//TEST(ContinuousRandomization , Randomize ) {

//    ifstream in;
//    in.open("bin/TestFiles/ContinuousRandomization.txt");

//    char firstLine[256];
//    in.getline(firstLine,256);
//    double initial_energy;
//    double final_energy;
//    double randomized_energy;
//    double randomized_energy_new;
//    double vcut;
//    double ecut;
//    string med;
//    string particleName;
//    double rnd;
//    cout.precision(16);

//    double max_diff = 0;
//    double diff = 0;
//    double sum_diff = 0;
//    int counter = 0;

//    while(in.good())
//    {
//        in>>rnd>>particleName>>med>>ecut>>vcut>>initial_energy>>final_energy>>randomized_energy;
//        counter++;
//        Medium *medium = new Medium(med,1.);
//        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
//        EnergyCutSettings *cut_settings = new EnergyCutSettings(ecut,vcut);

//        vector<CrossSections*> crosssections;

//        crosssections.resize(4);
//        crosssections.at(0) = new Ionization(particle, medium, cut_settings);
//        crosssections.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
//        crosssections.at(2) = new Photonuclear(particle, medium, cut_settings);
//        crosssections.at(3) = new Epairproduction(particle, medium, cut_settings);

//        ContinuousRandomization * cont = new ContinuousRandomization(particle,medium,crosssections);
//        randomized_energy_new = cont->Randomize(initial_energy,final_energy,rnd);

//        diff= abs(1 - randomized_energy_new/randomized_energy);
//        sum_diff += diff;

//        ASSERT_NEAR(randomized_energy_new, randomized_energy, 1e-1*randomized_energy);
//        cout<<randomized_energy_new<<"\t"<<randomized_energy<<"\t"<<diff<<endl;
//        if(diff>max_diff)
//        {
//            max_diff = diff;
//        }
//        for(unsigned int i = 0 ; i < crosssections.size() ; i++){

//            delete crosssections.at(i);
//        }

//        delete cut_settings;
//        delete medium;
//        delete particle;
//        delete cont;
//    }
//    cout<<"Maxi_diff "<<max_diff<<endl;
//    cout<<"Average diff "<<sum_diff/counter<<endl;

//}


TEST(ContinuousRandomization , Randomize_interpol ) {

    ifstream in;
    in.open("bin/TestFiles/ContinuousRandomization_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double initial_energy;
    double final_energy;
    double randomized_energy;
    double randomized_energy_new;
    double vcut;
    double ecut;
    string med;
    string particleName;
    double rnd;
    double energy_old;
    bool first = true;
    int i = -1;

    cout.precision(16);



    while(in.good())
    {
        if(first)in>>rnd>>particleName>>med>>ecut>>vcut>>initial_energy>>final_energy>>randomized_energy;
        first = false;
        energy_old = -1;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        EnergyCutSettings *cut_settings = new EnergyCutSettings(ecut,vcut);

        vector<CrossSections*> crosssections;

        crosssections.resize(4);
        crosssections.at(0) = new Ionization(particle, medium, cut_settings);
        crosssections.at(1) = new Bremsstrahlung(particle, medium, cut_settings);
        crosssections.at(2) = new Photonuclear(particle, medium, cut_settings);
        crosssections.at(3) = new Epairproduction(particle, medium, cut_settings);

        ContinuousRandomization * cont = new ContinuousRandomization(particle,medium,crosssections);
        for(unsigned int i=0 ; i<crosssections.size();i++)
        {
            crosssections.at(i)->EnableDEdxInterpolation();
            crosssections.at(i)->EnableDNdxInterpolation();
        }
        cont->EnableDE2dxInterpolation();
        cont->EnableDE2deInterpolation();

        while(energy_old < initial_energy)
        {

            energy_old = initial_energy;

            randomized_energy_new = cont->Randomize(initial_energy,final_energy,rnd);

            ASSERT_NEAR(randomized_energy_new, randomized_energy, 1e-1*randomized_energy);

            in>>rnd>>particleName>>med>>ecut>>vcut>>initial_energy>>final_energy>>randomized_energy;
        }

        for(unsigned int i = 0 ; i < crosssections.size() ; i++){

            delete crosssections.at(i);
        }

        delete cut_settings;
        delete medium;
        delete particle;
        delete cont;
    }

}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
