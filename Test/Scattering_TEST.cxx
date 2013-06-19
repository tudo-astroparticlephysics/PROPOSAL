#include "gtest/gtest.h"
#include "PROPOSAL/Scattering.h"
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
    /*
    EnergyCutSettings A;
    EnergyCutSettings B =A;

    EXPECT_TRUE(A==B);
    */

}

TEST(Assignment , Copyconstructor2 ) {
    /*
    EnergyCutSettings A(5000,0.1);
    EnergyCutSettings B(A);

    EXPECT_TRUE(A==B);
    */

}

TEST(Assignment , Operator ) {
   /*
    EnergyCutSettings A;
    EnergyCutSettings B(200,0.01);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);

    A.SetEcut(300);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);
    */
}



int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
