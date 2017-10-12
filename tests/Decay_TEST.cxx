
// #include <iostream>

#include "gtest/gtest.h"

#include "PROPOSAL/Decay.h"
#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Constants.h"


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

    Particle *particle = new Particle(ParticleType::MuMinus,position,direction,1e5,10);
    Decay *A = new Decay(particle);
    Decay *B = new Decay(particle);
    EXPECT_TRUE(*A==*B);

    Decay *C = new Decay();
    Decay *D = new Decay();

    EXPECT_TRUE(*C==*D);

    A->GetParticle()->SetEnergy(1e6);
    B->GetParticle()->SetEnergy(1e6);
    EXPECT_TRUE(*A==*B);

}

TEST(Comparison , Comparison_not_equal ) {
    direction.SetSphericalCoordinates(1,20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    Particle *particle = new Particle(ParticleType::MuMinus,position,direction,1e5,10);
    Particle *particle2 = new Particle(ParticleType::TauMinus,position,direction,1e5,10);
    Decay *A = new Decay(particle);
    Decay *B = new Decay(particle2);

    EXPECT_TRUE(*A!=*B);
    A->SetParticle(particle2);
    EXPECT_TRUE(*A==*B);

    A->SetMultiplier(1.2);
    EXPECT_TRUE(*A!=*B);
    B->SetMultiplier(1.2);

    EXPECT_TRUE(*A==*B);
    A->GetRootFinder()->SetMaxSteps(5);
    EXPECT_TRUE(*A!=*B);
}

TEST(Assignment , Copyconstructor ) {
    Decay A;
    Decay B =A;
    EXPECT_TRUE(A==B);

}

TEST(Assignment , Copyconstructor2 ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();
    Particle *particle = new Particle(ParticleType::MuMinus,position,direction,1e5,10);

    Decay A(particle);
    Decay B(A);

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Operator ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();
    Particle *particle = new Particle(ParticleType::MuMinus,position,direction,1e5,10);
    Particle *particle2 = new Particle(ParticleType::TauMinus,position,direction,1e5,10);
    Decay A(particle);
    Decay B(particle2);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);

    Particle *particle3 = new Particle(ParticleType::TauMinus,position,direction,1e5,10);
    Decay *C = new Decay(particle3);
    EXPECT_TRUE(A!=*C);

    A=*C;

    EXPECT_TRUE(A==*C);

}

TEST(Assignment , Swap ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    Particle *particle = new Particle(ParticleType::MuMinus,position,direction,1e5,10);
    Particle *particle2 = new Particle(ParticleType::MuMinus,position,direction,1e5,10);

    Decay A(particle);
    Decay B(particle2);
    EXPECT_TRUE(A==B);

    Particle *particle3 = new Particle(ParticleType::TauMinus,position,direction,1e5,10);
    Particle *particle4 = new Particle(ParticleType::TauMinus,position,direction,1e5,10);

    Decay *C = new Decay(particle3);
    Decay *D = new Decay(particle4);
    C->GetRootFinder()->SetMaxSteps(5);
    D->GetRootFinder()->SetMaxSteps(5);

    EXPECT_TRUE(*C==*D);

    A.swap(*C);

    EXPECT_TRUE(A==*D);
    EXPECT_TRUE(*C==B);


}

TEST(Decay , decay ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    ifstream in;
    in.open("bin/TestFiles/Decay_decay.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double energy;
    double decay;
    double decay_new;
    string particleName;

    cout.precision(16);


    while(in.good())
    {
        in>>particleName>>energy>>decay;

        Particle *particle = new Particle(Particle::GetTypeFromName(particleName),position,direction,1e5,10);
        particle->SetEnergy(energy);

        Decay* dec  = new Decay(particle);

        decay_new = dec->MakeDecay();

        ASSERT_NEAR(decay_new, decay, 1e-14*decay);

        delete particle;
        delete dec;
    }
}


TEST(Decay , ProductEnergy ) {
    direction.SetSphericalCoordinates(1,.20*PI/180.,20*PI/180.);
    direction.CalculateCartesianFromSpherical();

    ifstream in;
    in.open("bin/TestFiles/Decay_ProductEnergy.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double energy;
    double product_energy;
    double product_energy_new;
    string particleName;

    cout.precision(16);

    RndFromFile* Rand1 = new RndFromFile("bin/TestFiles/rnd.txt");
    double rnd1;

    RndFromFile* Rand2 = new RndFromFile("bin/TestFiles/rnd.txt");
    double rnd2 = Rand2->rnd();

    RndFromFile* Rand3 = new RndFromFile("bin/TestFiles/rnd.txt");
    double rnd3= Rand3->rnd();
    rnd3= Rand3->rnd();

    while(in.good())
    {
        in>>particleName>>energy>>product_energy;

        Particle *particle = new Particle(Particle::GetTypeFromName(particleName),position,direction,1e5,10);
        particle->SetEnergy(energy);

        rnd1= Rand1->rnd();
        rnd2= Rand2->rnd();
        rnd3= Rand3->rnd();

        Decay* dec  = new Decay(particle);

        product_energy_new = dec->CalculateProductEnergy(rnd1,rnd2,rnd3);

        ASSERT_NEAR(product_energy_new, product_energy, 1e-14*product_energy);

        delete particle;
        delete dec;
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
