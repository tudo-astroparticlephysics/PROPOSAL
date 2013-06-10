#include "gtest/gtest.h"
#include "PROPOSAL/Geometry.h"
#include <iostream>
#include <string>
#include <cmath>
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/MathModel.h"

using namespace std;

TEST(Comparison , Comparison_equal ) {
    Geometry A;
    Geometry B;
    EXPECT_TRUE(A==B);
    Geometry* C = new Geometry();
    Geometry* D = new Geometry();
    C->InitSphere(0.,1.,2.,10.,0.);
    D->InitSphere(0.,1.,2.,10.,0.);
    EXPECT_TRUE(*C==*D);
    D->InitSphere(0.,0.,0.,0.,0.);
    EXPECT_TRUE(A==*D);

}

TEST(Comparison , Comparison_not_equal ) {
    Geometry A;
    Geometry B;
    B.InitBox(1.,2.,3.,1.,3.,2.);
    EXPECT_TRUE(A!=B);
    Geometry* C = new Geometry();
    Geometry* D = new Geometry();
    C->InitCylinder(1.,2.,3.,4.,3.,3.);
    D->InitBox(1.,2.,3.,4.,5.,6.);
    EXPECT_TRUE(*C!=*D);
}

TEST(Assignment , Copyconstructor ) {
    Geometry A;
    Geometry B =A;

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Copyconstructor2 ) {
    Geometry A;
    A.InitBox(1.,2.,3.,4.,5.,6.);
    Geometry B(A);

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Operator ) {
    Geometry A;
    Geometry B;
    A.InitSphere(1.,2.,3.,4.,3.);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);
}

TEST(Assignment , Swap ) {
    Geometry A;
    Geometry B;
    EXPECT_TRUE(A==B);
    Geometry* C = new Geometry();
    Geometry* D = new Geometry();
    C->InitSphere(1.,2.,3.,4.,3.);
    D->InitSphere(1.,2.,3.,4.,3.);

    EXPECT_TRUE(*C==*D);

    A.swap(*C);
    EXPECT_TRUE(A==*D);
    EXPECT_TRUE(B==*C);


}

TEST(IsInside , Box ) {
    Geometry A;

    double x        =   0;
    double y        =   0;
    double z        =   0;
    double theta    =   0;
    double phi      =   0;

    double rnd_x;
    double rnd_y;
    double rnd_z;

    double width_x  =   10;
    double width_y  =   10;
    double hight    =   10;

    double x0   =   3;
    double y0   =   3;
    double z0   =   3;


    A.InitBox(x0,y0,z0,width_x,width_y,hight);

    Particle * particle = new Particle("mu",x,y,z,theta,phi,0,0);
    int counter =0;
    int counter2 =0;

    MathModel M;

    for(int i = 0; i<1000000; i++)
    {
        rnd_x   = M.RandomDouble();
        rnd_y   = M.RandomDouble();
        rnd_z   = M.RandomDouble();

        rnd_x   = M.RandomDouble();
        rnd_y   = M.RandomDouble();
        rnd_z   = M.RandomDouble();

        x   =   ( 2 * rnd_x  - 1)*width_x;
        y   =   ( 2 * rnd_y  - 1)*width_y;
        z   =   ( 2 * rnd_z  - 1)*hight;

        particle->SetX(x);
        particle->SetY(y);
        particle->SetZ(z);

        if( particle->GetX() > x0 - 0.5*width_x &&
            particle->GetX() < x0 + 0.5*width_x &&
            particle->GetY() > y0 - 0.5*width_y &&
            particle->GetY() < y0 + 0.5*width_y &&
            particle->GetZ() > z0 - 0.5*hight   &&
            particle->GetZ() < z0 + 0.5*hight )
        {
            counter++;
            EXPECT_TRUE(A.IsParticleInside(particle));
        }
        else
        {
            counter2++;
            EXPECT_FALSE(A.IsParticleInside(particle));
        }

    }
    cout<<counter<<"\t"<<counter2<<endl;
}


TEST(IsInside , Sphere ) {
    Geometry A;

    double x        =   0;
    double y        =   0;
    double z        =   0;
    double theta    =   0;
    double phi      =   0;

    double rnd_x;
    double rnd_y;
    double rnd_z;

    double radius       =   60;
    double inner_radius =   0;

    A.InitSphere(0.,0.,0.,radius,inner_radius);

    Particle * particle = new Particle("mu",x,y,z,theta,phi,0,0);
    int counter =0;
    int counter2 =0;

    MathModel M;

    for(int i = 0; i<10000000; i++)
    {
        rnd_x   = M.RandomDouble();
        rnd_y   = M.RandomDouble();
        rnd_z   = M.RandomDouble();

        rnd_x   = M.RandomDouble();
        rnd_y   = M.RandomDouble();
        rnd_z   = M.RandomDouble();

        x   =   ( 2 * rnd_x  - 1)*radius;
        y   =   ( 2 * rnd_y  - 1)*radius;
        z   =   ( 2 * rnd_z  - 1)*radius;

        particle->SetX(x);
        particle->SetY(y);
        particle->SetZ(z);

        if( sqrt( x*x + y*y + z*z) < radius && sqrt( x*x + y*y + z*z) > inner_radius)
        {
            counter++;
            EXPECT_TRUE(A.IsParticleInside(particle));
        }
        else
        {
            counter2++;
            EXPECT_FALSE(A.IsParticleInside(particle));
        }

    }
    cout<<counter<<"\t"<<counter2<<endl;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




