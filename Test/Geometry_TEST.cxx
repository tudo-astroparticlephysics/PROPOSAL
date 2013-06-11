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

    double rnd_x0;
    double rnd_y0;
    double rnd_z0;

    double rnd_theta;
    double rnd_phi;

    double width_x  =   10;
    double width_y  =   10;
    double height   =   10;

    double big_width_x  =   4*width_x;
    double big_width_y  =   4*width_y;
    double big_height   =   4*height;

    double x0   =   0;
    double y0   =   0;
    double z0   =   0;

    Particle * particle = new Particle("mu",x,y,z,theta,phi,0,0);
    int is_inside  =0;
    int is_outside =0;

    double volumia_ratio =0;

    MathModel M;
    int number_particles = 1e5;
    int number_volumina  = 1e2;



    for(int i = 0; i<number_volumina; i++)
    {

        // Chose the origin of the box-geometry
        // This box should be inside the big box in which the particle
        // will be located
        rnd_x0  = M.RandomDouble();
        rnd_y0  = M.RandomDouble();
        rnd_z0  = M.RandomDouble();

        x0   =   ( 2 * rnd_x  - 1)* 0.5 *( big_width_x - width_x );
        y0   =   ( 2 * rnd_y  - 1)* 0.5 *( big_width_y - width_y );
        z0   =   ( 2 * rnd_z  - 1)* 0.5 *( big_height  - height  );

        A.InitBox(x0,y0,z0,width_x,width_y,height);

        volumia_ratio = width_x*width_y*height /( big_width_x*big_width_y*big_height);
        for(int j = 0; j<number_particles; j++)
        {

            // Chose particle location and angle
            rnd_x   = M.RandomDouble();
            rnd_y   = M.RandomDouble();
            rnd_z   = M.RandomDouble();

            rnd_theta   = M.RandomDouble();
            rnd_phi     = M.RandomDouble();

            theta   = rnd_theta*180;
            phi     = rnd_phi*360;

            x   =   ( 2 * rnd_x  - 1)* 0.5 * big_width_x;
            y   =   ( 2 * rnd_y  - 1)* 0.5 * big_width_y;
            z   =   ( 2 * rnd_z  - 1)* 0.5 * big_height;

            particle->SetX(x);
            particle->SetY(y);
            particle->SetZ(z);
            particle->SetTheta(theta);
            particle->SetPhi(phi);

            // if this constraints are true the particle is inside the box geometry
            if( particle->GetX() > x0 - 0.5*width_x &&
                particle->GetX() < x0 + 0.5*width_x &&
                particle->GetY() > y0 - 0.5*width_y &&
                particle->GetY() < y0 + 0.5*width_y &&
                particle->GetZ() > z0 - 0.5*height   &&
                particle->GetZ() < z0 + 0.5*height )
            {
                is_inside++;
                EXPECT_TRUE(A.IsParticleInside(particle));
            }
            else
            {
                is_outside++;
                EXPECT_FALSE(A.IsParticleInside(particle));
            }
        }
        ASSERT_NEAR(1.*is_inside ,volumia_ratio*number_particles , 3*sqrt(volumia_ratio*number_particles) );
        is_inside   = 0;
        is_outside  = 0;
    }
    // Check what happens if particles are on the border of the box

    A.InitBox(0,0,0,width_x,width_y,height);

    // Particle is on the top surface.
    // Theta 0° - 90° means particle is moving outside
    // This should be treated as outside
    // Theta 90° - 180° means particle is moving inside (should be treated as inside)
    // The value of phi does not matter
    particle->SetX(0);
    particle->SetY(0);
    particle->SetZ(0.5*height);
    for(int i = 0; i<1e4; i++)
    {
        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        if(theta < 90)
            EXPECT_FALSE(A.IsParticleInside(particle));
        if(theta > 90)
            EXPECT_TRUE(A.IsParticleInside(particle));

    }

    //Make this test for every surface of the box

    // bottom
    particle->SetX(0);
    particle->SetY(0);
    particle->SetZ(-0.5*height);
    particle->SetPhi(0);
    for(int i = 0; i<1e4; i++)
    {
        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        if(theta > 90)
            EXPECT_FALSE(A.IsParticleInside(particle));
        if(theta < 90)
            EXPECT_TRUE(A.IsParticleInside(particle));

    }

    // Surface in positiv x direction
    particle->SetX(0.5*width_x);
    particle->SetY(0);
    particle->SetZ(0);
    for(int i = 0; i<1e4; i++)
    {
        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        // phi = 0 is in positive x direction
        if(phi < 90 || phi > 270)
            EXPECT_FALSE(A.IsParticleInside(particle));
        else
            EXPECT_TRUE(A.IsParticleInside(particle));

    }
    // Surface in negativ x direction
    particle->SetX(-0.5*width_x);
    particle->SetY(0);
    particle->SetZ(0);
    for(int i = 0; i<1e4; i++)
    {
        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        // phi = 0 is in positive x direction
        if(phi < 90 || phi > 270)
            EXPECT_TRUE(A.IsParticleInside(particle));
        else
            EXPECT_FALSE(A.IsParticleInside(particle));

    }
    // Surface in positiv y direction
    particle->SetX(0);
    particle->SetY(0.5*width_y);
    particle->SetZ(0);
    for(int i = 0; i<1e4; i++)
    {
        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        // phi = 0 is in positive x direction
        if(phi < 180)
            EXPECT_FALSE(A.IsParticleInside(particle));
        else
            EXPECT_TRUE(A.IsParticleInside(particle));

    }
    // Surface in negativ y direction
    particle->SetX(0);
    particle->SetY(-0.5*width_y);
    particle->SetZ(0);
    for(int i = 0; i<1e4; i++)
    {
        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        // phi = 0 is in positive x direction
        if(phi < 180)
            EXPECT_TRUE(A.IsParticleInside(particle));
        else
            EXPECT_FALSE(A.IsParticleInside(particle));

    }

    //For completness check one corner
    particle->SetX(0.5*width_x);
    particle->SetY(0.5*width_y);
    particle->SetZ(0.5*height);
    for(int i = 0; i<1e4; i++)
    {
        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        if(theta < 90 || phi <180 || phi > 270)
            EXPECT_FALSE(A.IsParticleInside(particle));
        else
            EXPECT_TRUE(A.IsParticleInside(particle));

    }

}


//TEST(IsInside , Cylinder ) {
//    Geometry A;

//    double x        =   0;
//    double y        =   0;
//    double z        =   0;
//    double theta    =   90;
//    double phi      =   0;

//    double rnd_x;
//    double rnd_y;
//    double rnd_z;

//    double radius       =   10;
//    double inner_radius =   0;
//    double hight        =   10;

//    double x0   =   0;
//    double y0   =   0;
//    double z0   =   0;


//    A.InitCylinder(x0,y0,z0,radius,inner_radius,hight);

//    Particle * particle = new Particle("mu",x,y,z,theta,phi,0,0);

//    int counter =0;
//    int counter2 =0;

//    MathModel M;
//    M.set_seed(3838);

//    for(int i = 0; i<1e7; i++)
//    {
//        rnd_x   = M.RandomDouble();
//        rnd_y   = M.RandomDouble();
//        rnd_z   = M.RandomDouble();

//        rnd_x   = M.RandomDouble();
//        rnd_y   = M.RandomDouble();
//        rnd_z   = M.RandomDouble();

//        x   =   ( 2 * rnd_x  - 1)*radius;
//        y   =   ( 2 * rnd_y  - 1)*radius;
//        z   =   ( 2 * rnd_z  - 1)*hight;

//        particle->SetX(x);
//        particle->SetY(y);
//        particle->SetZ(z);

////        if(A.IsParticleInside(particle))counter++;
////        else counter2++;
//        if( sqrt( x*x + y*y) < radius &&
//            sqrt( x*x + y*y) > inner_radius &&
//            z > z0 - 0.5*hight  &&
//            z < z0 + 0.5*hight )
//        {
//            counter++;
//            EXPECT_TRUE(A.IsParticleInside(particle));
//        }
//        else
//        {
//            counter2++;
//            EXPECT_FALSE(A.IsParticleInside(particle));
//        }

//    }
//    cout<<counter<<"\t"<<counter2<<endl;
//}

//TEST(IsInside , Sphere ) {
//    Geometry A;

//    double x        =   0;
//    double y        =   0;
//    double z        =   0;
//    double theta    =   0;
//    double phi      =   0;

//    double rnd_x;
//    double rnd_y;
//    double rnd_z;

//    double radius       =   60;
//    double inner_radius =   0;

//    A.InitSphere(0.,0.,0.,radius,inner_radius);

//    Particle * particle = new Particle("mu",x,y,z,theta,phi,0,0);
//    int counter =0;
//    int counter2 =0;

//    MathModel M;

//    for(int i = 0; i<10000000; i++)
//    {
//        rnd_x   = M.RandomDouble();
//        rnd_y   = M.RandomDouble();
//        rnd_z   = M.RandomDouble();

//        rnd_x   = M.RandomDouble();
//        rnd_y   = M.RandomDouble();
//        rnd_z   = M.RandomDouble();

//        x   =   ( 2 * rnd_x  - 1)*radius;
//        y   =   ( 2 * rnd_y  - 1)*radius;
//        z   =   ( 2 * rnd_z  - 1)*radius;

//        particle->SetX(x);
//        particle->SetY(y);
//        particle->SetZ(z);

//        if( sqrt( x*x + y*y + z*z) < radius && sqrt( x*x + y*y + z*z) > inner_radius)
//        {
//            counter++;
//            EXPECT_TRUE(A.IsParticleInside(particle));
//        }
//        else
//        {
//            counter2++;
//            EXPECT_FALSE(A.IsParticleInside(particle));
//        }

//    }
//    cout<<counter<<"\t"<<counter2<<endl;
//}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




