
// #include <iostream>
// #include <string>
// #include <cmath>

#include "gtest/gtest.h"

#include "PROPOSAL/Geometry.h"
#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/Constants.h"
// #include "PROPOSAL/PROPOSALParticle.h"

using namespace std;
using namespace PROPOSAL;

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

    PROPOSALParticle * particle = new PROPOSALParticle(ParticleType::MuMinus,x,y,z,theta,phi,0,0);
    int is_inside  =0;
    int is_outside =0;

    double volumia_ratio =0;

    MathModel M;
    int number_particles = 1e6;
    int number_volumina  = 1e1;



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

        // The values are divided by 100 to convert the units...
        // Init functions expects m but here everthing is in cm
        A.InitBox(x0/100,y0/100,z0/100,width_x/100,width_y/100,height/100);

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

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    A.InitBox(0,0,0,width_x/100,width_y/100,height/100);

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


TEST(IsInside , Cylinder ) {
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
    double rnd_inner_radius;

    double radius       =   10;
    double inner_radius =   0;
    double height       =   10;

    double big_width_x  =   4*radius;
    double big_width_y  =   4*radius;
    double big_height   =   4*height;

    double x0   =   0;
    double y0   =   0;
    double z0   =   0;

    PROPOSALParticle * particle = new PROPOSALParticle(ParticleType::MuMinus,x,y,z,theta,phi,0,0);
    int is_inside  =0;
    int is_outside =0;

    double volumia_ratio =0;

    MathModel M;
    int number_particles = 1e6;
    int number_volumina  = 1e1;

    for(int i = 0; i<number_volumina; i++)
    {

        // Chose the origin of the cylinder-geometry
        // This cylinder should be inside the big box in which the particle
        // will be located
        rnd_x0  = M.RandomDouble();
        rnd_y0  = M.RandomDouble();
        rnd_z0  = M.RandomDouble();

        x0   =   ( 2 * rnd_x0  - 1)* ( 0.5 * big_width_x - radius );
        y0   =   ( 2 * rnd_y0  - 1)* ( 0.5 * big_width_y - radius );
        z0   =   ( 2 * rnd_z0  - 1)* 0.5 *( big_height  - height  );

        x0   =   0;
        y0   =   0;
        z0   =   0;

        rnd_inner_radius    = M.RandomDouble();

        inner_radius    = radius*rnd_inner_radius;

        // The values are divided by 100 to convert the units...
        // Init functions expects m but here everthing is in cm
        A.InitCylinder(x0/100,y0/100,z0/100,radius/100,inner_radius/100,height/100);

        volumia_ratio = height*PI*( pow(radius,2) - pow(inner_radius,2) )
                /( big_width_x*big_width_y*big_height);
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

            // if this constraints are true the particle is inside the cylinder geometry
            if( sqrt( pow( (x-x0),2 ) + pow( (y-y0),2 ) ) < radius &&
                sqrt( pow( (x-x0),2 ) + pow( (y-y0),2 ) ) > inner_radius &&
                z > z0 - 0.5*height  &&
                z < z0 + 0.5*height )
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

    // Test borders

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    A.InitSphere(0,0,0,radius/100,0);

    z=0;
    particle->SetZ(z);

    double cos;
    double dir_vec_x;
    double dir_vec_y;
    double dir_vec_z;

    int excluded =0 ;
    for(int i = 0; i<1e4; i++)
    {
        rnd_x = M.RandomDouble();

        x   =   radius * rnd_x;
        y   =   radius *sqrt(1 - rnd_x*rnd_x);


        particle->SetX(x);
        particle->SetY(y);

        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        dir_vec_x = particle->GetCosPhi()*particle->GetSinTheta();
        dir_vec_y = particle->GetSinPhi()*particle->GetSinTheta();
        dir_vec_z = particle->GetCosTheta();

        //cosine of angle between direction vector and position vector
        cos = (-x * dir_vec_x -y*dir_vec_y - z* dir_vec_z) / radius;

        if(cos < 1 && cos > 0)
            EXPECT_TRUE(A.IsParticleInside(particle));
        else
            EXPECT_FALSE(A.IsParticleInside(particle));

    }

    // Particle is on the top surface.
    // Theta 0° - 90° means particle is moving outside
    // This should be treated as outside
    // Theta 90° - 180° means particle is moving inside (should be treated as inside)
    // The value of phi does not matter
    particle->SetX(0);
    particle->SetY(0);
    particle->SetZ(0.5*height);

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    A.InitCylinder(0,0,0,radius/100,0,height/100);

    for(int i = 0; i<1e4; i++)
    {
        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);
        // Computer precision controll
        if(x*x+y*y -inner_radius*inner_radius ==0 )
        {
            if(theta < 90)
                EXPECT_FALSE(A.IsParticleInside(particle));
            if(theta > 90)
                EXPECT_TRUE(A.IsParticleInside(particle));
        }

    }

    //Make this test for every surface of the box

    // bottom
    particle->SetX(0);
    particle->SetY(0);
    particle->SetZ(-0.5*height);
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

    // Test inner border
    inner_radius    =   5;

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    A.InitCylinder(0,0,0,radius/100,inner_radius/100,height/100);

    z=0;
    particle->SetZ(z);

    excluded=0;

    for(int i = 0; i<1e4; i++)
    {
        rnd_x = M.RandomDouble();

        x   =   inner_radius * rnd_x;
        y   =   inner_radius *sqrt(1 - rnd_x*rnd_x);


        particle->SetX(x);
        particle->SetY(y);

        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        dir_vec_x = particle->GetCosPhi()*particle->GetSinTheta();
        dir_vec_y = particle->GetSinPhi()*particle->GetSinTheta();
        dir_vec_z = particle->GetCosTheta();

        //cosine of angle between direction vector and position vector
        cos = (-x * dir_vec_x -y*dir_vec_y - z* dir_vec_z) / radius;

        if(cos < 1 && cos > 0)
            EXPECT_FALSE(A.IsParticleInside(particle));
        else
            EXPECT_TRUE(A.IsParticleInside(particle));

    }

}

TEST(IsInside , Sphere ) {
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
    double rnd_inner_radius;

    double rnd_theta;
    double rnd_phi;

    double radius       =   10;
    double inner_radius =   0;

    double big_width_x  =   4*radius;
    double big_width_y  =   4*radius;
    double big_height   =   4*radius;

    double x0   =   0;
    double y0   =   0;
    double z0   =   0;

    PROPOSALParticle * particle = new PROPOSALParticle(ParticleType::MuMinus,x,y,z,theta,phi,0,0);
    double is_inside  =0;
    int is_outside =0;

    double volumia_ratio =0;

    MathModel M;
    int number_particles = 1e6;
    int number_volumina  = 1e1;

    Geometry A;

    for(int i = 0; i<number_volumina; i++)
    {
        // Chose the origin of the box-geometry
        // This box should be inside the big box in which the particle
        // will be located
        rnd_x0  = M.RandomDouble();
        rnd_y0  = M.RandomDouble();
        rnd_z0  = M.RandomDouble();

        rnd_inner_radius    = M.RandomDouble();

        x0   =   ( 2 * rnd_x0  - 1)* ( 0.5 * big_width_x - radius );
        y0   =   ( 2 * rnd_y0  - 1)* ( 0.5 * big_width_y - radius );
        z0   =   ( 2 * rnd_z0  - 1)* ( 0.5 * big_height  - radius );

        inner_radius    = radius*rnd_inner_radius;

        // The values are divided by 100 to convert the units...
        // Init functions expects m but here everthing is in cm
        A.InitSphere(x0/100,y0/100,z0/100,radius/100,inner_radius/100);

        volumia_ratio = (4./3.*PI* ( pow(radius ,3)  - pow(inner_radius ,3) ) )
                        /( big_width_x*big_width_y*big_height);

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

            if( sqrt( pow( (x-x0),2 ) + pow( (y-y0),2 ) + pow( (z-z0),2 ) ) < radius &&
                sqrt( pow( (x-x0),2 ) + pow( (y-y0),2 ) + pow( (z-z0),2 ) ) > inner_radius )

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
        ASSERT_NEAR(is_inside ,volumia_ratio*number_particles , 3*sqrt(volumia_ratio*number_particles) );
        is_inside   = 0;
        is_outside  = 0;
    }

    // Test borders

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    A.InitSphere(0,0,0,radius/100,0);

    z=0;
    particle->SetZ(z);

    double cos;
    double dir_vec_x;
    double dir_vec_y;
    double dir_vec_z;

    int excluded =0 ;
    for(int i = 0; i<1e4; i++)
    {
        rnd_x = M.RandomDouble();

        x   =   radius * rnd_x;
        y   =   radius *sqrt(1 - rnd_x*rnd_x);


        particle->SetX(x);
        particle->SetY(y);

        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        dir_vec_x = particle->GetCosPhi()*particle->GetSinTheta();
        dir_vec_y = particle->GetSinPhi()*particle->GetSinTheta();
        dir_vec_z = particle->GetCosTheta();

        //cosine of angle between direction vector and position vector
        cos = (-x * dir_vec_x -y*dir_vec_y - z* dir_vec_z) / radius;

        if(cos < 1 && cos > 0)
            EXPECT_TRUE(A.IsParticleInside(particle));
        else
            EXPECT_FALSE(A.IsParticleInside(particle));

    }


    // Test inner border
    inner_radius    =   5;

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    A.InitSphere(0,0,0,radius/100,inner_radius/100);

    z=0;
    particle->SetZ(z);

    excluded=0;

    for(int i = 0; i<1e4; i++)
    {
        rnd_x = M.RandomDouble();

        x   =   inner_radius * rnd_x;
        y   =   inner_radius *sqrt(1 - rnd_x*rnd_x);


        particle->SetX(x);
        particle->SetY(y);

        rnd_theta   = M.RandomDouble();
        rnd_phi     = M.RandomDouble();

        theta   = rnd_theta*180;
        phi     = rnd_phi*360;
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        dir_vec_x = particle->GetCosPhi()*particle->GetSinTheta();
        dir_vec_y = particle->GetSinPhi()*particle->GetSinTheta();
        dir_vec_z = particle->GetCosTheta();

        //cosine of angle between direction vector and position vector
        cos = (-x * dir_vec_x -y*dir_vec_y - z* dir_vec_z) / radius;

        if(cos < 1 && cos > 0)
            EXPECT_FALSE(A.IsParticleInside(particle));
        else
            EXPECT_TRUE(A.IsParticleInside(particle));

    }
}



TEST(DistanceTo , Sphere ) {
    double x        =   0;
    double y        =   0;
    double z        =   0;
    double theta    =   0;
    double phi      =   0;

    double radius           =   10;
    double inner_radius     =   0;
    double particle_radius  =   0;

    double rnd_phi;
    double rnd_theta;
    double rnd_inner_radius;


    pair<double,double> distance;

    PROPOSALParticle * particle = new PROPOSALParticle(ParticleType::MuMinus,x,y,z,theta,phi,0,0);


    MathModel M;
    int number_particles = 1e5;

    Geometry A;
    cout.precision(16);

    for(int i = 0; i < 11 ; i++)
    {

        particle_radius = 2. + i*2.;

        for(int j = 0; j<number_particles; j++)
        {

            rnd_inner_radius    = M.RandomDouble();
            inner_radius        = radius*rnd_inner_radius;

            // The values are divided by 100 to convert the units...
            // Init functions expects m but here everthing is in cm
            A.InitSphere(0,0,0,radius/100,inner_radius/100);

            rnd_phi             = M.RandomDouble();
            rnd_theta           = M.RandomDouble();

            phi     =   rnd_phi*2*PI;
            theta   =   rnd_theta*PI;


            // Chose particle location and angle

            x   =   -1*particle_radius*cos(phi)*sin(theta);
            y   =   -1*particle_radius*sin(phi)*sin(theta);
            z   =   -1*particle_radius*cos(theta);

            theta   = 180/PI*theta;
            phi     = 180/PI*phi;


            particle->SetX(x);
            particle->SetY(y);
            particle->SetZ(z);
            particle->SetTheta(theta);
            particle->SetPhi(phi);

            distance    =   A.DistanceToBorder(particle);

            if(particle_radius < radius && particle_radius > inner_radius)
            {
                EXPECT_EQ( distance.second, -1.);
                ASSERT_NEAR(distance.first, particle_radius-inner_radius,1e-8*(particle_radius-inner_radius));
            }
            if(particle_radius <= inner_radius)
            {
                ASSERT_NEAR(distance.first, particle_radius+inner_radius,1e-8*(particle_radius+inner_radius));
                ASSERT_NEAR(distance.second, particle_radius+radius,1e-8*(particle_radius+radius));
            }
            if(particle_radius > radius)
            {
                ASSERT_NEAR(distance.first, particle_radius-radius,1e-8*(particle_radius-radius));
                ASSERT_NEAR(distance.second, particle_radius-inner_radius,1e-8*(particle_radius-inner_radius));
            }
            if(particle_radius == radius)
            {
                EXPECT_EQ(distance.second, -1.);
                ASSERT_NEAR(distance.first, particle_radius-inner_radius,1e-8*(particle_radius-inner_radius));
            }

            if(particle_radius >= radius)
            {
                x   = -1*x;
                y   = -1*y;
                z   = -1*z;
                particle->SetX(x);
                particle->SetY(y);
                particle->SetZ(z);

                //Now the particle is moving away from the sphere so we expect no intersection
                distance    =   A.DistanceToBorder(particle);
                EXPECT_EQ(distance.first, -1.);
                EXPECT_EQ(distance.second, -1.);


            }
            if(particle_radius > 20)
            {

                theta   = PI/180*theta;
                phi     = PI/180*phi;

                x   =   -1*inner_radius*cos(phi)*sin(theta);
                y   =   -1*inner_radius*sin(phi)*sin(theta);
                z   =   -1*inner_radius*cos(theta);

                theta   = 180/PI*theta;
                phi     = 180/PI*phi;

                particle->SetX(x);
                particle->SetY(y);
                particle->SetZ(z);
                particle->SetTheta(theta);
                particle->SetPhi(phi);

                distance    =   A.DistanceToBorder(particle);
                ASSERT_NEAR(distance.first, 2*inner_radius,1e-8*(2*inner_radius));
                ASSERT_NEAR(distance.second, inner_radius+radius,1e-8*(inner_radius + radius));

            }
        }
    }
}


TEST(DistanceTo , Cylinder ) {
    double x        =   0;
    double y        =   0;
    double z        =   0;
    double theta    =   0;
    double phi      =   0;

    double height           =   10;
    double radius           =   10;
    double inner_radius     =   0;
    double particle_radius  =   0;

    double rnd_phi;
    double rnd_inner_radius;


    pair<double,double> distance;

    PROPOSALParticle * particle = new PROPOSALParticle(ParticleType::MuMinus,x,y,z,theta,phi,0,0);


    MathModel M;
    int number_particles = 1e5;

    Geometry A;
    cout.precision(16);

    for(int i = 0; i < 10 ; i++)
    {

        particle_radius = 2. + i*2.;

        for(int j = 0; j<number_particles; j++)
        {

            rnd_inner_radius    = M.RandomDouble();
            inner_radius        = radius*rnd_inner_radius;

            // The values are divided by 100 to convert the units...
            // Init functions expects m but here everthing is in cm
            A.InitCylinder(0,0,0,radius/100,inner_radius/100,height/100);

            rnd_phi             = M.RandomDouble();

            phi     =   rnd_phi*2*PI;
            theta   =   0.5*PI;


            // Chose particle location and angle

            x   =   -1*particle_radius*cos(phi)*sin(theta);
            y   =   -1*particle_radius*sin(phi)*sin(theta);
            z   =   -1*0.5*height*cos(theta);

            theta   = 180/PI*theta;
            phi     = 180/PI*phi;

            //if(phi<0)phi=360+phi;

            particle->SetX(x);
            particle->SetY(y);
            particle->SetZ(z);
            particle->SetTheta(theta);
            particle->SetPhi(phi);

            distance    =   A.DistanceToBorder(particle);

            if(particle_radius < radius && particle_radius > inner_radius)
            {
                EXPECT_EQ( distance.second, -1.);
                ASSERT_NEAR(distance.first, particle_radius-inner_radius,1e-8*(particle_radius-inner_radius));
            }
            if(particle_radius <= inner_radius)
            {
                ASSERT_NEAR(distance.first, particle_radius+inner_radius,1e-8*(particle_radius+inner_radius));
                ASSERT_NEAR(distance.second, particle_radius+radius,1e-8*(particle_radius+radius));
            }
            if(particle_radius > radius)
            {
                ASSERT_NEAR(distance.first, particle_radius-radius,1e-8*(particle_radius-radius));
                ASSERT_NEAR(distance.second, particle_radius-inner_radius,1e-8*(particle_radius-inner_radius));
            }
            if(particle_radius == radius)
            {
               // cout<<"TEST "<<distance.first<<"\t"<<distance.second<<"\t"<<inner_radius<<endl;
                EXPECT_EQ(distance.second, -1);
                ASSERT_NEAR(distance.first, particle_radius-inner_radius,1e-8*(particle_radius-inner_radius));
            }

            if(particle_radius >= radius)
            {
                x   = -1*x;
                y   = -1*y;
                z   = -1*z;
                particle->SetX(x);
                particle->SetY(y);
                particle->SetZ(z);

                //Now the particle is moving away from the sphere so we expect no intersection
                distance    =   A.DistanceToBorder(particle);
                EXPECT_EQ(distance.first, -1.);
                EXPECT_EQ(distance.second, -1.);


            }
            if(particle_radius > 20)
            {
                particle_radius =   inner_radius;

                x   =   -1*particle_radius*cos(phi)*sin(theta);
                y   =   -1*particle_radius*sin(phi)*sin(theta);
                z   =   -1*particle_radius*cos(theta);

                particle->SetX(x);
                particle->SetY(y);
                particle->SetZ(z);
                particle->SetTheta(theta);
                particle->SetPhi(phi);

                distance    =   A.DistanceToBorder(particle);
                ASSERT_NEAR(distance.first, 2*inner_radius,1e-8*(2*inner_radius));
                ASSERT_NEAR(distance.second, inner_radius+radius,1e-8*(inner_radius + radius));

            }
        }
    }

    //One test for inner_radius =0
    inner_radius    =   0;

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    A.InitCylinder(0,0,0,radius/100,inner_radius/100,height/100);

    // Chose particle location and angle

    x   =   0;
    y   =   0;
    z   =   height +10;

    particle->SetX(x);
    particle->SetY(y);
    particle->SetZ(z);
    particle->SetTheta(180);
    particle->SetPhi(0);

    distance    =   A.DistanceToBorder(particle);

    ASSERT_NEAR(distance.first,z-0.5*height,1e-8*(z-0.5*height));
    ASSERT_NEAR(distance.second,z+0.5*height,1e-8*(z+0.5*height));

    double rnd_alpha;
    double alpha;

    inner_radius    =   6;

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    A.InitCylinder(0,0,0,radius/100,inner_radius/100,height/100);

    for(int j = 0; j<number_particles; j++)
    {

        rnd_alpha = M.RandomDouble();

        rnd_phi   = M.RandomDouble();

        alpha   =   0.3*PI*rnd_alpha;


        // Chose particle location and angle

        x   =   0;
        y   =   0;
        z   =   height +0.5*height;

        phi     = 360*rnd_phi;
        theta   = 180 - 180/PI*alpha;

        particle->SetX(x);
        particle->SetY(y);
        particle->SetZ(z);
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        double dist1 = inner_radius/sin(alpha);
        double dist2 = radius/sin(alpha);

        distance    =   A.DistanceToBorder(particle);

        //  case 1 throught inner cylinder => no intersection
        //  ___  x    ___
        // |   | |   |   |
        // |   | |   |   |
        // |   | |   |   |
        // |   | |   |   |
        // |   | |   |   |
        // |___| |   |___|
        //       |
        if( alpha < atan(inner_radius/(z+0.5*height)) )
        {
            EXPECT_EQ(distance.first,-1);
            EXPECT_EQ(distance.second,-1);
        }
        //  case 2 first inner cylinder then bottom surface
        //  ___  x      ___
        // |   |  \    |   |
        // |   |   \   |   |
        // |   |    \  |   |
        // |   |     \ |   |
        // |   |      *|   |
        // |___|       |\ _|
        //               *

        else if( alpha < atan(radius/(z+0.5*height)) )
        {
            dist2 = (z+0.5*height)/cos(alpha);
            ASSERT_NEAR(distance.first,dist1,1e-8*(dist1));
            ASSERT_NEAR(distance.second,dist2,1e-8*(dist2));
        }
        //  case 3 first inner cylinder then outer cylinder
        //  ___     x   ___
        // |   |     \ |   |
        // |   |      \|   |
        // |   |       *   |
        // |   |       |\  |
        // |   |       | \ |
        // |   |       |  \|
        // |   |       |   *
        // |___|       |___|\
        //

        else if( alpha < atan(inner_radius/height) )
        {
            ASSERT_NEAR(distance.first,dist1,1e-8*(dist1));
            ASSERT_NEAR(distance.second,dist2,1e-8*(dist2));
        }
        //  case 4 first upper surface then outer cylinder
        //            x
        //             \
        //  ___         \__
        // |   |       | * |
        // |   |       |  \|
        // |   |       |   *
        // |   |       |   |\
        // |   |       |   |
        // |   |       |   |
        // |___|       |___|
        //
        else if(alpha < atan(radius/height))
        {
            dist1 = height/cos(alpha);
            ASSERT_NEAR(distance.first,dist1,1e-8*(dist1));
            ASSERT_NEAR(distance.second,dist2,1e-8*(dist2));
        }
        //  case 5  no intersection
        //      x_____________
        //  ___      ___
        // |   |    |   |
        // |   |    |   |
        // |   |    |   |
        // |   |    |   |
        // |   |    |   |
        // |___|    |___|
        //
        else
        {
            EXPECT_EQ(distance.first,-1);
            EXPECT_EQ(distance.second,-1);
        }
    }
}


TEST(DistanceTo , Box ) {
    double x        =   0;
    double y        =   0;
    double z        =   0;
    double theta    =   0;
    double phi      =   0;

    double width    =   10;
    double height   =   width;

    double rnd_phi;
    double rnd_theta;


    double dist;
    double dist1;
    double dist2;

    pair<double,double> distance;

    PROPOSALParticle * particle = new PROPOSALParticle(ParticleType::MuMinus,x,y,z,theta,phi,0,0);


    MathModel M;
    int number_particles = 1e5;

    Geometry A;
    cout.precision(16);

    for(int j = 0; j<number_particles; j++)
    {
        rnd_phi =   M.RandomDouble();

        // The values are divided by 100 to convert the units...
        // Init functions expects m but here everthing is in cm
        A.InitBox(0,0,0,width/100,width/100,height/100);

        phi     =   rnd_phi*2*PI;
        theta   =   0.5*PI;


        // Chose particle location and angle


        theta   = 180/PI*theta;
        phi     = 180/PI*phi;

        particle->SetX(x);
        particle->SetY(y);
        particle->SetZ(z);
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        distance    =   A.DistanceToBorder(particle);

        if( phi < 45)
        {
            dist = 0.5*width/cos(phi/180*PI);
            EXPECT_EQ( distance.second, -1.);
            ASSERT_NEAR(distance.first, dist,1e-8*(dist));
        }
        else if( phi < 90)
        {
            phi = 90 - phi;
            dist = 0.5*width/cos(phi/180*PI);
            EXPECT_EQ( distance.second, -1.);
            ASSERT_NEAR(distance.first, dist,1e-8*(dist));
        }
        else if( phi < 135)
        {
            phi = phi - 90;
            dist = 0.5*width/cos(phi/180*PI);
            EXPECT_EQ( distance.second, -1.);
            ASSERT_NEAR(distance.first, dist,1e-8*(dist));
        }
        else if( phi < 180)
        {
            phi = 180 - phi;
            dist = 0.5*width/cos(phi/180*PI);
            EXPECT_EQ( distance.second, -1.);
            ASSERT_NEAR(distance.first, dist,1e-8*(dist));
        }
        else if( phi < 225)
        {
            phi = phi - 180;
            dist = 0.5*width/cos(phi/180*PI);
            EXPECT_EQ( distance.second, -1.);
            ASSERT_NEAR(distance.first, dist,1e-8*(dist));
        }
        else if( phi < 270)
        {
            phi = 270 - phi;
            dist = 0.5*width/cos(phi/180*PI);
            EXPECT_EQ( distance.second, -1.);
            ASSERT_NEAR(distance.first, dist,1e-8*(dist));
        }
        else if( phi < 315)
        {
            phi = phi -270;
            dist = 0.5*width/cos(phi/180*PI);
            EXPECT_EQ( distance.second, -1.);
            ASSERT_NEAR(distance.first, dist,1e-8*(dist));
        }
        else if( phi < 360)
        {
            phi = 360 - phi;
            dist = 0.5*width/cos(phi/180*PI);
            EXPECT_EQ( distance.second, -1.);
            ASSERT_NEAR(distance.first, dist,1e-8*(dist));
        }

    }
    for(int i = 0; i < number_particles;i++)
    {
        rnd_phi =   M.RandomDouble();

        phi     =   rnd_phi*0.5*PI;
        theta   =   0.5*PI;


        // Chose particle location and angle

        theta   = 180/PI*theta;
        phi     = 180/PI*phi;

        x= -1*width;
        y=0;
        z=0;

        particle->SetX(x);
        particle->SetY(y);
        particle->SetZ(z);
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        phi =   phi/180*PI;
        distance    =   A.DistanceToBorder(particle);

        //                       ________________           z|
        //                      |               |            |
        //                      |               |            |_____
        //                      |               |                  x
        //     x----------------*---------------*--------->
        //                      |               |
        //                      |               |
        //                      |               |
        //                      |_______________|

        if( phi < atan(0.5*width/(0.5*width-x)) )
        {
            dist1   =   (-x-0.5*width)/cos(phi);
            dist2   =   (0.5*width-x)/cos(phi);
            ASSERT_NEAR(distance.first, dist1,1e-8*(dist1));
            ASSERT_NEAR(distance.second, dist2,1e-8*(dist2));

        }
        //
        //                          ^
        //                         /                        z|
        //                       _*_____________             |
        //                      |/              |            |_____
        //                      *               |                  x
        //                     /|               |
        //                    / |               |
        //                   /  |               |
        //                  x   |               |
        //                      |               |
        //                      |_______________|

        else if( phi < atan(width*0.5/(-x -0.5*width)))
        {
            dist1   =   (-x-0.5*width)/cos(phi);
            dist2   =   0.5*width/sin(phi);
            ASSERT_NEAR(distance.first, dist1,1e-8*(dist1));
            ASSERT_NEAR(distance.second, dist2,1e-8*(dist2));

        }
        //                       ^
        //                      /
        //                     / _______________
        //                    / |               |  z|
        //                   /  |               |   |
        //                  /   |               |   |_____
        //                 /    |               |        x
        //                /     |               |
        //               x      |               |
        //                      |               |
        //                      |_______________|

        else
        {
            EXPECT_EQ(distance.first,-1);
            EXPECT_EQ(distance.second,-1);
        }
    }

    //and one test for z surfaces
    for(int i = 0; i < number_particles;i++)
    {
        rnd_theta =   M.RandomDouble();

        phi     =   0;
        theta   =   rnd_theta*0.5*PI;;


        // Chose particle location and angle

        theta   = 180/PI*theta;
        phi     = 180/PI*phi;

        x=0;
        y=0;
        z=-1*height;

        particle->SetX(x);
        particle->SetY(y);
        particle->SetZ(z);
        particle->SetTheta(theta);
        particle->SetPhi(phi);

        theta =   theta/180*PI;
        distance    =   A.DistanceToBorder(particle);

        //                       ________________       x|
        //                      |               |        |
        //                      |               |        |_____
        //                      |               |              z
        //     x----------------*---------------*--------->
        //                      |               |
        //                      |               |
        //                      |               |
        //                      |_______________|

        if( theta < atan(0.5*height/(0.5*height-z)) )
        {
            dist1   =   (-z-0.5*height)/cos(theta);
            dist2   =   (0.5*height-z)/cos(theta);
            ASSERT_NEAR(distance.first, dist1,1e-8*(dist1));
            ASSERT_NEAR(distance.second, dist2,1e-8*(dist2));

        }
        //
        //                          ^
        //                         /                    x|
        //                       _*_____________         |
        //                      |/              |        |_____
        //                      *               |             z
        //                     /|               |
        //                    / |               |
        //                   /  |               |
        //                  x   |               |
        //                      |               |
        //                      |_______________|

        else if( theta < atan(height*0.5/(-z -0.5*height)))
        {
            dist1   =   (-z-0.5*height)/cos(theta);
            dist2   =   0.5*height/sin(theta);
            ASSERT_NEAR(distance.first, dist1,1e-8*(dist1));
            ASSERT_NEAR(distance.second, dist2,1e-8*(dist2));
        }
        //                       ^
        //                      /
        //                     / _______________        x|
        //                    / |               |        |
        //                   /  |               |        |_____
        //                  /   |               |             z
        //                 /    |               |
        //                /     |               |
        //               x      |               |
        //                      |               |
        //                      |_______________|

        else
        {
            EXPECT_EQ(distance.first,-1);
            EXPECT_EQ(distance.second,-1);
        }
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




