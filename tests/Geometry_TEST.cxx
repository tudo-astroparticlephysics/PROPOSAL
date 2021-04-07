
#include <iostream>
#include "gtest/gtest.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/Sphere.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/math/Spherical3D.h"

using namespace PROPOSAL;

Cartesian3D posi = Cartesian3D(0,0,0);
double rad_i = 0;
double rad_o = 1e3;
double x_len = 1;
double y_len = 1;
double z_len = 1;

TEST(Comparison, Comparison_equal)
{
    Sphere A(posi, rad_o, rad_i);
    Sphere B(posi, rad_o, rad_i);
    EXPECT_TRUE(A == B);

    Sphere* C = new Sphere(posi, rad_o, rad_i);
    Sphere* D = new Sphere(posi, rad_o, rad_i);
    EXPECT_TRUE(*C == *D);

    delete C;
    delete D;
}

TEST(Comparison, Comparison_not_equal)
{
    Sphere A(posi, rad_o, rad_i);
    Box B(posi, x_len, y_len, z_len);
    EXPECT_TRUE(A != B);

    Sphere* C = new Sphere(posi, rad_o, rad_i);
    Box* D    = new Box(posi, x_len, y_len, z_len);
    EXPECT_TRUE(*C != *D);

    delete C;
    delete D;
}

TEST(Assignment, Copyconstructor)
{
    Sphere A(posi, rad_o, rad_i);
    Sphere B = A;

    EXPECT_TRUE(A == B);

    Geometry* C = new Sphere(posi, rad_o, rad_i);
    Geometry* D = new Box(posi, x_len, y_len, z_len);

    *D = *C;

    EXPECT_FALSE(*C == *D);

    auto E = std::make_shared<Sphere>(Cartesian3D(1.0, 0.0, 0.0), 20.0, 10.0);
    auto F = ::std::make_shared<Sphere>(posi, rad_o, rad_i);

    *F = *E;

    EXPECT_TRUE(*E == *F);
}

TEST(Assignment, Copyconstructor2)
{
    Sphere A(posi, rad_o, rad_i);
    Sphere B(A);

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Operator)
{
    Sphere A(posi, rad_o, rad_i);
    Sphere B(Cartesian3D(0.0, 0.0, 0.0), 2.0, 1.0);

    EXPECT_TRUE(A != B);

    B = A;

    EXPECT_TRUE(A == B);
}

TEST(DistanceToClosestApproach, Method)
{
    Sphere A(posi, rad_o, rad_i);
    Cartesian3D position = Cartesian3D(1, -1, 0);
    Cartesian3D direction = Cartesian3D(0, 1, 0);

    double distance_closest_approach = A.DistanceToClosestApproach(position, direction);

    ASSERT_NEAR(distance_closest_approach, 1., 1e-9);
}

TEST(IsInside, Box)
{
    Cartesian3D particle_position(0, 0, 0);
    Spherical3D particle_direction(0, 0, 0);

    double rnd_x;
    double rnd_y;
    double rnd_z;

    double rnd_x0;
    double rnd_y0;
    double rnd_z0;

    double rnd_theta;
    double rnd_phi;

    double width_x = 10;
    double width_y = 10;
    double height  = 10;

    double big_width_x = 4 * width_x;
    double big_width_y = 4 * width_y;
    double big_height  = 4 * height;

    Cartesian3D position_geometry(0, 0, 0);

    int is_inside  = 0;
    int is_outside = 0;

    double volumia_ratio = 0;

    // MathModel M;
    int number_particles = 1e6;
    int number_volumina  = 1e1;

    for (int i = 0; i < number_volumina; i++)
    {

        // Chose the origin of the box-geometry
        // This box should be inside the big box in which the particle
        // will be located
        rnd_x0 = RandomGenerator::Get().RandomDouble();
        rnd_y0 = RandomGenerator::Get().RandomDouble();
        rnd_z0 = RandomGenerator::Get().RandomDouble();

        position_geometry.SetCoordinates({(2 * rnd_x0 - 1) * 0.5 * (big_width_x - width_x),
                                          (2 * rnd_y0 - 1) * 0.5 * (big_width_y - width_y),
                                          (2 * rnd_z0 - 1) * 0.5 * (big_height - height)});

        Box A(position_geometry, width_x, width_y, height);

        volumia_ratio = width_x * width_y * height / (big_width_x * big_width_y * big_height);
        for (int j = 0; j < number_particles; j++)
        {

            // Chose particle location and angle
            rnd_x = RandomGenerator::Get().RandomDouble();
            rnd_y = RandomGenerator::Get().RandomDouble();
            rnd_z = RandomGenerator::Get().RandomDouble();

            rnd_theta = RandomGenerator::Get().RandomDouble();
            rnd_phi   = RandomGenerator::Get().RandomDouble();

            particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});
            particle_position.SetCoordinates({(2 * rnd_x - 1) * 0.5 * big_width_x,
                                              (2 * rnd_y - 1) * 0.5 * big_width_y,
                                              (2 * rnd_z - 1) * 0.5 * big_height});

            // if this constraints are true the particle is inside the box geometry
            if (particle_position.GetX() > position_geometry.GetX() - 0.5 * width_x &&
                particle_position.GetX() < position_geometry.GetX() + 0.5 * width_x &&
                particle_position.GetY() > position_geometry.GetY() - 0.5 * width_y &&
                particle_position.GetY() < position_geometry.GetY() + 0.5 * width_y &&
                particle_position.GetZ() > position_geometry.GetZ() - 0.5 * height &&
                particle_position.GetZ() < position_geometry.GetZ() + 0.5 * height)
            {
                is_inside++;
                EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
            } else
            {
                is_outside++;
                EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
            }
        }
        ASSERT_NEAR(1. * is_inside, volumia_ratio * number_particles, 3 * sqrt(volumia_ratio * number_particles));
        is_inside  = 0;
        is_outside = 0;
    }
    // Check what happens if particles are on the border of the box

    Box A(Cartesian3D(0, 0, 0), width_x, width_y, height);

    // Particle is on the top surface.
    // Theta 0° - 90° means particle is moving outside
    // This should be treated as outside
    // Theta 90° - 180° means particle is moving inside (should be treated as inside)
    // The value of phi does not matter
    particle_position.SetCoordinates({0, 0, 0.5 * height});
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

        if (particle_direction.GetZenith() < 0.5 * PI) {
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        }
        if (particle_direction.GetZenith() > 0.5 * PI) {
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
        }
    }

    // Make this test for every surface of the box

    // bottom
    particle_position.SetCoordinates({0, 0, -0.5 * height});
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

        if (particle_direction.GetZenith() > 0.5 * PI) {
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        }
        if (particle_direction.GetZenith() < 0.5 * PI) {
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
        }
    }

    // Surface in positiv x direction
    particle_position.SetCoordinates({0.5 * width_x, 0, 0});
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

        // phi = 0 is in positive x direction
        if (particle_direction.GetAzimuth() < 0.5 * PI ||
            particle_direction.GetAzimuth() > 1.5 * PI)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
    }
    // Surface in negativ x direction
    particle_position.SetCoordinates({-0.5 * width_x, 0, 0});
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

        // phi = 0 is in positive x direction
        if (particle_direction.GetAzimuth() < 0.5 * PI ||
            particle_direction.GetAzimuth() > 1.5 * PI)
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
    }
    // Surface in positiv y direction
    particle_position.SetCoordinates({0, 0.5 * width_y, 0});
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

        // phi = 0 is in positive x direction
        if (particle_direction.GetAzimuth() < PI)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
    }
    // Surface in negativ y direction
    particle_position.SetCoordinates({0, -0.5 * width_y, 0});
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

        // phi = 0 is in positive x direction
        if (particle_direction.GetAzimuth() < PI)
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
    }

    // For completness check one corner
    particle_position.SetCoordinates({0.5 * width_x, 0.5 * width_y, 0.5 * height});
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

        if (particle_direction.GetZenith() < 0.5 * PI ||
            particle_direction.GetAzimuth() < PI ||
            particle_direction.GetAzimuth() > 1.5 * PI)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
    }
}

TEST(IsInside, Cylinder)
{

    Cartesian3D particle_position(0, 0, 0);
    Spherical3D particle_direction(0, 0, 0);

    double rnd_x;
    double rnd_y;
    double rnd_z;

    double rnd_x0;
    double rnd_y0;
    double rnd_z0;

    double rnd_theta;
    double rnd_phi;
    double rnd_inner_radius;

    double radius       = 10;
    double inner_radius = 0;
    double height       = 10;

    double big_width_x = 4 * radius;
    double big_width_y = 4 * radius;
    double big_height  = 4 * height;

    Cartesian3D position_geometry(0, 0, 0);

    int is_inside  = 0;
    int is_outside = 0;

    double volumia_ratio = 0;

    // MathModel M;
    int number_particles = 1e6;
    int number_volumina  = 1e1;

    for (int i = 0; i < number_volumina; i++)
    {

        // Chose the origin of the cylinder-geometry
        // This cylinder should be inside the big box in which the particle
        // will be located
        rnd_x0 = RandomGenerator::Get().RandomDouble();
        rnd_y0 = RandomGenerator::Get().RandomDouble();
        rnd_z0 = RandomGenerator::Get().RandomDouble();

        position_geometry.SetCoordinates({(2 * rnd_x0 - 1) * (0.5 * big_width_x - radius),
                                          (2 * rnd_y0 - 1) * (0.5 * big_width_y - radius),
                                          (2 * rnd_z0 - 1) * 0.5 * (big_height - height)});

        // position_geometry.SetCartesianCoordinates(0,0,0);

        rnd_inner_radius = RandomGenerator::Get().RandomDouble();

        inner_radius = radius * rnd_inner_radius;

        Cylinder A(position_geometry, height, radius, inner_radius);

        volumia_ratio =
            height * PI * (pow(radius, 2) - pow(inner_radius, 2)) / (big_width_x * big_width_y * big_height);
        for (int j = 0; j < number_particles; j++)
        {

            // Chose particle location and angle
            rnd_x = RandomGenerator::Get().RandomDouble();
            rnd_y = RandomGenerator::Get().RandomDouble();
            rnd_z = RandomGenerator::Get().RandomDouble();

            rnd_theta = RandomGenerator::Get().RandomDouble();
            rnd_phi   = RandomGenerator::Get().RandomDouble();

            particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

            particle_position.SetCoordinates({(2 * rnd_x - 1) * 0.5 * big_width_x,
                                              (2 * rnd_y - 1) * 0.5 * big_width_y,
                                              (2 * rnd_z - 1) * 0.5 * big_height});

            // if this constraints are true the particle is inside the cylinder geometry
            if (sqrt(pow((particle_position.GetX() - position_geometry.GetX()), 2) +
                     pow((particle_position.GetY() - position_geometry.GetY()), 2)) < radius &&
                sqrt(pow((particle_position.GetX() - position_geometry.GetX()), 2) +
                     pow((particle_position.GetY() - position_geometry.GetY()), 2)) > inner_radius &&
                particle_position.GetZ() > position_geometry.GetZ() - 0.5 * height &&
                particle_position.GetZ() < position_geometry.GetZ() + 0.5 * height)
            {

                is_inside++;

                EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
            } else
            {
                is_outside++;
                EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
            }
        }
        ASSERT_NEAR(1. * is_inside, volumia_ratio * number_particles, 3 * sqrt(volumia_ratio * number_particles));
        is_inside  = 0;
        is_outside = 0;
    }

    // Test borders

    Sphere B(Cartesian3D(0, 0, 0), radius, 0);

    double cos;

    for (int i = 0; i < 1e4; i++)
    {
        rnd_x = RandomGenerator::Get().RandomDouble();

        particle_position.SetCoordinates({radius * rnd_x, radius * sqrt(1 - rnd_x * rnd_x), 0});

        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

        // cosine of angle between direction vector and position vector
        cos = -(particle_position * particle_direction) / radius;

        if (cos < 1 && cos > 0)
            EXPECT_TRUE(B.IsInside(particle_position, particle_direction));
        else
            EXPECT_FALSE(B.IsInside(particle_position, particle_direction));
    }

    // Particle is on the top surface.
    // Theta 0° - 90° means particle is moving outside
    // This should be treated as outside
    // Theta 90° - 180° means particle is moving inside (should be treated as inside)
    // The value of phi does not matter
    particle_position.SetCoordinates({0, 0, 0.5 * height});

    Cylinder C(Cartesian3D(0, 0, 0), height, radius, 0);

    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});
        // Computer precision controll
        if (particle_position.GetX() * particle_position.GetX() + particle_position.GetY() * particle_position.GetY() -
                inner_radius * inner_radius ==
            0)
        {
            if (particle_direction.GetZenith() < PI / 2.) {
                EXPECT_FALSE(C.IsInside(particle_position, particle_direction));
            }
            if (particle_direction.GetZenith() > PI / 2.) {
                EXPECT_TRUE(C.IsInside(particle_position, particle_direction));
            }
        }
    }

    // Make this test for every surface of the box

    // bottom
    particle_position.SetCoordinates({0, 0, -0.5 * height});
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

        if (particle_direction.GetZenith() > PI / 2.) {
            EXPECT_FALSE(C.IsInside(particle_position, particle_direction));
        }
        if (particle_direction.GetZenith() < PI / 2.) {
            EXPECT_TRUE(C.IsInside(particle_position, particle_direction));
        }
    }

    // Test inner border
    inner_radius = 5;

    Cylinder A(Cartesian3D(0, 0, 0), height, radius, inner_radius);

    for (int i = 0; i < 1e4; i++)
    {
        rnd_x = RandomGenerator::Get().RandomDouble();

        particle_position.SetCoordinates({inner_radius * rnd_x, inner_radius * sqrt(1 - rnd_x * rnd_x), 0});

        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates(Spherical3D(1, rnd_phi * 2 * PI, rnd_theta * PI).GetCartesianCoordinates());

        // cosine of angle between direction vector and position vector
        cos = -(particle_position * particle_direction) / radius;

        if (cos < 1 && cos > 0)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
    }
}

TEST(IsInside, Sphere)
{

    Cartesian3D particle_position(0, 0, 0);
    Cartesian3D particle_direction(0, 0, 0);

    double rnd_x;
    double rnd_y;
    double rnd_z;

    double rnd_x0;
    double rnd_y0;
    double rnd_z0;
    double rnd_inner_radius;

    double rnd_theta;
    double rnd_phi;

    double radius       = 10;
    double inner_radius = 0;

    double big_width_x = 4 * radius;
    double big_width_y = 4 * radius;
    double big_height  = 4 * radius;

    Cartesian3D position_geometry(0, 0, 0);

    int is_inside  = 0;
    int is_outside = 0;

    double volumia_ratio = 0;

    // MathModel M;
    int number_particles = 1e6;
    int number_volumina  = 1e1;

    for (int i = 0; i < number_volumina; i++)
    {
        // Chose the origin of the box-geometry
        // This box should be inside the big box in which the particle
        // will be located
        rnd_x0 = RandomGenerator::Get().RandomDouble();
        rnd_y0 = RandomGenerator::Get().RandomDouble();
        rnd_z0 = RandomGenerator::Get().RandomDouble();

        rnd_inner_radius = RandomGenerator::Get().RandomDouble();

        position_geometry.SetCoordinates({(2 * rnd_x0 - 1) * (0.5 * big_width_x - radius),
                                          (2 * rnd_y0 - 1) * (0.5 * big_width_y - radius),
                                          (2 * rnd_z0 - 1) * (0.5 * big_height - radius)});

        inner_radius = radius * rnd_inner_radius;

        Sphere A(position_geometry, radius, inner_radius);

        volumia_ratio =
            (4. / 3. * PI * (pow(radius, 3) - pow(inner_radius, 3))) / (big_width_x * big_width_y * big_height);

        for (int j = 0; j < number_particles; j++)
        {

            // Chose particle location and angle
            rnd_x = RandomGenerator::Get().RandomDouble();
            rnd_y = RandomGenerator::Get().RandomDouble();
            rnd_z = RandomGenerator::Get().RandomDouble();

            rnd_theta = RandomGenerator::Get().RandomDouble();
            rnd_phi   = RandomGenerator::Get().RandomDouble();

            particle_direction.SetCoordinates(Spherical3D(1, rnd_phi * 2 * PI, rnd_theta * PI).GetCartesianCoordinates());

            particle_position.SetCoordinates({(2 * rnd_x - 1) * 0.5 * big_width_x,
                                              (2 * rnd_y - 1) * 0.5 * big_width_y,
                                              (2 * rnd_z - 1) * 0.5 * big_height});

            if ((particle_position - position_geometry).magnitude() < radius &&
                (particle_position - position_geometry).magnitude() > inner_radius)
            {
                is_inside++;
                EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
            } else
            {
                is_outside++;
                EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
            }
        }
        ASSERT_NEAR(1. * is_inside, volumia_ratio * number_particles, 3 * sqrt(volumia_ratio * number_particles));
        is_inside  = 0;
        is_outside = 0;
    }

    // Test borders

    Sphere A(Cartesian3D(0, 0, 0), radius);

    double cos;

    for (int i = 0; i < 1e4; i++)
    {
        rnd_x = RandomGenerator::Get().RandomDouble();

        particle_position.SetCoordinates({radius * rnd_x, radius * sqrt(1 - rnd_x * rnd_x), 0});

        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates(Spherical3D(1, rnd_phi * 2 * PI, rnd_theta * PI).GetCartesianCoordinates());

        // cosine of angle between direction vector and position vector
        cos = -(particle_position * particle_direction) / radius;

        if (cos < 1 && cos > 0)
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
    }

    // Test inner border
    inner_radius = 5;

    Sphere B(Cartesian3D(0, 0, 0), radius, inner_radius);

    for (int i = 0; i < 1e4; i++)
    {
        rnd_x = RandomGenerator::Get().RandomDouble();

        particle_position.SetCoordinates({inner_radius * rnd_x, inner_radius * sqrt(1 - rnd_x * rnd_x), 0});

        rnd_theta = RandomGenerator::Get().RandomDouble();
        rnd_phi   = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates(Spherical3D(1, rnd_phi * 2 * PI, rnd_theta * PI).GetCartesianCoordinates());

        // cosine of angle between direction vector and position vector
        cos = -(particle_position * particle_direction) / radius;

        if (cos < 1 && cos > 0)
            EXPECT_FALSE(B.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(B.IsInside(particle_position, particle_direction));
    }
}

TEST(DistanceTo, Sphere)
{
    Cartesian3D particle_position(0, 0, 0);
    Spherical3D particle_direction(0, 0, 0);

    double radius          = 10;
    double inner_radius    = 0;
    double particle_radius = 0;

    double rnd_phi;
    double rnd_theta;
    double rnd_inner_radius;

    std::pair<double, double> distance;

    // MathModel M;
    int number_particles = 1e5;

    std::cout.precision(16);

    for (int i = 0; i < 11; i++)
    {

        particle_radius = 2. + i * 2.;

        for (int j = 0; j < number_particles; j++)
        {

            rnd_inner_radius = RandomGenerator::Get().RandomDouble();
            inner_radius     = radius * rnd_inner_radius;

            Sphere A(Cartesian3D(0, 0, 0), radius, inner_radius);

            rnd_phi   = RandomGenerator::Get().RandomDouble();
            rnd_theta = RandomGenerator::Get().RandomDouble();

            particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, rnd_theta * PI});

            // Chose particle location and angle
            particle_position.SetCoordinates(Spherical3D(particle_radius, rnd_phi * 2 * PI, rnd_theta * PI).GetCartesianCoordinates());
            particle_position = -particle_position;

            distance = A.DistanceToBorder(particle_position, particle_direction);

            if (particle_radius < radius && particle_radius > inner_radius)
            {
                EXPECT_EQ(distance.second, -1.);
                ASSERT_NEAR(distance.first, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }
            if (particle_radius <= inner_radius)
            {
                ASSERT_NEAR(distance.first, particle_radius + inner_radius, 1e-8 * (particle_radius + inner_radius));
                ASSERT_NEAR(distance.second, particle_radius + radius, 1e-8 * (particle_radius + radius));
            }
            if (particle_radius > radius)
            {
                ASSERT_NEAR(distance.first, particle_radius - radius, 1e-8 * (particle_radius - radius));
                ASSERT_NEAR(distance.second, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }
            if (particle_radius == radius)
            {
                EXPECT_EQ(distance.second, -1.);
                ASSERT_NEAR(distance.first, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }

            if (particle_radius >= radius)
            {
                particle_position = -1 * particle_position;

                // Now the particle is moving away from the sphere so we expect no intersection
                distance = A.DistanceToBorder(particle_position, particle_direction);
                EXPECT_EQ(distance.first, -1.);
                EXPECT_EQ(distance.second, -1.);
            }
            if (particle_radius > 20)
            {
                particle_position.SetCoordinates(Spherical3D(inner_radius, rnd_phi * 2 * PI, rnd_theta * PI).GetCartesianCoordinates());
                particle_position = -particle_position;

                distance = A.DistanceToBorder(particle_position, particle_direction);
                ASSERT_NEAR(distance.first, 2 * inner_radius, 1e-8 * (2 * inner_radius));
                ASSERT_NEAR(distance.second, inner_radius + radius, 1e-8 * (inner_radius + radius));
            }
        }
    }
}

//
TEST(DistanceTo, Cylinder)
{
    Cartesian3D particle_position(0, 0, 0);
    Spherical3D particle_direction(0, 0, 0);

    double height          = 10;
    double radius          = 10;
    double inner_radius    = 0;
    double particle_radius = 0;

    double rnd_phi;
    double rnd_inner_radius;

    double z;

    std::pair<double, double> distance;

    // MathModel M;
    int number_particles = 1e5;

    std::cout.precision(16);

    for (int i = 0; i < 10; i++)
    {

        particle_radius = 2. + i * 2.;

        for (int j = 0; j < number_particles; j++)
        {

            rnd_inner_radius = RandomGenerator::Get().RandomDouble();
            inner_radius     = radius * rnd_inner_radius;

            Cylinder A(Cartesian3D(0, 0, 0), height, radius, inner_radius);

            rnd_phi = RandomGenerator::Get().RandomDouble();

            // Chose particle location and angle
            particle_position.SetCoordinates(Spherical3D(1, rnd_phi * 2 * PI, 0.5 * PI).GetCartesianCoordinates());
            particle_position.SetCoordinates({particle_radius * particle_position.GetX(),
                                              particle_radius * particle_position.GetY(),
                                              0.5 * height * particle_position.GetZ()});
            particle_position = -particle_position;

            particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, 0.5 * PI});

            distance = A.DistanceToBorder(particle_position, particle_direction);

            if (particle_radius < radius && particle_radius > inner_radius)
            {
                EXPECT_EQ(distance.second, -1.);
                ASSERT_NEAR(distance.first, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }
            if (particle_radius <= inner_radius)
            {
                ASSERT_NEAR(distance.first, particle_radius + inner_radius, 1e-8 * (particle_radius + inner_radius));
                ASSERT_NEAR(distance.second, particle_radius + radius, 1e-8 * (particle_radius + radius));
            }
            if (particle_radius > radius)
            {
                ASSERT_NEAR(distance.first, particle_radius - radius, 1e-8 * (particle_radius - radius));
                ASSERT_NEAR(distance.second, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }
            if (particle_radius == radius)
            {
                EXPECT_EQ(distance.second, -1);
                ASSERT_NEAR(distance.first, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }

            if (particle_radius >= radius)
            {
                particle_position = -particle_position;

                // Now the particle is moving away from the sphere so we expect no intersection
                distance = A.DistanceToBorder(particle_position, particle_direction);
                EXPECT_EQ(distance.first, -1.);
                EXPECT_EQ(distance.second, -1.);
            }
            if (particle_radius > 20)
            {
                particle_radius = inner_radius;

                particle_position.SetCoordinates(Spherical3D(particle_radius, rnd_phi * 2 * PI, 0.5 * PI).GetCartesianCoordinates());
                particle_position = -particle_position;

                distance = A.DistanceToBorder(particle_position, particle_direction);
                ASSERT_NEAR(distance.first, 2 * inner_radius, 1e-8 * (2 * inner_radius));
                ASSERT_NEAR(distance.second, inner_radius + radius, 1e-8 * (inner_radius + radius));
            }
        }
    }

    // One test for inner_radius =0
    inner_radius = 0;

    Cylinder B(Cartesian3D(0, 0, 0), height, radius, inner_radius);

    // Chose particle location and angle

    particle_position.SetCoordinates({0, 0, height + 10});

    particle_direction.SetCoordinates({1, 0, PI});

    z = particle_position.GetZ();

    distance = B.DistanceToBorder(particle_position, particle_direction);

    ASSERT_NEAR(distance.first, z - 0.5 * height, 1e-8 * (z - 0.5 * height));
    ASSERT_NEAR(distance.second, z + 0.5 * height, 1e-8 * (z + 0.5 * height));

    double rnd_alpha;
    double alpha;

    inner_radius = 6;

    Cylinder A(Cartesian3D(0, 0, 0), height, radius, inner_radius);

    for (int j = 0; j < number_particles; j++)
    {

        rnd_alpha = RandomGenerator::Get().RandomDouble();

        rnd_phi = RandomGenerator::Get().RandomDouble();

        alpha = 0.3 * PI * rnd_alpha;

        // Chose particle location and angle

        particle_position.SetCoordinates({0, 0, height + 0.5 * height});
        z = particle_position.GetZ();

        particle_direction.SetCoordinates({1, rnd_phi, PI - alpha});

        double dist1 = inner_radius / std::sin(alpha);
        double dist2 = radius / std::sin(alpha);

        distance = A.DistanceToBorder(particle_position, particle_direction);

        /*  case 1 throught inner cylinder => no intersection
            ___  x    ___
           |   | |   |   |
           |   | |   |   |
           |   | |   |   |
           |   | |   |   |
           |   | |   |   |
           |___| |   |___|
                 |
        */
        if (alpha < std::atan(inner_radius / (z + 0.5 * height)))
        {
            EXPECT_EQ(distance.first, -1);
            EXPECT_EQ(distance.second, -1);
        }
        /*  case 2 first inner cylinder then bottom surface
            ___  x      ___
           |   |  \    |   |
           |   |   \   |   |
           |   |    \  |   |
           |   |     \ |   |
           |   |      *|   |
           |___|       |\ _|
                         *
        */

        else if (alpha < std::atan(radius / (z + 0.5 * height)))
        {
            dist2 = (z + 0.5 * height) / std::cos(alpha);
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));
        }
        /*  case 3 first inner cylinder then outer cylinder
            ___     x   ___
           |   |     \ |   |
           |   |      \|   |
           |   |       *   |
           |   |       |\  |
           |   |       | \ |
           |   |       |  \|
           |   |       |   *
           |___|       |___|\
        */

        else if (alpha < std::atan(inner_radius / height))
        {
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));
        }
        /*  case 4 first upper surface then outer cylinder
                      x
                       \
            ___         \__
           |   |       | * |
           |   |       |  \|
           |   |       |   *
           |   |       |   |\
           |   |       |   |
           |   |       |   |
           |___|       |___|
        */
        else if (alpha < std::atan(radius / height))
        {
            dist1 = height / std::cos(alpha);
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));
        }
        /*  case 5  no intersection
                x_____________
            ___      ___
           |   |    |   |
           |   |    |   |
           |   |    |   |
           |   |    |   |
           |   |    |   |
           |___|    |___|
        */
        else
        {
            EXPECT_EQ(distance.first, -1);
            EXPECT_EQ(distance.second, -1);
        }
    }
}

TEST(DistanceTo, Box)
{
    Cartesian3D particle_position(0, 0, 0);
    Spherical3D particle_direction(0, 0, 0);

    double width  = 10;
    double height = width;

    double rnd_phi;
    double rnd_theta;

    double phi;

    double dist;
    double dist1;
    double dist2;

    std::pair<double, double> distance;

    // MathModel M;
    int number_particles = 1e5;

    std::cout.precision(16);

    for (int j = 0; j < number_particles; j++)
    {
        rnd_phi = RandomGenerator::Get().RandomDouble();

        Box A(Cartesian3D(0, 0, 0), width, width, height);

        particle_direction.SetCoordinates({1, rnd_phi * 2 * PI, 0.5 * PI});

        distance = A.DistanceToBorder(particle_position, particle_direction);

        phi = particle_direction.GetAzimuth() * 180. / PI;
        if (phi < 45)
        {
            dist = 0.5 * width / std::cos(phi / 180 * PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 90)
        {
            phi  = 90 - phi;
            dist = 0.5 * width / std::cos(phi / 180 * PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 135)
        {
            phi  = phi - 90;
            dist = 0.5 * width / std::cos(phi / 180 * PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 180)
        {
            phi  = 180 - phi;
            dist = 0.5 * width / std::cos(phi / 180 * PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 225)
        {
            phi  = phi - 180;
            dist = 0.5 * width / std::cos(phi / 180 * PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 270)
        {
            phi  = 270 - phi;
            dist = 0.5 * width / std::cos(phi / 180 * PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 315)
        {
            phi  = phi - 270;
            dist = 0.5 * width / std::cos(phi / 180 * PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 360)
        {
            phi  = 360 - phi;
            dist = 0.5 * width / std::cos(phi / 180 * PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        }
    }

    Box A(Cartesian3D(0, 0, 0), width, width, height);

    for (int i = 0; i < number_particles; i++)
    {
        rnd_phi = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, rnd_phi * 0.5 * PI, 0.5 * PI});
        particle_position.SetCoordinates({-1 * width, 0, 0});

        phi      = particle_direction.GetAzimuth();
        distance = A.DistanceToBorder(particle_position, particle_direction);

        //                       ________________           z|
        //                      |               |            |
        //                      |               |            |_____
        //                      |               |                  x
        //     x----------------*---------------*--------->
        //                      |               |
        //                      |               |
        //                      |               |
        //                      |_______________|

        if (phi < std::atan(0.5 * width / (0.5 * width - particle_position.GetX())))
        {
            dist1 = (-particle_position.GetX() - 0.5 * width) / std::cos(phi);
            dist2 = (0.5 * width - particle_position.GetX()) / std::cos(phi);
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));

        }
        /*
                                   ^
                                  /                        z|
                                _*_____________             |
                               |/              |            |_____
                               *               |                  x
                              /|               |
                             / |               |
                            /  |               |
                           x   |               |
                               |               |
                               |_______________|
        */

        else if (phi < std::atan(width * 0.5 / (-particle_position.GetX() - 0.5 * width)))
        {
            dist1 = (-particle_position.GetX() - 0.5 * width) / std::cos(phi);
            dist2 = 0.5 * width / std::sin(phi);
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));

        }
        /*                       ^
                                /
                               / _______________
                              / |               |  z|
                             /  |               |   |
                            /   |               |   |_____
                           /    |               |        x
                          /     |               |
                         x      |               |
                                |               |
                                |_______________|
        */

        else
        {
            EXPECT_EQ(distance.first, -1);
            EXPECT_EQ(distance.second, -1);
        }
    }

    // and one test for z surfaces
    for (int i = 0; i < number_particles; i++)
    {
        rnd_theta = RandomGenerator::Get().RandomDouble();

        particle_direction.SetCoordinates({1, 0, rnd_theta * 0.5 * PI});

        particle_position.SetCoordinates({0, 0, -1 * height});

        distance = A.DistanceToBorder(particle_position, particle_direction);

        /*                       ________________       x|
                                |               |        |
                                |               |        |_____
                                |               |              z
               x----------------*---------------*--------->
                                |               |
                                |               |
                                |               |
                                |_______________|
        */

        if (particle_direction.GetZenith() < std::atan(0.5 * height / (0.5 * height - particle_position.GetZ())))
        {
            dist1 = (-particle_position.GetZ() - 0.5 * height) / std::cos(particle_direction.GetZenith());
            dist2 = (0.5 * height - particle_position.GetZ()) / std::cos(particle_direction.GetZenith());
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));

        }
        /*
                                    ^
                                   /                    x|
                                 _*_____________         |
                                |/              |        |_____
                                *               |             z
                               /|               |
                              / |               |
                             /  |               |
                            x   |               |
                                |               |
                                |_______________|
        */

        else if (particle_direction.GetZenith() < std::atan(height * 0.5 / (-particle_position.GetZ() - 0.5 * height)))
        {
            dist1 = (-particle_position.GetZ() - 0.5 * height) / std::cos(particle_direction.GetZenith());
            dist2 = 0.5 * height / std::sin(particle_direction.GetZenith());
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));
        }
        /*                       ^
                                /
                               / _______________        x|
                              / |               |        |
                             /  |               |        |_____
                            /   |               |             z
                           /    |               |
                          /     |               |
                         x      |               |
                                |               |
                                |_______________|
        */

        else
        {
            EXPECT_EQ(distance.first, -1);
            EXPECT_EQ(distance.second, -1);
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
