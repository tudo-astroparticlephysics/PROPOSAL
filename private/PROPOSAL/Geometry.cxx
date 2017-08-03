/*
 * Geometry.cxx
 *
 *  Created on: 05.06.2013
 *      Author: koehne
 */

#include <cmath>
// #include <algorithm>
// #include <vector>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Geometry.h"
#include "PROPOSAL/Output.h"

using namespace std;
using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// OStreams
// ------------------------------------------------------------------------- //

namespace PROPOSAL {

ostream& operator<<(ostream& os, Geometry const& geometry)
{
    os << "--------Geometry( " << &geometry << " )--------" << endl;
    os << "\t" << geometry.name_ << endl;
    os << "\tOrigin (x,y,z):\t" << geometry.x0_ << "\t" << geometry.y0_ << "\t"
       << geometry.z0_ << endl;
    os << "\tHirarchy:\t" << geometry.hirarchy_ << endl;

    geometry.print(os);

    return os;
}

} // namespace PROPOSAL

// ------------------------------------------------------------------------- //
// Geometry
// ------------------------------------------------------------------------- //

Geometry::Geometry(std::string name)
    : x0_(0.)
    , y0_(0.)
    , z0_(0.)
    , name_(name)
    , hirarchy_(0)
{
    // Nothing to do here
}

Geometry::Geometry(std::string name, double x, double y, double z)
    : x0_(100.0 * x)
    , y0_(100.0 * y)
    , z0_(100.0 * z)
    , name_(name)
    , hirarchy_(0)
{
    // Nothing to do here
}

Geometry::Geometry(const Geometry& geometry)
    : x0_(geometry.x0_)
    , y0_(geometry.y0_)
    , z0_(geometry.z0_)
    , name_(geometry.name_)
    , hirarchy_(geometry.hirarchy_)
{
    // Nothing to do here
}

// ------------------------------------------------------------------------- //
void Geometry::swap(Geometry& geometry)
{
    using std::swap;

    swap(x0_, geometry.x0_);
    swap(y0_, geometry.y0_);
    swap(z0_, geometry.z0_);
    swap(hirarchy_, geometry.hirarchy_);

    name_.swap(geometry.name_);
}

// ------------------------------------------------------------------------- //
Geometry& Geometry::operator=(const Geometry& geometry)
{
    if (this != &geometry) {
        x0_ = geometry.x0_;
        y0_ = geometry.y0_;
        z0_ = geometry.z0_;
        name_ = geometry.name_;
        hirarchy_ = geometry.hirarchy_;
    }

    return *this;
}

// ------------------------------------------------------------------------- //
bool Geometry::operator==(const Geometry& geometry) const
{
    if (x0_ != geometry.x0_)
        return false;
    else if (y0_ != geometry.y0_)
        return false;
    else if (z0_ != geometry.z0_)
        return false;
    else if (name_.compare(geometry.name_) != 0)
        return false;
    else if (hirarchy_ != geometry.hirarchy_)
        return false;
    else
        return this->compare(geometry);
}

// ------------------------------------------------------------------------- //
bool Geometry::operator!=(const Geometry& geometry) const
{
    return !(*this == geometry);
}

// ------------------------------------------------------------------------- //
// Member functions
// ------------------------------------------------------------------------- //

bool Geometry::IsInside(PROPOSALParticle* particle)
{
    bool is_inside = false;

    pair<double, double> dist = DistanceToBorder(particle);

    if (dist.first > 0 && dist.second < 0) {
        is_inside = true;
    }
    return is_inside;
}

// ------------------------------------------------------------------------- //
bool Geometry::IsInfront(PROPOSALParticle* particle)
{
    bool is_infront = false;

    pair<double, double> dist = DistanceToBorder(particle);

    if (dist.first > 0 && dist.second > 0) {
        is_infront = true;
    }
    return is_infront;
}

// ------------------------------------------------------------------------- //
bool Geometry::IsBehind(PROPOSALParticle* particle)
{
    bool is_behind = false;

    pair<double, double> dist = DistanceToBorder(particle);

    if (dist.first < 0 && dist.second < 0) {
        is_behind = true;
    }
    return is_behind;
}

// ------------------------------------------------------------------------- //
double Geometry::DistanceToClosestApproach(PROPOSALParticle* particle)
{
    double dir_vec_x = particle->GetCosPhi() * particle->GetSinTheta();
    double dir_vec_y = particle->GetSinPhi() * particle->GetSinTheta();
    double dir_vec_z = particle->GetCosTheta();

    double distance = (x0_ - particle->GetX()) * dir_vec_x +
                      (y0_ - particle->GetY()) * dir_vec_y +
                      (z0_ - particle->GetZ()) * dir_vec_z;
    return distance;
}

/******************************************************************************
*                                   Sphere                                    *
******************************************************************************/

Sphere::Sphere()
    : Geometry("Sphere")
    , radius_(0.0)
    , inner_radius_(0.0)
{
    // Do nothing here
}

Sphere::Sphere(double x0,
               double y0,
               double z0,
               double radius,
               double inner_radius)
    : Geometry("Sphere", x0, y0, z0)
    , radius_(100.0 * radius)
    , inner_radius_(100.0 * inner_radius)
{
    if (inner_radius_ > radius_) {
        log_error("Inner radius %f is greater then radius %f (will be swaped)",
                  inner_radius_,
                  radius_);
        std::swap(inner_radius_, radius_);
    }
    if (inner_radius_ == radius_) {
        log_error("Warning: Inner radius %f == radius %f (Volume is 0)",
                  inner_radius_,
                  radius_);
    }
}

Sphere::Sphere(const Sphere& sphere)
    : Geometry(sphere)
    , radius_(sphere.radius_)
    , inner_radius_(sphere.inner_radius_)
{
    // Nothing to do here
}

// ------------------------------------------------------------------------- //
void Sphere::swap(Geometry& geometry)
{
    using std::swap;

    Sphere* sphere = dynamic_cast<Sphere*>(&geometry);
    if (!sphere) {
        log_warn("Cannot swap Sphere!");
        return;
    }

    Geometry::swap(*sphere);

    swap(inner_radius_, sphere->inner_radius_);
    swap(radius_, sphere->radius_);
}

//------------------------------------------------------------------------- //
Sphere& Sphere::operator=(const Geometry& geometry)
{
    if (this != &geometry) {
        const Sphere* sphere = dynamic_cast<const Sphere*>(&geometry);
        if (!sphere) {
            log_warn("Cannot assign Sphere!");
            return *this;
        }

        Sphere tmp(*sphere);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Sphere::compare(const Geometry& geometry) const
{
    const Sphere* sphere = dynamic_cast<const Sphere*>(&geometry);

    if (!sphere)
        return false;
    else if (inner_radius_ != sphere->inner_radius_)
        return false;
    else if (radius_ != sphere->radius_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
void Sphere::print(std::ostream& os) const
{
    os << "\tRadius: " << radius_ << "\tInner radius: " << inner_radius_
       << endl;
    os << "------------------------------------";
}

// ------------------------------------------------------------------------- //
pair<double, double> Sphere::DistanceToBorder(PROPOSALParticle* particle)
{
    // Calculate intersection of particle trajectory and the sphere
    // sphere (x1 + x0)^2 + (x2 + y0)^2 + (x3 + z0)^2 = radius^2
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph
    // *sinth, sinph *sinth , costh)
    // Insert and transform leads to C * t^2 + B * t + A = 0
    // length of direction vector =1 => C = 1
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle
    // trajectory)

    double A, B, t1, t2;
    double dir_vec_x = particle->GetCosPhi() * particle->GetSinTheta();
    double dir_vec_y = particle->GetSinPhi() * particle->GetSinTheta();
    double dir_vec_z = particle->GetCosTheta();

    pair<double, double> distance;

    double determinant;

    A = pow((particle->GetX() - x0_), 2) + pow((particle->GetY() - y0_), 2) +
        pow((particle->GetZ() - z0_), 2) - pow(radius_, 2);

    B = 2 * ((particle->GetX() - x0_) * dir_vec_x +
             (particle->GetY() - y0_) * dir_vec_y +
             (particle->GetZ() - z0_) * dir_vec_z);

    determinant = pow(B / 2, 2) - A;

    if (determinant > 0) // determinant == 0 (boundery point) is ignored
    {
        t1 = -1 * B / 2 + sqrt(determinant);
        t2 = -1 * B / 2 - sqrt(determinant);

        // Computer precision controll
        if (t1 > 0 && t1 < GEOMETRY_PRECISION)
            t1 = 0;
        if (t2 > 0 && t2 < GEOMETRY_PRECISION)
            t2 = 0;

        // (-1/-1) sphere is behind particle or particle is on border but moving
        // outside
        // ( dist_1 / dist_2 ) sphere is infront of the particle
        // ( dist_1 / -1 ) particle is inside the sphere or on border and moving
        // inside
        if (t1 <= 0)
            distance.first = -1;

        else
            distance.first = t1;

        if (t2 <= 0)
            distance.second = -1;

        else
            distance.second = t2;

        // distance.first should be the smaller one
        if (distance.first < 0)
            std::swap(distance.first, distance.second);
        if (distance.first > 0 && distance.second > 0) {
            if (distance.second < distance.first) {
                std::swap(distance.first, distance.second);
            }
        }

    } else // particle trajectory does not have an intersection with the sphere
    {
        distance.first = -1;
        distance.second = -1;
    }

    // No intersection so we don't have to check the distance to the inner
    // sphere
    // if there is any
    if (distance.first < 0 && distance.second < 0)
        return distance;

    // This sqhere might be hollow and we have to check if the inner border is
    // reached before.
    // So we caluculate the intersection with the inner sphere.

    if (inner_radius_ > 0) {
        A = pow((particle->GetX() - x0_), 2) +
            pow((particle->GetY() - y0_), 2) +
            pow((particle->GetZ() - z0_), 2) - pow(inner_radius_, 2);

        B = 2 * ((particle->GetX() - x0_) * dir_vec_x +
                 (particle->GetY() - y0_) * dir_vec_y +
                 (particle->GetZ() - z0_) * dir_vec_z);

        determinant = pow(B / 2, 2) - A;

        if (determinant > 0) // determinant == 0 (boundery point) is ignored
        {
            t1 = -1 * B / 2 + sqrt(determinant);
            t2 = -1 * B / 2 - sqrt(determinant);

            // Computer precision controll
            if (t1 > 0 && t1 < GEOMETRY_PRECISION)
                t1 = 0;
            if (t2 > 0 && t2 < GEOMETRY_PRECISION)
                t2 = 0;

            // Ok we have an intersection with the inner sphere

            // If distance.first and distance.second are positive this means
            // the sphere is infornt of the particle. So the first distance
            // ( intersection with the outer border) does not change
            // but the second distance has to be updated (intersection with the
            // inner border)
            if (distance.first > 0 && distance.second > 0) {
                if (t1 > 0) {
                    if (t1 < distance.second)
                        distance.second = t1;
                }
                if (t2 > 0) {
                    if (t2 < distance.second)
                        distance.second = t2;
                }
            } else // The particle is inside the outer sphere
            {
                // The inner cylinder is infront of the particle trajectory
                // distance.first has to be updated
                if (t1 > 0 && t2 > 0) {
                    if (t1 < t2)
                        distance.first = t1;
                    else
                        distance.first = t2;
                }
                // The particle is inside the inner sphere
                // this means distance.second becomes distanc.first
                // and distance.first beomces distance to intersection with
                // the inner sphere in direction of the particle trajectory
                if ((t1 > 0 && t2 < 0) || (t2 > 0 && t1 < 0)) {
                    std::swap(distance.first, distance.second);
                    if (t1 > 0)
                        distance.first = t1;
                    else
                        distance.first = t2;
                }
                // Now we have to check if the particle is on the border of
                // the inner sphere
                if (t1 == 0) {
                    // The particle is moving into the inner sphere
                    if (t2 > 0) {
                        std::swap(distance.first, distance.second);
                        distance.first = t2;
                    }
                    // if not we don't have to update distance.first
                }
                if (t2 == 0) {
                    // The particle is moving into the inner sphere
                    if (t1 > 0) {
                        std::swap(distance.first, distance.second);
                        distance.first = t1;
                    }
                    // if not we don't have to update distance.first
                }
            }
        }
    }
    // Make a computer precision controll!
    // This is necessary cause due to numerical effects it meight be happen
    // that a particle which is located on a gemoetry border is treated as
    // inside
    // or outside

    if (distance.first < GEOMETRY_PRECISION)
        distance.first = -1;
    if (distance.second < GEOMETRY_PRECISION)
        distance.second = -1;
    if (distance.first < 0)
        std::swap(distance.first, distance.second);

    return distance;
}

/******************************************************************************
*                                    Box                                      *
******************************************************************************/

Box::Box()
    : Geometry("Box")
    , x_(0.0)
    , y_(0.0)
    , z_(0.0)
{
    // Do nothing here
}

Box::Box(double x0, double y0, double z0, double x, double y, double z)
    : Geometry("Box", x0, y0, z0)
    , x_(100.0 * x)
    , y_(100.0 * y)
    , z_(100.0 * z)
{
    // Do nothing here
}

Box::Box(const Box& box)
    : Geometry(box)
    , x_(box.x_)
    , y_(box.y_)
    , z_(box.z_)
{
    // Nothing to do here
}

// ------------------------------------------------------------------------- //
void Box::swap(Geometry& geometry)
{
    using std::swap;

    Box* box = dynamic_cast<Box*>(&geometry);
    if (!box) {
        log_warn("Cannot swap Box!");
        return;
    }

    Geometry::swap(*box);

    swap(x_, box->x_);
    swap(y_, box->y_);
    swap(z_, box->z_);
}

//------------------------------------------------------------------------- //
Box& Box::operator=(const Geometry& geometry)
{
    if (this != &geometry) {
        const Box* box = dynamic_cast<const Box*>(&geometry);
        if (!box) {
            log_warn("Cannot assign Sphere!");
            return *this;
        }

        Box tmp(*box);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Box::compare(const Geometry& geometry) const
{
    const Box* box = dynamic_cast<const Box*>(&geometry);

    if (!box)
        return false;
    else if (x_ != box->x_)
        return false;
    else if (y_ != box->y_)
        return false;
    else if (z_ != box->z_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
void Box::print(std::ostream& os) const
{
    os << "\tWidth_x: " << x_ << "\tWidth_y " << y_ << "\tHeight: " << z_
       << endl;
    os << "------------------------------------";
}

// ------------------------------------------------------------------------- //
pair<double, double> Box::DistanceToBorder(PROPOSALParticle* particle)
{
    // Calculate intersection of particle trajectory and the box
    // Surface of the box is defined by six planes:
    // E1: x1   =   x0_ + 0.5*x
    // E2: x1   =   x0_ - 0.5*x
    // E3: x2   =   y0_ + 0.5*y
    // E4: x2   =   y0_ - 0.5*y
    // E5: x3   =   z0_ + 0.5*z
    // E6: x3   =   z0_ - 0.5*z
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph
    // *sinth, sinph *sinth , costh)
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle
    // trajectory)

    double dir_vec_x = particle->GetCosPhi() * particle->GetSinTheta();
    double dir_vec_y = particle->GetSinPhi() * particle->GetSinTheta();
    double dir_vec_z = particle->GetCosTheta();

    pair<double, double> distance;
    double t;
    double intersection_x;
    double intersection_y;
    double intersection_z;

    vector<double> dist;

    // intersection with E1
    if (dir_vec_x != 0) // if dir_vec == 0 particle trajectory is parallel to E1
    {
        t = (x0_ + 0.5 * x_ - particle->GetX()) / dir_vec_x;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_y = particle->GetY() + t * dir_vec_y;
            intersection_z = particle->GetZ() + t * dir_vec_z;
            if (intersection_y >= (y0_ - 0.5 * y_) &&
                intersection_y <= (y0_ + 0.5 * y_) &&
                intersection_z >= (z0_ - 0.5 * z_) &&
                intersection_z <= (z0_ + 0.5 * z_)) {
                dist.push_back(t);
            }
        }
    }

    // intersection with E2
    if (dir_vec_x != 0) // if dir_vec == 0 particle trajectory is parallel to E2
    {
        t = (x0_ - 0.5 * x_ - particle->GetX()) / dir_vec_x;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_y = particle->GetY() + t * dir_vec_y;
            intersection_z = particle->GetZ() + t * dir_vec_z;
            if (intersection_y >= (y0_ - 0.5 * y_) &&
                intersection_y <= (y0_ + 0.5 * y_) &&
                intersection_z >= (z0_ - 0.5 * z_) &&
                intersection_z <= (z0_ + 0.5 * z_)) {
                dist.push_back(t);
            }
        }
    }

    // intersection with E3
    if (dir_vec_y != 0) // if dir_vec == 0 particle trajectory is parallel to E3
    {
        t = (y0_ + 0.5 * y_ - particle->GetY()) / dir_vec_y;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_x = particle->GetX() + t * dir_vec_x;
            intersection_z = particle->GetZ() + t * dir_vec_z;
            if (intersection_x >= (x0_ - 0.5 * x_) &&
                intersection_x <= (x0_ + 0.5 * x_) &&
                intersection_z >= (z0_ - 0.5 * z_) &&
                intersection_z <= (z0_ + 0.5 * z_)) {
                dist.push_back(t);
            }
        }
    }

    // intersection with E4
    if (dir_vec_y != 0) // if dir_vec == 0 particle trajectory is parallel to E4
    {
        t = (y0_ - 0.5 * y_ - particle->GetY()) / dir_vec_y;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_x = particle->GetX() + t * dir_vec_x;
            intersection_z = particle->GetZ() + t * dir_vec_z;
            if (intersection_x >= (x0_ - 0.5 * x_) &&
                intersection_x <= (x0_ + 0.5 * x_) &&
                intersection_z >= (z0_ - 0.5 * z_) &&
                intersection_z <= (z0_ + 0.5 * z_)) {
                dist.push_back(t);
            }
        }
    }

    // intersection with E5
    if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E5
    {
        t = (z0_ + 0.5 * z_ - particle->GetZ()) / dir_vec_z;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_x = particle->GetX() + t * dir_vec_x;
            intersection_y = particle->GetY() + t * dir_vec_y;
            if (intersection_x >= (x0_ - 0.5 * x_) &&
                intersection_x <= (x0_ + 0.5 * x_) &&
                intersection_y >= (y0_ - 0.5 * y_) &&
                intersection_y <= (y0_ + 0.5 * y_)) {
                dist.push_back(t);
            }
        }
    }

    // intersection with E6
    if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E6
    {
        t = (z0_ - 0.5 * z_ - particle->GetZ()) / dir_vec_z;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_x = particle->GetX() + t * dir_vec_x;
            intersection_y = particle->GetY() + t * dir_vec_y;
            if (intersection_x >= (x0_ - 0.5 * x_) &&
                intersection_x <= (x0_ + 0.5 * x_) &&
                intersection_y >= (y0_ - 0.5 * y_) &&
                intersection_y <= (y0_ + 0.5 * y_)) {
                dist.push_back(t);
            }
        }
    }

    if (dist.size() < 1) // No intersection with the box
    {
        distance.first = -1;
        distance.second = -1;
    } else if (dist.size() == 1) // Particle is inside the box and we have one
                                 // intersection in direction of the particle
                                 // trajectory
    {
        distance.first = dist.at(0);
        distance.second = -1;
    } else if (dist.size() == 2) // Particle is outside and the box is infront
                                 // of the particle trajectory ( two
                                 // intersections).
    {
        distance.first = dist.at(0);
        distance.second = dist.at(1);
        if (distance.second < distance.first) {
            std::swap(distance.first, distance.second);
        }

    } else {
        log_error("This point should nerver be reached... (-1/-1) is returned");

        distance.first = -1;
        distance.second = -1;
    }

    // Make a computer precision controll!
    // This is necessary cause due to numerical effects it meight be happen
    // that a particle which is located on a gemoetry border is treated as
    // inside
    // or outside
    if (distance.first < GEOMETRY_PRECISION)
        distance.first = -1;
    if (distance.second < GEOMETRY_PRECISION)
        distance.second = -1;
    if (distance.first < 0)
        std::swap(distance.first, distance.second);

    return distance;
}

/******************************************************************************
*                                  Cylinder                                   *
******************************************************************************/

Cylinder::Cylinder()
    : Geometry("Cylinder")
    , radius_(0.0)
    , inner_radius_(0.0)
    , z_(0.0)
{
    // Do nothing here
}

Cylinder::Cylinder(double x0,
                   double y0,
                   double z0,
                   double radius,
                   double inner_radius,
                   double z)
    : Geometry("Cylinder", x0, y0, z0)
    , radius_(100 * radius)
    , inner_radius_(100 * inner_radius)
    , z_(100 * z)
{
    if (inner_radius_ > radius_) {
        log_error("Inner radius %f is greater then radius %f (will be swaped)",
                  inner_radius_,
                  radius_);
        std::swap(inner_radius_, radius_);
    }
    if (inner_radius_ == radius_) {
        log_error("Warning: Inner radius %f == radius %f (Volume is 0)",
                  inner_radius_,
                  radius_);
    }
}

Cylinder::Cylinder(const Cylinder& cylinder)
    : Geometry(cylinder)
    , radius_(cylinder.radius_)
    , inner_radius_(cylinder.inner_radius_)
    , z_(cylinder.z_)
{
    // Nothing to do here
}

// ------------------------------------------------------------------------- //
void Cylinder::swap(Geometry& geometry)
{
    using std::swap;

    Cylinder* cylinder = dynamic_cast<Cylinder*>(&geometry);
    if (!cylinder) {
        log_warn("Cannot swap Cylinder!");
        return;
    }

    Geometry::swap(*cylinder);

    swap(inner_radius_, cylinder->inner_radius_);
    swap(radius_, cylinder->radius_);
    swap(z_, cylinder->z_);
}

//------------------------------------------------------------------------- //
Cylinder& Cylinder::operator=(const Geometry& geometry)
{
    if (this != &geometry) {
        const Cylinder* cylinder = dynamic_cast<const Cylinder*>(&geometry);
        if (!cylinder) {
            log_warn("Cannot assign Sphere!");
            return *this;
        }
        Cylinder tmp(*cylinder);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Cylinder::compare(const Geometry& geometry) const
{
    const Cylinder* cylinder = dynamic_cast<const Cylinder*>(&geometry);

    if (!cylinder)
        return false;
    else if (inner_radius_ != cylinder->inner_radius_)
        return false;
    else if (radius_ != cylinder->radius_)
        return false;
    else if (z_ != cylinder->z_)
        return false;
    else
        return true;
}

void Cylinder::print(std::ostream& os) const
{
    os << "\tRadius: " << radius_ << "\tInnner radius: " << inner_radius_
       << " Height: " << z_ << endl;
    os << "------------------------------------";
}

// ------------------------------------------------------------------------- //
pair<double, double> Cylinder::DistanceToBorder(PROPOSALParticle* particle)
{
    // Calculate intersection of particle trajectory and the cylinder
    // cylinder barrel (x1 + x0)^2 + (x2 + y0)^2  = radius^2 [ z0_-0.5*z_ <
    // particle->z <z0_ - 0.5*z_ ]
    // top/bottom surface:
    // E1: x3   =   z0_ + 0.5*z
    // E2: x3   =   z0_ - 0.5*z
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph
    // *sinth, sinph *sinth , costh)
    // Insert and transform leads to C * t^2 + B * t + A = 0
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle
    // trajectory)

    // (-1/-1) cylinder is behind particle or particle is on border but moving
    // outside
    // ( dist_1 / dist_2 ) cylinder is infront of the particle
    // ( dist_1 / -1 ) particle is inside the cylinder or on border and moving
    // inside

    double A, B, C, t1, t2, t;
    double dir_vec_x = particle->GetCosPhi() * particle->GetSinTheta();
    double dir_vec_y = particle->GetSinPhi() * particle->GetSinTheta();
    double dir_vec_z = particle->GetCosTheta();

    double determinant;

    double intersection_x;
    double intersection_y;
    double intersection_z;

    vector<double> dist;

    pair<double, double> distance;

    if (!(dir_vec_x == 0 && dir_vec_y == 0)) // Otherwise the particle
                                             // trajectory is parallel to
                                             // cylinder barrel
    {

        A = pow((particle->GetX() - x0_), 2) +
            pow((particle->GetY() - y0_), 2) - pow(radius_, 2);

        B = 2 * ((particle->GetX() - x0_) * dir_vec_x +
                 (particle->GetY() - y0_) * dir_vec_y);

        C = dir_vec_x * dir_vec_x + dir_vec_y * dir_vec_y;

        B /= C;
        A /= C;

        determinant = pow(B / 2, 2) - A;

        if (determinant > 0) // determinant == 0 (boundery point) is ignored
        {
            t1 = -1 * B / 2 + sqrt(determinant);
            t2 = -1 * B / 2 - sqrt(determinant);

            // Computer precision controll
            if (t1 > 0 && t1 < GEOMETRY_PRECISION)
                t1 = 0;
            if (t2 > 0 && t2 < GEOMETRY_PRECISION)
                t2 = 0;

            if (t1 > 0) {
                intersection_z = particle->GetZ() + t1 * dir_vec_z;
                // is inside the borders
                if (intersection_z > z0_ - 0.5 * z_ &&
                    intersection_z < z0_ + 0.5 * z_) {
                    dist.push_back(t1);
                }
            }

            if (t2 > 0) {
                intersection_z = particle->GetZ() + t2 * dir_vec_z;
                // is inside the borders
                if (intersection_z > z0_ - 0.5 * z_ &&
                    intersection_z < z0_ + 0.5 * z_) {
                    dist.push_back(t2);
                }
            }
        }
    }

    // if we have found already to intersections we don't have to check for
    // intersections
    // with top or bottom surface
    if (dist.size() < 2) {
        // intersection with E1
        if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel
                            // to E1 (Should not happen)
        {
            t = (z0_ + 0.5 * z_ - particle->GetZ()) / dir_vec_z;
            // Computer precision controll
            if (t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            if (t > 0) // Interection is in particle trajectory direction
            {
                intersection_x = particle->GetX() + t * dir_vec_x;
                intersection_y = particle->GetY() + t * dir_vec_y;

                if (sqrt(pow((intersection_x - x0_), 2) +
                         pow((intersection_y - y0_), 2)) <= radius_ &&
                    sqrt(pow((intersection_x - x0_), 2) +
                         pow((intersection_y - y0_), 2)) >= inner_radius_) {
                    dist.push_back(t);
                }
            }
        }

        // intersection with E2
        if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel
                            // to E2 (Should not happen)
        {
            t = (z0_ - 0.5 * z_ - particle->GetZ()) / dir_vec_z;

            // Computer precision controll
            if (t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            if (t > 0) // Interection is in particle trajectory direction
            {
                intersection_x = particle->GetX() + t * dir_vec_x;
                intersection_y = particle->GetY() + t * dir_vec_y;

                if (sqrt(pow((intersection_x - x0_), 2) +
                         pow((intersection_y - y0_), 2)) <= radius_ &&
                    sqrt(pow((intersection_x - x0_), 2) +
                         pow((intersection_y - y0_), 2)) >= inner_radius_) {
                    dist.push_back(t);
                }
            }
        }
    }
    // No intersection with the outer cylinder
    if (dist.size() < 1) {
        distance.first = -1;
        distance.second = -1;
        //    return distance;
    } else if (dist.size() == 1) // particle is inside the cylinder
    {
        distance.first = dist.at(0);
        distance.second = -1;

    } else if (dist.size() == 2) // cylinder is infront of the particle
    {
        distance.first = dist.at(0);
        distance.second = dist.at(1);

        if (distance.second < distance.first) {
            std::swap(distance.first, distance.second);
        }

    } else {
        log_error("This point should never be reached");
    }
    // This cylinder might be hollow and we have to check if the inner border is
    // reached before.
    // So we caluculate the intersection with the inner cylinder.

    if (inner_radius_ > 0) {
        if (!(dir_vec_x == 0 && dir_vec_y == 0)) {

            A = pow((particle->GetX() - x0_), 2) +
                pow((particle->GetY() - y0_), 2) - pow(inner_radius_, 2);

            B = 2 * ((particle->GetX() - x0_) * dir_vec_x +
                     (particle->GetY() - y0_) * dir_vec_y);

            C = dir_vec_x * dir_vec_x + dir_vec_y * dir_vec_y;

            B /= C;
            A /= C;

            determinant = pow(B / 2, 2) - A;

            if (determinant > 0) // determinant == 0 (boundery point) is ignored
            {
                t1 = -1 * B / 2 + sqrt(determinant);
                t2 = -1 * B / 2 - sqrt(determinant);

                // Computer precision controll
                if (t1 > 0 && t1 < GEOMETRY_PRECISION)
                    t1 = 0;
                if (t2 > 0 && t2 < GEOMETRY_PRECISION)
                    t2 = 0;

                // Ok we have a intersection with the inner cylinder

                // If distance.first and distance.second are positive this means
                // the cylinder is infornt of the particle. So the first
                // distance
                // ( intersection with the outer border) does not change
                // but the second distance has to be updated (intersection with
                // the inner border)
                if (distance.first > 0 && distance.second > 0) {
                    if (t1 > 0) {
                        intersection_z = particle->GetZ() + t1 * dir_vec_z;
                        // is inside the borders
                        if (intersection_z > z0_ - 0.5 * z_ &&
                            intersection_z < z0_ + 0.5 * z_) {
                            if (t1 < distance.second)
                                distance.second = t1;
                        }
                    }

                    if (t2 > 0) {
                        intersection_z = particle->GetZ() + t2 * dir_vec_z;
                        // is inside the borders
                        if (intersection_z > z0_ - 0.5 * z_ &&
                            intersection_z < z0_ + 0.5 * z_) {
                            if (t2 < distance.second)
                                distance.second = t2;
                        }
                    }
                }
                // The particle trajectory hits the inner cylinder first
                // (particle flys through the hole and hits inner cylinder
                // barrel first)
                //  ___  ^     ___
                // |   |  \   |   |
                // |   |   \  |   |
                // |   |    \ |   |
                // |   |     \|   |
                // |   |      *   |
                // |   |      |\  |
                // |   |      | x |
                // |   |      |   |
                // |___|      |___|
                else if (distance.first < 0 && distance.second < 0) {
                    if (t1 > 0) {
                        intersection_z = particle->GetZ() + t1 * dir_vec_z;
                        // is inside the borders
                        if (intersection_z > z0_ - 0.5 * z_ &&
                            intersection_z < z0_ + 0.5 * z_) {
                            distance.first = t1;
                        }
                    }
                    if (t2 > 0) {
                        intersection_z = particle->GetZ() + t2 * dir_vec_z;
                        // is inside the borders
                        if (intersection_z > z0_ - 0.5 * z_ &&
                            intersection_z < z0_ + 0.5 * z_) {
                            distance.first = t2;
                        }
                    }
                }
                // The particle is inside the outer cylinder
                // or particle hits bottom or top surface and then hits
                // the inner cylinder barrel
                else {
                    // The inner cylinder is infront of the particle trajectory
                    // distance.first has to be updated
                    if (t1 > 0 && t2 > 0) {
                        //  _____        _____
                        // |     |      |     |
                        // |     |      |     |
                        // |     |      |     |
                        // |     |      |     |
                        // |     |      |     |
                        // | x---*----->|     |
                        // |     |      |     |
                        // |     |      |     |
                        // |_____|      |_____|
                        //
                        if (particle->GetZ() >= z0_ - 0.5 * z_ &&
                            particle->GetZ() <= z0_ + 0.5 * z_ &&
                            sqrt(pow((particle->GetX() - x0_), 2) +
                                 pow((particle->GetY() - y0_), 2)) <=
                              radius_ + GEOMETRY_PRECISION &&
                            sqrt(pow((particle->GetX() - x0_), 2) +
                                 pow((particle->GetY() - y0_), 2)) >=
                              inner_radius_ - GEOMETRY_PRECISION) {
                            if (t1 < distance.first) {
                                intersection_z =
                                  particle->GetZ() + t1 * dir_vec_z;
                                // is inside the borders
                                if (intersection_z > z0_ - 0.5 * z_ &&
                                    intersection_z < z0_ + 0.5 * z_) {
                                    // This case means particle is inside in
                                    // hits
                                    // inner cylinder first
                                    distance.first = t1;
                                }
                            }
                            if (t2 < distance.first) {
                                intersection_z =
                                  particle->GetZ() + t2 * dir_vec_z;
                                // is inside the borders
                                if (intersection_z > z0_ - 0.5 * z_ &&
                                    intersection_z < z0_ + 0.5 * z_) {
                                    // This case means particle is inside in
                                    // hits
                                    // inner cylinder first
                                    distance.first = t2;
                                }
                            }
                        }
                        //               ^
                        //  _____       / _____
                        // |     |     / |     |
                        // |     |    /  |     |
                        // |     |   /   |     |
                        // |     |  /    |     |
                        // |     | /     |     |
                        // |     |/      |     |
                        // |     *       |     |
                        // |    /|       |     |
                        // |___/_|       |_____|
                        //    *
                        //   x
                        else {
                            intersection_z = particle->GetZ() + t1 * dir_vec_z;
                            // is inside the borders
                            if (intersection_z > z0_ - 0.5 * z_ &&
                                intersection_z < z0_ + 0.5 * z_) {
                                // This case means particle is inside in hits
                                // inner cylinder first
                                distance.second = t1;
                            }

                            if (distance.second < 0) {
                                intersection_z =
                                  particle->GetZ() + t2 * dir_vec_z;
                                // is inside the borders
                                if (intersection_z > z0_ - 0.5 * z_ &&
                                    intersection_z < z0_ + 0.5 * z_) {
                                    // This case means particle is inside in
                                    // hits
                                    // inner cylinder first
                                    distance.second = t2;
                                }
                            }
                        }
                    }
                    // The particle is inside the inner cylinder
                    // this means distance.second becomes distanc.first
                    // and distance.first beomces distance to intersection with
                    // the inner cylinder in direction of the particle
                    // trajectory
                    //  _____        _____
                    // |     |      |     |
                    // |     |      |     |
                    // |     |      |     |
                    // |     |      |     |
                    // |     |      |     |
                    // |     |  x---*-----*------->
                    // |     |      |     |
                    // |     |      |     |
                    // |_____|      |_____|
                    //
                    if ((t1 > 0 && t2 < 0) || (t2 > 0 && t1 < 0)) {
                        if (t1 > 0) {
                            intersection_z = particle->GetZ() + t1 * dir_vec_z;
                            // is inside the borders
                            if (intersection_z > z0_ - 0.5 * z_ &&
                                intersection_z < z0_ + 0.5 * z_) {
                                std::swap(distance.first, distance.second);
                                distance.first = t1;
                            }
                        } else {
                            intersection_z = particle->GetZ() + t2 * dir_vec_z;
                            // is inside the borders
                            if (intersection_z > z0_ - 0.5 * z_ &&
                                intersection_z < z0_ + 0.5 * z_) {
                                std::swap(distance.first, distance.second);
                                distance.first = t2;
                            }
                        }
                    }
                    // Now we have to check if the particle is on the border of
                    // the inner sphere
                    if (t1 == 0) {
                        // The particle is moving into the inner cylinder
                        if (t2 > 0) {
                            intersection_z = particle->GetZ() + t2 * dir_vec_z;
                            // is inside the borders
                            if (intersection_z > z0_ - 0.5 * z_ &&
                                intersection_z < z0_ + 0.5 * z_) {
                                std::swap(distance.first, distance.second);
                                distance.first = t2;
                            }
                        }
                        // if not we don't have to update distance.first
                    }
                    if (t2 == 0) {
                        // The particle is moving into the inner sphere
                        if (t1 > 0) {
                            intersection_z = particle->GetZ() + t1 * dir_vec_z;
                            // is inside the borders
                            if (intersection_z > z0_ - 0.5 * z_ &&
                                intersection_z < z0_ + 0.5 * z_) {
                                std::swap(distance.first, distance.second);
                                distance.first = t1;
                            }
                        }
                        // if not we don't have to update distance.first
                    }
                }
            }
        }
    }

    // Make a computer precision controll!
    // This is necessary cause due to numerical effects it meight be happen
    // that a particle which is located on a gemoetry border is treated as
    // inside
    // or outside

    if (distance.first < GEOMETRY_PRECISION)
        distance.first = -1;
    if (distance.second < GEOMETRY_PRECISION)
        distance.second = -1;
    if (distance.first < 0)
        std::swap(distance.first, distance.second);

    return distance;
}
