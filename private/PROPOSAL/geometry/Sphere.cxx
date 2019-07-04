#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/geometry/Sphere.h"

using namespace PROPOSAL;

Sphere::Sphere()
    : Geometry("Sphere")
    , radius_(0.0)
    , inner_radius_(0.0)
{
    // Do nothing here
}

Sphere::Sphere(const Vector3D position, double radius, double inner_radius)
    : Geometry("Sphere", position)
    , radius_(100.0 * radius)
    , inner_radius_(100.0 * inner_radius)
{
    if (inner_radius_ > radius_)
    {
        log_error("Inner radius %f is greater then radius %f (will be swaped)", inner_radius_, radius_);
        std::swap(inner_radius_, radius_);
    }
    if (inner_radius_ == radius_)
    {
        log_error("Warning: Inner radius %f == radius %f (Volume is 0)", inner_radius_, radius_);
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
    Sphere* sphere = dynamic_cast<Sphere*>(&geometry);
    if (!sphere)
    {
        log_warn("Cannot swap Sphere!");
        return;
    }

    Geometry::swap(*sphere);

    std::swap(inner_radius_, sphere->inner_radius_);
    std::swap(radius_, sphere->radius_);
}

//------------------------------------------------------------------------- //
Sphere& Sphere::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const Sphere* sphere = dynamic_cast<const Sphere*>(&geometry);
        if (!sphere)
        {
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
    os << "Radius: " << radius_ << "\tInner radius: " << inner_radius_ << '\n';
}

// ------------------------------------------------------------------------- //
std::pair<double, double> Sphere::DistanceToBorder(const Vector3D& position, const Vector3D& direction)
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

    double A, B, t1, t2, difference_length_squared;

    std::pair<double, double> distance;

    double determinant;

    difference_length_squared = std::pow((position - position_).magnitude(), 2);
    A                         = difference_length_squared - radius_ * radius_;

    B = scalar_product(position - position_, direction);

    determinant = B * B - A;

    if (determinant > 0) // determinant == 0 (boundery point) is ignored
    {
        t1 = -1 * B + std::sqrt(determinant);
        t2 = -1 * B - std::sqrt(determinant);

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
        if (distance.first > 0 && distance.second > 0)
        {
            if (distance.second < distance.first)
            {
                std::swap(distance.first, distance.second);
            }
        }

    } else // particle trajectory does not have an intersection with the sphere
    {
        distance.first  = -1;
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

    if (inner_radius_ > 0)
    {
        A = difference_length_squared - inner_radius_ * inner_radius_;

        determinant = B * B - A;

        if (determinant > 0) // determinant == 0 (boundery point) is ignored
        {
            t1 = -1 * B + std::sqrt(determinant);
            t2 = -1 * B - std::sqrt(determinant);

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
            if (distance.first > 0 && distance.second > 0)
            {
                if (t1 > 0)
                {
                    if (t1 < distance.second)
                        distance.second = t1;
                }
                if (t2 > 0)
                {
                    if (t2 < distance.second)
                        distance.second = t2;
                }
            } else // The particle is inside the outer sphere
            {
                // The inner cylinder is infront of the particle trajectory
                // distance.first has to be updated
                if (t1 > 0 && t2 > 0)
                {
                    if (t1 < t2)
                        distance.first = t1;
                    else
                        distance.first = t2;
                }
                // The particle is inside the inner sphere
                // this means distance.second becomes distanc.first
                // and distance.first beomces distance to intersection with
                // the inner sphere in direction of the particle trajectory
                if ((t1 > 0 && t2 < 0) || (t2 > 0 && t1 < 0))
                {
                    std::swap(distance.first, distance.second);
                    if (t1 > 0)
                        distance.first = t1;
                    else
                        distance.first = t2;
                }
                // Now we have to check if the particle is on the border of
                // the inner sphere
                if (t1 == 0)
                {
                    // The particle is moving into the inner sphere
                    if (t2 > 0)
                    {
                        std::swap(distance.first, distance.second);
                        distance.first = t2;
                    }
                    // if not we don't have to update distance.first
                }
                if (t2 == 0)
                {
                    // The particle is moving into the inner sphere
                    if (t1 > 0)
                    {
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
