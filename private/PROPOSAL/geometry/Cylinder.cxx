
#include <cmath>
#include <vector>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/geometry/Cylinder.h"

using namespace PROPOSAL;

Cylinder::Cylinder()
    : Geometry("Cylinder")
    , radius_(0.0)
    , inner_radius_(0.0)
    , z_(0.0)
{
    // Do nothing here
}

Cylinder::Cylinder(const Vector3D position, double radius, double inner_radius, double z)
    : Geometry("Cylinder", position)
    , radius_(100 * radius)
    , inner_radius_(100 * inner_radius)
    , z_(100 * z)
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
    Cylinder* cylinder = dynamic_cast<Cylinder*>(&geometry);
    if (!cylinder)
    {
        log_warn("Cannot swap Cylinder!");
        return;
    }

    Geometry::swap(*cylinder);

    std::swap(inner_radius_, cylinder->inner_radius_);
    std::swap(radius_, cylinder->radius_);
    std::swap(z_, cylinder->z_);
}

//------------------------------------------------------------------------- //
Cylinder& Cylinder::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const Cylinder* cylinder = dynamic_cast<const Cylinder*>(&geometry);
        if (!cylinder)
        {
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
    os << "Radius: " << radius_ << "\tInnner radius: " << inner_radius_ << " Height: " << z_ << '\n';
}

// ------------------------------------------------------------------------- //
std::pair<double, double> Cylinder::DistanceToBorder(const Vector3D& position, const Vector3D& direction)
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
    double dir_vec_x = direction.GetX();
    double dir_vec_y = direction.GetY();
    double dir_vec_z = direction.GetZ();

    double determinant;

    double intersection_x;
    double intersection_y;
    double intersection_z;

    std::vector<double> dist;

    std::pair<double, double> distance;

    double z_calc_pos = position_.GetZ() + 0.5 * z_;
    double z_calc_neg = position_.GetZ() - 0.5 * z_;

    if (!(dir_vec_x == 0 && dir_vec_y == 0)) // Otherwise the particle
                                             // trajectory is parallel to
                                             // cylinder barrel
    {

        A = std::pow((position.GetX() - position_.GetX()), 2) +
            std::pow((position.GetY() - position_.GetY()), 2) -
            radius_*radius_;

        B = 2 * ((position.GetX() - position_.GetX()) * dir_vec_x + (position.GetY() - position_.GetY()) * dir_vec_y);

        C = dir_vec_x * dir_vec_x + dir_vec_y * dir_vec_y;

        B /= C;
        A /= C;

        determinant = 0.25 * B*B - A;

        if (determinant > 0) // determinant == 0 (boundery point) is ignored
        {
            t1 = -1 * B / 2 + std::sqrt(determinant);
            t2 = -1 * B / 2 - std::sqrt(determinant);

            // Computer precision controll
            if (t1 > 0 && t1 < GEOMETRY_PRECISION)
                t1 = 0;
            if (t2 > 0 && t2 < GEOMETRY_PRECISION)
                t2 = 0;

            if (t1 > 0)
            {
                intersection_z = position.GetZ() + t1 * dir_vec_z;
                // is inside the borders
                if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                {
                    dist.push_back(t1);
                }
            }

            if (t2 > 0)
            {
                intersection_z = position.GetZ() + t2 * dir_vec_z;
                // is inside the borders
                if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                {
                    dist.push_back(t2);
                }
            }
        }
    }

    // if we have found already to intersections we don't have to check for
    // intersections
    // with top or bottom surface
    if (dist.size() < 2)
    {
        // intersection with E1
        if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel
                            // to E1 (Should not happen)
        {
            t = (z_calc_pos - position.GetZ()) / dir_vec_z;
            // Computer precision controll
            if (t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            if (t > 0) // Interection is in particle trajectory direction
            {
                intersection_x = position.GetX() + t * dir_vec_x;
                intersection_y = position.GetY() + t * dir_vec_y;

                if (std::sqrt(std::pow((intersection_x - position_.GetX()), 2) +
                    std::pow((intersection_y - position_.GetY()), 2)) <=
                        radius_ &&
                    std::sqrt(std::pow((intersection_x - position_.GetX()), 2) +
                    std::pow((intersection_y - position_.GetY()), 2)) >=
                        inner_radius_)
                {
                    dist.push_back(t);
                }
            }
        }

        // intersection with E2
        if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel
                            // to E2 (Should not happen)
        {
            t = (z_calc_neg - position.GetZ()) / dir_vec_z;

            // Computer precision controll
            if (t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            if (t > 0) // Interection is in particle trajectory direction
            {
                intersection_x = position.GetX() + t * dir_vec_x;
                intersection_y = position.GetY() + t * dir_vec_y;

                if (std::sqrt(std::pow((intersection_x - position_.GetX()), 2) +
                    std::pow((intersection_y - position_.GetY()), 2)) <=
                        radius_ &&
                    std::sqrt(std::pow((intersection_x - position_.GetX()), 2) +
                    std::pow((intersection_y - position_.GetY()), 2)) >=
                        inner_radius_)
                {
                    dist.push_back(t);
                }
            }
        }
    }
    // No intersection with the outer cylinder
    if (dist.size() < 1)
    {
        distance.first  = -1;
        distance.second = -1;
        //    return distance;
    } else if (dist.size() == 1) // particle is inside the cylinder
    {
        distance.first  = dist.at(0);
        distance.second = -1;

    } else if (dist.size() == 2) // cylinder is infront of the particle
    {
        distance.first  = dist.at(0);
        distance.second = dist.at(1);

        if (distance.second < distance.first)
        {
            std::swap(distance.first, distance.second);
        }

    } else
    {
        log_error("This point should never be reached");
    }
    // This cylinder might be hollow and we have to check if the inner border is
    // reached before.
    // So we caluculate the intersection with the inner cylinder.

    if (inner_radius_ > 0)
    {
        if (!(dir_vec_x == 0 && dir_vec_y == 0))
        {

            A = std::pow((position.GetX() - position_.GetX()), 2) +
                std::pow((position.GetY() - position_.GetY()), 2) -
                inner_radius_*inner_radius_;

            B = 2 *
                ((position.GetX() - position_.GetX()) * dir_vec_x + (position.GetY() - position_.GetY()) * dir_vec_y);

            C = dir_vec_x * dir_vec_x + dir_vec_y * dir_vec_y;

            B /= C;
            A /= C;

            determinant = 0.25 * B*B - A;

            if (determinant > 0) // determinant == 0 (boundery point) is ignored
            {
                t1 = -1 * B / 2 + std::sqrt(determinant);
                t2 = -1 * B / 2 - std::sqrt(determinant);

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
                if (distance.first > 0 && distance.second > 0)
                {
                    if (t1 > 0)
                    {
                        intersection_z = position.GetZ() + t1 * dir_vec_z;
                        // is inside the borders
                        if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                        {
                            if (t1 < distance.second)
                                distance.second = t1;
                        }
                    }

                    if (t2 > 0)
                    {
                        intersection_z = position.GetZ() + t2 * dir_vec_z;
                        // is inside the borders
                        if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                        {
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
                else if (distance.first < 0 && distance.second < 0)
                {
                    if (t1 > 0)
                    {
                        intersection_z = position.GetZ() + t1 * dir_vec_z;
                        // is inside the borders
                        if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                        {
                            distance.first = t1;
                        }
                    }
                    if (t2 > 0)
                    {
                        intersection_z = position.GetZ() + t2 * dir_vec_z;
                        // is inside the borders
                        if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                        {
                            distance.first = t2;
                        }
                    }
                }
                // The particle is inside the outer cylinder
                // or particle hits bottom or top surface and then hits
                // the inner cylinder barrel
                else
                {
                    // The inner cylinder is infront of the particle trajectory
                    // distance.first has to be updated
                    if (t1 > 0 && t2 > 0)
                    {
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
                        if (position.GetZ() >= z_calc_neg && position.GetZ() <= z_calc_pos &&
                            std::sqrt(std::pow((position.GetX() - position_.GetX()), 2) +
                            std::pow((position.GetY() - position_.GetY()), 2)) <=
                                radius_ + GEOMETRY_PRECISION &&
                            std::sqrt(std::pow((position.GetX() - position_.GetX()), 2) +
                            std::pow((position.GetY() - position_.GetY()), 2)) >=
                                inner_radius_ - GEOMETRY_PRECISION)
                        {
                            if (t1 < distance.first)
                            {
                                intersection_z = position.GetZ() + t1 * dir_vec_z;
                                // is inside the borders
                                if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                                {
                                    // This case means particle is inside in
                                    // hits
                                    // inner cylinder first
                                    distance.first = t1;
                                }
                            }
                            if (t2 < distance.first)
                            {
                                intersection_z = position.GetZ() + t2 * dir_vec_z;
                                // is inside the borders
                                if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                                {
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
                        else
                        {
                            intersection_z = position.GetZ() + t1 * dir_vec_z;
                            // is inside the borders
                            if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                            {
                                // This case means particle is inside in hits
                                // inner cylinder first
                                distance.second = t1;
                            }

                            if (distance.second < 0)
                            {
                                intersection_z = position.GetZ() + t2 * dir_vec_z;
                                // is inside the borders
                                if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                                {
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
                    if ((t1 > 0 && t2 < 0) || (t2 > 0 && t1 < 0))
                    {
                        if (t1 > 0)
                        {
                            intersection_z = position.GetZ() + t1 * dir_vec_z;
                            // is inside the borders
                            if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                            {
                                std::swap(distance.first, distance.second);
                                distance.first = t1;
                            }
                        } else
                        {
                            intersection_z = position.GetZ() + t2 * dir_vec_z;
                            // is inside the borders
                            if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                            {
                                std::swap(distance.first, distance.second);
                                distance.first = t2;
                            }
                        }
                    }
                    // Now we have to check if the particle is on the border of
                    // the inner sphere
                    if (t1 == 0)
                    {
                        // The particle is moving into the inner cylinder
                        if (t2 > 0)
                        {
                            intersection_z = position.GetZ() + t2 * dir_vec_z;
                            // is inside the borders
                            if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                            {
                                std::swap(distance.first, distance.second);
                                distance.first = t2;
                            }
                        }
                        // if not we don't have to update distance.first
                    }
                    if (t2 == 0)
                    {
                        // The particle is moving into the inner sphere
                        if (t1 > 0)
                        {
                            intersection_z = position.GetZ() + t1 * dir_vec_z;
                            // is inside the borders
                            if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                            {
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
