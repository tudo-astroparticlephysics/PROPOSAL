
#include <cmath>
#include <vector>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/geometry/Cylinder.h"

using namespace PROPOSAL;


Cylinder::Cylinder(const Vector3D& position, double z, double radius, double inner_radius)
    : Geometry("Cylinder", position)
    , radius_(radius)
    , inner_radius_(inner_radius)
    , z_(z)
{
    Logging::Get("proposal.geometry")->info("Order of function parameters for Cylinder constructor has been changed in vesion 7.");
    if (inner_radius_ > radius_)
    {
        Logging::Get("proposal.geometry")->error("Inner radius {} is greater then radius {} (will be swaped)", inner_radius_, radius_);
        std::swap(inner_radius_, radius_);
    }
    if (inner_radius_ == radius_)
    {
        Logging::Get("proposal.geometry")->error("Warning: Inner radius {} == radius {} (Volume is 0)", inner_radius_, radius_);
    }
}


Cylinder::Cylinder(const nlohmann::json& config)
    : Geometry(config)
{
    assert(config["outer_radius"].is_number());
    assert(config["height"].is_number());


    radius_ = config["outer_radius"].get<double>();
    inner_radius_ = config.value("inner_radius", 0);
    z_ = config["height"].get<double>();

    assert(inner_radius_>=0);
    assert(radius_>inner_radius_);
    assert(z_>0);
}

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

std::pair<double, double> Cylinder::DistanceToBorder(const Vector3D& position, const Vector3D& direction) const
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
    auto dir_vec = Cartesian3D(direction);
    auto pos_vec = Cartesian3D(position);

    double determinant;

    double intersection_x;
    double intersection_y;
    double intersection_z;

    std::vector<double> dist;

    std::pair<double, double> distance;

    double z_calc_pos = position_.GetZ() + 0.5 * z_;
    double z_calc_neg = position_.GetZ() - 0.5 * z_;

    if (!(dir_vec.GetX() == 0 && dir_vec.GetY() == 0)) // Otherwise the particle
                                             // trajectory is parallel to
                                             // cylinder barrel
    {

        A = std::pow((pos_vec.GetX() - position_.GetX()), 2) +
            std::pow((pos_vec.GetY() - position_.GetY()), 2) -
            radius_*radius_;

        B = 2 * ((pos_vec.GetX() - position_.GetX()) * dir_vec.GetX() + (pos_vec.GetY() - position_.GetY()) * dir_vec.GetY());

        C = dir_vec.GetX() * dir_vec.GetX() + dir_vec.GetY() * dir_vec.GetY();

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
                intersection_z = pos_vec.GetZ() + t1 * dir_vec.GetZ();
                // is inside the borders
                if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                {
                    dist.push_back(t1);
                }
            }

            if (t2 > 0)
            {
                intersection_z = pos_vec.GetZ() + t2 * dir_vec.GetZ();
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
        if (dir_vec.GetZ() != 0) // if dir_vec == 0 particle trajectory is parallel
                            // to E1 (Should not happen)
        {
            t = (z_calc_pos - pos_vec.GetZ()) / dir_vec.GetZ();
            // Computer precision controll
            if (t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            if (t > 0) // Interection is in particle trajectory direction
            {
                intersection_x = pos_vec.GetX() + t * dir_vec.GetX();
                intersection_y = pos_vec.GetY() + t * dir_vec.GetY();

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
        if (dir_vec.GetZ() != 0) // if dir_vec == 0 particle trajectory is parallel
                            // to E2 (Should not happen)
        {
            t = (z_calc_neg - pos_vec.GetZ()) / dir_vec.GetZ();

            // Computer precision controll
            if (t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            if (t > 0) // Interection is in particle trajectory direction
            {
                intersection_x = pos_vec.GetX() + t * dir_vec.GetX();
                intersection_y = pos_vec.GetY() + t * dir_vec.GetY();

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
        Logging::Get("proposal.geometry")->error("This point should never be reached");
    }
    // This cylinder might be hollow and we have to check if the inner border is
    // reached before.
    // So we caluculate the intersection with the inner cylinder.

    if (inner_radius_ > 0)
    {
        if (!(dir_vec.GetX() == 0 && dir_vec.GetY() == 0))
        {

            A = std::pow((pos_vec.GetX() - position_.GetX()), 2) +
                std::pow((pos_vec.GetY() - position_.GetY()), 2) -
                inner_radius_*inner_radius_;

            B = 2 *
                ((pos_vec.GetX() - position_.GetX()) * dir_vec.GetX() + (pos_vec.GetY() - position_.GetY()) * dir_vec.GetY());

            C = dir_vec.GetX() * dir_vec.GetX() + dir_vec.GetY() * dir_vec.GetY();

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
                        intersection_z = pos_vec.GetZ() + t1 * dir_vec.GetZ();
                        // is inside the borders
                        if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                        {
                            if (t1 < distance.second)
                                distance.second = t1;
                        }
                    }

                    if (t2 > 0)
                    {
                        intersection_z = pos_vec.GetZ() + t2 * dir_vec.GetZ();
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
                        intersection_z = pos_vec.GetZ() + t1 * dir_vec.GetZ();
                        // is inside the borders
                        if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                        {
                            distance.first = t1;
                        }
                    }
                    if (t2 > 0)
                    {
                        intersection_z = pos_vec.GetZ() + t2 * dir_vec.GetZ();
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
                        if (pos_vec.GetZ() >= z_calc_neg && pos_vec.GetZ() <= z_calc_pos &&
                            std::sqrt(std::pow((pos_vec.GetX() - position_.GetX()), 2) +
                            std::pow((pos_vec.GetY() - position_.GetY()), 2)) <=
                                radius_ + GEOMETRY_PRECISION &&
                            std::sqrt(std::pow((pos_vec.GetX() - position_.GetX()), 2) +
                            std::pow((pos_vec.GetY() - position_.GetY()), 2)) >=
                                inner_radius_ - GEOMETRY_PRECISION)
                        {
                            if (t1 < distance.first)
                            {
                                intersection_z = pos_vec.GetZ() + t1 * dir_vec.GetZ();
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
                                intersection_z = pos_vec.GetZ() + t2 * dir_vec.GetZ();
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
                            intersection_z = pos_vec.GetZ() + t1 * dir_vec.GetZ();
                            // is inside the borders
                            if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                            {
                                // This case means particle is inside in hits
                                // inner cylinder first
                                distance.second = t1;
                            }

                            if (distance.second < 0)
                            {
                                intersection_z = pos_vec.GetZ() + t2 * dir_vec.GetZ();
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
                            intersection_z = pos_vec.GetZ() + t1 * dir_vec.GetZ();
                            // is inside the borders
                            if (intersection_z > z_calc_neg && intersection_z < z_calc_pos)
                            {
                                std::swap(distance.first, distance.second);
                                distance.first = t1;
                            }
                        } else
                        {
                            intersection_z = pos_vec.GetZ() + t2 * dir_vec.GetZ();
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
                            intersection_z = pos_vec.GetZ() + t2 * dir_vec.GetZ();
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
                            intersection_z = pos_vec.GetZ() + t1 * dir_vec.GetZ();
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
