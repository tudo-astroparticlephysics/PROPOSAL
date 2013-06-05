/**
 * This class contains geometry information
 * for each ProcessCollection
 *
 * @author Jan-Hendrik Köhne
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "PROPOSAL/Particle.h"

class Geometry
{
private:

    double x0_;             //!< x-coordinate of origin ( center of box, cylinder, sphere)
    double y0_;             //!< y-coordinate of origin ( center of box, cylinder, sphere)
    double z0_;             //!< z-coordinate of origin ( center of box, cylinder, sphere)

    double inner_radius_;   //!< for spherical shells or hollow cylinder (0 for sphere / cylinder)
    double radius_;         //!< the radius of the sphere/ cylinder

    double x_;              //!< width of box in x-direction
    double y_;              //!< width of box in y-direction
    double z_;              //!< height of box/cylinder



public:

    //Constructors

    Geometry();
    Geometry(const Geometry&);
    Geometry& operator=(const Geometry&);
    bool operator==(const Geometry &geometry) const;
    bool operator!=(const Geometry &geometry) const;

//----------------------------------------------------------------------------//
    //Memberfunctions

    bool IsParticleInside(Particle* particle);

//----------------------------------------------------------------------------//

    double DistanceToBorder(Particle* particle);

//----------------------------------------------------------------------------//
    void swap(Geometry &geometry);

//----------------------------------------------------------------------------//
    //Getter
    double GetX0() const {
        return x0_;
    }

    double GetY0() const {
        return y0_;
    }

    double GetZ0() const {
        return z0_;
    }

    double GetX() const {
        return x_;
    }

    double GetY() const {
        return y_;
    }

    double GetZ() const {
        return z_;
    }

    double GetInnerRadius() const {
        return inner_radius_;
    }

    double GetRadius() const {
        return radius_;
    }

//----------------------------------------------------------------------------//
    //Setter
    void SetX0(double x0);
    void SetY0(double y0);
    void SetZ0(double z0);
    void SetX(double x);
    void SetY(double y);
    void SetZ(double z);
    void SetInnerRadius(double inner_radius_);
    void SetRadius(double radius);
//----------------------------------------------------------------------------//
    //Destructor
    ~Geometry();
};

#endif // GEOMETRY_H
