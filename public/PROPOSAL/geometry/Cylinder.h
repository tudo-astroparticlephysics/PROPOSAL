
#pragma once

#include "PROPOSAL/geometry/Geometry.h"

namespace PROPOSAL
{

class Cylinder: public Geometry
{
    public:

        Cylinder();
        Cylinder(Vector3D position, double radius, double inner_radius, double z);
        Cylinder(const Cylinder&);

        Geometry* clone() const { return new Cylinder(*this); };
        static Geometry* create() { return new Cylinder(); }
        void swap(Geometry&);

        virtual ~Cylinder() {}

        // Operators
        Cylinder& operator=(const Geometry&);

        // Methods
        std::pair<double,double> DistanceToBorder(Vector3D& position, Vector3D& direction);

        // Getter & Setter
        double GetInnerRadius() const { return inner_radius_; }
        double GetRadius() const { return radius_; }
        double GetZ() const { return z_; }

        void SetInnerRadius(double inner_radius) { inner_radius_ = inner_radius; };
        void SetRadius(double radius)            { radius_ = radius; };
        void SetZ(double z)                      { z_ = z; };

    private:

        bool compare(const Geometry&) const;
        void print(std::ostream&) const;

        double radius_;         //!< the radius of the sphere/ cylinder
        double inner_radius_;   //!< for spherical shells or hollow cylinder (0 for sphere / cylinder)
        double z_;              //!< height of box/cylinder
};

} /* PROPOSAL */
