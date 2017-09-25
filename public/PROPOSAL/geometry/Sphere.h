
#pragma once

#include "PROPOSAL/geometry/Geometry.h"

namespace PROPOSAL
{

class Sphere: public Geometry
{
    public:

        Sphere();
        Sphere(const Vector3D position, double radius, double inner_radius);
        Sphere(const Sphere&);

        Geometry* clone() const { return new Sphere(*this); };
        static Geometry* create() { return new Sphere(); }
        void swap(Geometry&);

        virtual ~Sphere() {}

        // Operators
        Sphere& operator=(const Geometry&);

        // Methods
        std::pair<double,double> DistanceToBorder(const Vector3D& position, const Vector3D& direction);

        // Getter & Setter
        double GetInnerRadius() const { return inner_radius_; }
        double GetRadius() const { return radius_; }

        void SetInnerRadius(double inner_radius) { inner_radius_ = inner_radius; };
        void SetRadius(double radius) { radius_ = radius; };

    private:

        bool compare(const Geometry&) const;
        void print(std::ostream&) const;

        double radius_;         //!< the radius of the sphere/ cylinder
        double inner_radius_;   //!< for spherical shells or hollow cylinder (0 for sphere / cylinder)
};

} /* PROPOSAL */
