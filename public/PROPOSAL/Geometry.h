/**
 * This class contains geometry information
 * for each ProcessCollection
 *
 * @author Jan-Hendrik Köhne
 */
#pragma once

#ifndef GEOMETRY_H
#define GEOMETRY_H

// #include <string>
// #include <utility>
// #include <iostream>

#include "PROPOSAL/PROPOSALParticle.h"


namespace PROPOSAL{

class Geometry
{
    public:

        Geometry(std::string); //TODO(mario): Maybe not useful Tue 2017/08/01
        Geometry(std::string, double x, double y, double z);
        Geometry(const Geometry&);

        virtual Geometry* clone() const = 0; // virtual constructor idiom (used for deep copies)
        virtual void swap(Geometry&);

        virtual ~Geometry() {};

        // Operators
        virtual Geometry& operator=(const Geometry&);
        bool operator==(const Geometry &geometry) const;
        bool operator!=(const Geometry &geometry) const;
        friend std::ostream& operator<<(std::ostream&, Geometry const&);


        // ----------------------------------------------------------------- //
        // Member functions
        // ----------------------------------------------------------------- //

        bool IsInside(PROPOSALParticle* particle);


        bool IsInfront(PROPOSALParticle* particle);


        bool IsBehind(PROPOSALParticle* particle);

        /*!
         * This function calculates the distance of the particle position
         * to the border of the geometry in direction of the particle trajectory.
         * If the particle trajectory does not have an intersection with the geometry
         * (-1 /-1) is returned
         * If the particle trajectory has two intersections (dist_1 /dist_2) is returned
         * If the particle has one intersection (dist_1 /-1) is returned
         * (one intersection means one intersection in direction of the particle trajectory
         * and one in the opposite direction. Cause we are not intersted in this one. it is set to -1)
         * Note: If the particle is on the geometry border this is not treated as an intersection
         * A particle on the geometry border which moves inside has one intersection,
         * a particle on the geometry border which moves outside has no intersection.
         * Distances smaller then GEOMETRY_PRECISION (1e-9) are also set to -1
         */
        virtual std::pair<double,double> DistanceToBorder(PROPOSALParticle* particle) = 0;

        /*!
         * Calculates the distance to the closest approch to the geometry center
         */
        double DistanceToClosestApproach(PROPOSALParticle* particle);


        // void swap(Geometry &geometry);

        // ----------------------------------------------------------------- //
        // Getter & Setter
        // ----------------------------------------------------------------- //

        double GetX0() const { return x0_; }
        double GetY0() const { return y0_; }
        double GetZ0() const { return z0_; }


        std::string GetName() const    { return name_; }
        unsigned int GetHirarchy() const { return hirarchy_; }

        void SetX0(double x0) { x0_ = x0; };
        void SetY0(double y0) { y0_ = y0; };
        void SetZ0(double z0) { z0_ = z0; };

        void SetHirarchy(unsigned int hirarchy) { hirarchy_ = hirarchy; };

    protected:

        // Implemented in child classes to be able to use equality operator
        virtual bool compare(const Geometry&) const = 0;
        virtual void print(std::ostream&) const = 0;

        double x0_;             //!< x-coordinate of origin ( center of box, cylinder, sphere)
        double y0_;             //!< y-coordinate of origin ( center of box, cylinder, sphere)
        double z0_;             //!< z-coordinate of origin ( center of box, cylinder, sphere)

        std::string name_;    //!< "box" , "cylinder" , "sphere" (sphere and cylinder might be hollow)

        unsigned int hirarchy_; //!< adds a hirarchy of geometry objects to allow crossing geometries

};


// ----------------------------------------------------------------------------
/// @brief Class Sphere
// ----------------------------------------------------------------------------
class Sphere: public Geometry
{
    public:

        Sphere();
        Sphere(double x0, double y0, double z0, double radius, double inner_radius);
        Sphere(const Sphere&);

        Geometry* clone() const { return new Sphere(*this); };
        void swap(Geometry&);

        virtual ~Sphere() {}

        // Operators
        Sphere& operator=(const Geometry&);

        // Methods
        std::pair<double,double> DistanceToBorder(PROPOSALParticle* particle);

        // Getter & Setter
        double GetInnerRadius() const { return inner_radius_; }
        double GetRadius() const      { return radius_; }

        void SetInnerRadius(double inner_radius) { inner_radius_ = inner_radius; };
        void SetRadius(double radius)            { radius_ = radius; };

    private:

        bool compare(const Geometry&) const;
        void print(std::ostream&) const;

        double radius_;         //!< the radius of the sphere/ cylinder
        double inner_radius_;   //!< for spherical shells or hollow cylinder (0 for sphere / cylinder)
};

// ----------------------------------------------------------------------------
/// @brief Class Box
// ----------------------------------------------------------------------------
class Box: public Geometry
{
    public:

        Box();
        Box(double x0, double y0, double z0, double x, double y, double z);
        Box(const Box&);

        Geometry* clone() const { return new Box(*this); };
        void swap(Geometry&);

        virtual ~Box() {}

        // Operators
        Box& operator=(const Geometry&);

        // Methods
        std::pair<double,double> DistanceToBorder(PROPOSALParticle* particle);

        // Getter & Setter
        double GetX() const { return x_; }
        double GetY() const { return y_; }
        double GetZ() const { return z_; }

        void SetX(double x) { x_ = x; };
        void SetY(double y) { y_ = y; };
        void SetZ(double z) { z_ = z; };


    private:

        bool compare(const Geometry&) const;
        void print(std::ostream&) const;


        double x_;              //!< width of box in x-direction
        double y_;              //!< width of box in y-direction
        double z_;              //!< width of box in z-direction

};

// ----------------------------------------------------------------------------
/// @brief Class Cylinder
// ----------------------------------------------------------------------------
class Cylinder: public Geometry
{
    public:

        Cylinder();
        Cylinder(double x0, double y0, double z0, double radius, double inner_radius, double z);
        Cylinder(const Cylinder&);

        Geometry* clone() const { return new Cylinder(*this); };
        void swap(Geometry&);

        virtual ~Cylinder() {}

        // Operators
        Cylinder& operator=(const Geometry&);

        // Methods
        std::pair<double,double> DistanceToBorder(PROPOSALParticle* particle);

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

}

#endif // GEOMETRY_H
