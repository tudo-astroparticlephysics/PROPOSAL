/**
 * This class contains geometry information
 * for each ProcessCollection
 *
 * @author Jan-Hendrik Köhne
 */
#pragma once

// #include <string>
// #include <utility>
#include <map>
#include <iostream>

#include "PROPOSAL/math/Vector3D.h"

namespace PROPOSAL{

class Geometry
{
    public:

        Geometry(const std::string); //TODO(mario): Maybe not useful Tue 2017/08/01
        Geometry(const std::string, const Vector3D position);
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

        bool IsInside(const Vector3D& position, const Vector3D& direction);


        bool IsInfront(const Vector3D& position, const Vector3D& direction);


        bool IsBehind(const Vector3D& position, const Vector3D& direction);

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
        virtual std::pair<double,double> DistanceToBorder(const Vector3D& position, const Vector3D& direction) = 0;

        /*!
         * Calculates the distance to the closest approch to the geometry center
         */
        double DistanceToClosestApproach(const Vector3D& position, const Vector3D& direction);


        // void swap(Geometry &geometry);

        // ----------------------------------------------------------------- //
        // Getter & Setter
        // ----------------------------------------------------------------- //

        Vector3D GetPosition() const { return position_; }

        std::string GetName() const { return name_; }
        unsigned int GetHirarchy() const { return hirarchy_; }

        void SetPosition(const Vector3D& position) { position_ = position; };

        void SetHirarchy(unsigned int hirarchy) { hirarchy_ = hirarchy; };

    protected:

        // Implemented in child classes to be able to use equality operator
        virtual bool compare(const Geometry&) const = 0;
        virtual void print(std::ostream&) const = 0;

        Vector3D position_;    //!< x,y,z-coordinate of origin ( center of box, cylinder, sphere)

        std::string name_;      //!< "box" , "cylinder" , "sphere" (sphere and cylinder might be hollow)

        unsigned int hirarchy_; //!< adds a hirarchy of geometry objects to allow crossing geometries

};

}
