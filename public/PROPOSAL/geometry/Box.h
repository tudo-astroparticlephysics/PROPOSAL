
#pragma once

#include "PROPOSAL/geometry/Geometry.h"

namespace PROPOSAL
{

class Box: public Geometry
{
    public:

        Box();
        Box(const Vector3D position, double x, double y, double z);
        Box(const Box&);

        Geometry* clone() const { return new Box(*this); };
        static Geometry* create() { return new Box(); }
        void swap(Geometry&);

        virtual ~Box() {}

        // Operators
        Box& operator=(const Geometry&);

        // Methods
        std::pair<double,double> DistanceToBorder(const Vector3D& position, const Vector3D& direction);

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

} /* PROPOSAL */
