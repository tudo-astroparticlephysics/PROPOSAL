#pragma once                                
#include <functional>                                                        

#include "PROPOSAL/math/Vector3D.h"

using namespace PROPOSAL;
                                                    
class Axis;

class Density_distr
{                                                                                                   
    public:                                
        Density_distr();
        Density_distr(const Axis& axis);
        Density_distr(const Density_distr&);

        virtual ~Density_distr() {};

        virtual Density_distr* clone() const = 0;

        virtual double Correct(Vector3D xi, Vector3D direction, double res) const = 0;
        virtual double Integrate(Vector3D xi, Vector3D direction, double l) const = 0;
        virtual double Calculate(Vector3D xi, Vector3D direction, double distance) const = 0;
        virtual double GetCorrection(Vector3D x) const = 0;

    protected:
        Axis* axis_;
};

class Axis
{
    public:

        Axis();
        Axis(Vector3D fp0, Vector3D fAxis);

        virtual Axis* clone() const = 0;

        virtual double GetDepth(Vector3D xi) const = 0;
        virtual double GetEffectiveDistance(Vector3D xi, Vector3D direction) const = 0;
        
        Vector3D GetAxis() { return fAxis_; };
        Vector3D GetFp0() { return fp0_; };

    protected:
        Vector3D fAxis_;
        Vector3D fp0_;
};

class RadialAxis : public Axis
{
    public:
        RadialAxis();
        RadialAxis(Vector3D fAxis, Vector3D fp0);
        
        RadialAxis* clone() const {return new RadialAxis(*this);};

        double GetDepth(Vector3D xi) const override;
        double GetEffectiveDistance(Vector3D xi, Vector3D direction) const override;
};

class CartesianAxis : public Axis
{
    public:
        CartesianAxis();
        CartesianAxis(Vector3D fAxis, Vector3D fp0);
        
        CartesianAxis* clone() const {return new CartesianAxis(*this);};

        double GetDepth(Vector3D xi) const override;
        double GetEffectiveDistance(Vector3D xi, Vector3D direction) const override;
};
