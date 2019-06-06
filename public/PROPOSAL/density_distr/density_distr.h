#pragma once                                
#include <functional>                                                        

#include "PROPOSAL/math/Vector3D.h"

using namespace PROPOSAL;
                                                    
class Density_distr
{                                                                                                   
    public:                                
        Density_distr();
        Density_distr(Vector3D fAxis, Vector3D fp0, std::function<double(double)> density_distribution);

        ~Density_distr();

        virtual Density_distr* clone() const = 0;
                                
        virtual double Integrate(Vector3D xi, Vector3D res)=0;

        double GetDepth(Vector3D xi);
                                        
        std::function<double(double)> GetDensityDistribution();
        Vector3D GetAxis();
        Vector3D GetFp0();
                                
    protected:                                
        Vector3D xi_;
        Vector3D const fAxis_;
        Vector3D const fp0_;
        
        std::function<double(double)> density_distribution_;
};
