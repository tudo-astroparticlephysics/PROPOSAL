#pragma once                                
#include <functional>                                                        

#include "PROPOSAL/math/Vector3D.h"

using namespace PROPOSAL;
                                                    
class Density_distr
{                                                                                                   
    public:                                
        Density_distr(double rho, std::function<double(double)> density_distribution);                                                             
        ~Density_distr();

        virtual Density_distr* clone() const = 0;
                                
        virtual double Integrate(double x_i, double res)=0;                                
        virtual double Integrate(double res)=0;                                
    
        void SetXi(double);                                
                                        
        std::function<double(double)> GetDensityDistribution();
        double GetXi();
        double GetRho();
                                
    private:                                
        double xi_;
        double rho_;
        
        std::function<double(double)> density_distribution_;
};
