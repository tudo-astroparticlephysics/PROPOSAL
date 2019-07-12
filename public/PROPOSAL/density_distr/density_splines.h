#pragma once                                
#include <functional>
#include <vector>
#include <iostream>
#include "PROPOSAL/density_distr/density_polynomial.h" 
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/math/MathMethods.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class Density_splines : public Density_distr
{                                                                                                   
    public:                                
        Density_splines(const Axis&, 
                        const std::vector<SplineCoefficients>&);
        Density_splines(const Density_splines&);

        Density_splines* clone() const { return new Density_splines(*this); };

        double Correct(Vector3D xi, Vector3D direction, double res) const override;
        double Integrate(Vector3D xi, Vector3D direction, double l) const override;
        double Calculate(Vector3D xi, Vector3D direction, double distance) const override;
        double GetCorrection(Vector3D xi) const override;
        
        double Helper_function(Vector3D xi, Vector3D direction, double res, double l) const ;
        double helper_function(Vector3D xi, Vector3D direction, double res, double l) const ;

    protected:
        std::vector<Density_polynomial*> dens_polynom_;
        std::vector<double> definition_area_;
        
        std::function<double(double)> density_distribution;
        std::function<double(double)> antiderived_density_distribution;
};

