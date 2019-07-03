#pragma once                                
#include <functional>
#include <vector>
#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/math/Vector3D.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Polynom      %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class Polynom
{
    public:
        Polynom(std::vector<double> coefficients);
        Polynom(const Polynom&);

        Polynom* clone() const { return new Polynom(*this); };

        double evaluate(double x);

        Polynom GetDerivative();
        Polynom GetAntiderivative(double constant);

        std::function<double(double)> GetFunction();

    protected:
        std::function<double(double)> function_;
        int N;
        double* coeff;
};

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class Density_polynomial : Density_distr
{                                                                                                   
    public:                                
        Density_polynomial(const Axis&, const Polynom&);
        Density_polynomial(const Density_polynomial&);

        Density_polynomial* clone() const { return new Density_polynomial(*this); };

        double Correct(Vector3D xi, Vector3D direction, double res) const override;
        double Integrate(Vector3D xi, Vector3D direction, double l) const override;
        double Calculate(Vector3D xi, Vector3D direction, double distance) const override;
        

    protected:
        Polynom polynom_;

        std::function<double(double)> density_distribution;
        std::function<double(double)> antiderived_density_distribution;
};

