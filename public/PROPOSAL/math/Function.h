
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


#pragma once

#include <functional>
#include <vector>
#include <fstream>

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Polynom      %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace PROPOSAL {

class Polynom {
   public:
    Polynom(std::vector<double> coefficients);
    Polynom(const Polynom&);

    ~Polynom() {};
    Polynom& operator=(const Polynom& poly) = default;

    Polynom* clone() const { return new Polynom(*this); };

    bool operator==(const Polynom& polynom) const;
    bool operator!=(const Polynom& polynom) const;

    double evaluate(double x);
    void shift(double x);

    Polynom GetDerivative();
    Polynom GetAntiderivative(double constant);
    std::vector<double> GetCoefficient() const;

    friend std::ostream& operator<<(std::ostream& os, const Polynom& p);

    std::function<double(double)> GetFunction();

   protected:
    int N_;
    double* coeff_;
};

}  // namespace PROPOSAL
