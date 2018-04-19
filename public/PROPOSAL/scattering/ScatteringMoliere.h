
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

#include <vector>

#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {

class Medium;

class ScatteringMoliere : public Scattering
{
public:
    // constructor
    ScatteringMoliere(Particle&, const Medium&);
    ScatteringMoliere(Particle&, const ScatteringMoliere&);
    ScatteringMoliere(const ScatteringMoliere&);
    ~ScatteringMoliere();

    Scattering* clone() const { return new ScatteringMoliere(*this); }
    virtual Scattering* clone(Particle& particle, const Utility& utility) const
    {
        (void)utility;
        return new ScatteringMoliere(particle, *this);
    }

private:
    ScatteringMoliere& operator=(const ScatteringMoliere&); // Undefined & not allowed

    bool compare(const Scattering&) const;

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef);

    const Medium* medium_;

    int numComp_; // number of components in medium
    double ZSq_A_average_;
    std::vector<double> Zi_;        // nuclear charge of different components
    std::vector<double> weight_ZZ_; // mass weights of different components time Z^2
    double weight_ZZ_sum_;          // inverse of sum of mass weights of different components time Z^2
    int max_weight_index_;          // index of the maximium of mass weights of different components

    // scattering parameters
    double chiCSq_; // characteristic angle² in rad²
    std::vector<double> B_;

    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//

    double f1M(double x);
    double f2M(double x);

    double f(double theta);

    double F1M(double x);
    double F2M(double x);

    double F(double theta);

    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//

    double GetRandom(double pre_factor);
};
} // namespace PROPOSAL
