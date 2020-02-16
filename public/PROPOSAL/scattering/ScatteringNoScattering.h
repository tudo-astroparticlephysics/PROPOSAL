
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

#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {

class Medium;

/**
 * \brief This class provides the scattering routine provided by moliere.
 *
 * More precise scattering angles will be added soon.
 */
class ScatteringNoScattering : public Scattering
{
public:
    ScatteringNoScattering(const ParticleDef&, std::shared_ptr<const Medium>);
    ScatteringNoScattering(const ParticleDef&, const ScatteringNoScattering&);
    ScatteringNoScattering(const ScatteringNoScattering&);
    ~ScatteringNoScattering();

    virtual Scattering* clone() const override { return new ScatteringNoScattering(*this); }
    virtual Scattering* clone(const ParticleDef& particle_def, const Utility& utility) const override
    {
        (void)utility;
        return new ScatteringNoScattering(particle_def, *this);
    }

private:
    ScatteringNoScattering& operator=(const ScatteringNoScattering&); // Undefined & not allowed

    bool compare(const Scattering&) const override;
    void print(std::ostream&) const override;

    RandomAngles CalculateRandomAngle(double dr, double ei, double ef, const Vector3D& pos, double rnd1, double rnd2, double rnd3, double rnd4) override;

    std::shared_ptr<const Medium> medium_;
};

} // namespace PROPOSAL
