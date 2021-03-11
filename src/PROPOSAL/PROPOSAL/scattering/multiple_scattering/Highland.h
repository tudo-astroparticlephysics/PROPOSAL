
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

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/scattering/multiple_scattering/Parametrization.h"
#include <array>

/**
 * \brief This class provides the scattering routine provided by moliere.
 *
 * More precise scattering angles will be added soon.
 */
namespace PROPOSAL {
namespace multiple_scattering {
    class Highland : public Parametrization {
    protected:
        double charge;
        double radiation_length;

        bool compare(const Parametrization&) const override;
        void print(std::ostream&) const override;

    public:
        Highland(const ParticleDef&, Medium const&);

        std::unique_ptr<Parametrization> clone() const override
        {
            return std::unique_ptr<Parametrization>(
                std::make_unique<Highland>(*this));
        }

        ScatteringOffset CalculateRandomAngle(double grammage, double ei,
            double ef, const std::array<double, 4>& rnd) override;
        virtual double CalculateTheta0(double grammage, double ei, double ef);
    };
} // namespace multiple_scattering

template <typename... Args> inline auto make_highland(Args... args)
{
    return std::unique_ptr<multiple_scattering::Parametrization>(
        new multiple_scattering::Highland(std::forward<Args>(args)...));
}
} // namespace PROPOSAL
