
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

#include <array>
#include <vector>

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/scattering/multiple_scattering/Parametrization.h"

namespace PROPOSAL {
struct ParticleDef;

namespace multiple_scattering {

    class Moliere : public Parametrization {
        bool compare(const Parametrization&) const override;
        void print(std::ostream&) const override;

        std::unique_ptr<Parametrization> clone() const final
        {
            return std::unique_ptr<Parametrization>(
                std::make_unique<Moliere>(*this));
        }

        int numComp_; // number of components in medium
        double ZSq_A_average_;
        std::vector<double> Zi_; // nuclear charge of different components
        std::vector<double>
            weight_ZZ_;        // mass weights of different components time Z^2
        double weight_ZZ_sum_; // inverse of sum of mass weights of different
                               // components time Z^2
        int max_weight_index_; // index of the maximium of mass weights of
                               // different components

        // scattering parameters
        double chiCSq_; // characteristic angle² in rad²
        std::vector<double> B_;

        double f1M(double x);
        double f2M(double x);

        double f(double theta);

        double F1M(double x);
        double F2M(double x);

        double F(double theta);

        double GetRandom(double pre_factor, double rnd);

    public:
        // constructor
        Moliere(const ParticleDef&, Medium const&);

        ScatteringOffset CalculateRandomAngle(double grammage, double ei,
            double ef, const std::array<double, 4>& rnd) override;
    };
} // namespace multiple_scattering

template <typename... Args> inline auto make_moliere(Args... args)
{
    return std::unique_ptr<multiple_scattering::Parametrization>(
        new multiple_scattering::Moliere(std::forward<Args>(args)...));
}

} // namespace PROPOSAL
