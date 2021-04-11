
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

#include <nlohmann/json.hpp>
#include <ostream>

namespace PROPOSAL {
namespace crosssection {
    struct KinematicLimits;
}

/**
 * \brief Class containing the energy cut such as ecut and vcut.
 *
 * The cuts ecut and vcut to be used must satisfy \f$e_{cut}\geq 0\f$ and
 * \f$0 \leq v_{cut} \leq 1\f$; If both values satisfy these
 * inequalities, the lower from \f$ E \cdot v_{cut} \f$ and \f$ e_{cut}\f$
 * will be used. If only one value satisfies the inequalities,
 * only that one value is used. If both values are outside these intervals,
 * \f$v_{cut}=1\f$ is assumed.
 */

class EnergyCutSettings {
    double ecut_;
    double vcut_;
    bool continuous_randomization_;

public:
    EnergyCutSettings(double, double, bool continuous_randomization = false);
    EnergyCutSettings(const nlohmann::json&);

    bool operator==(const EnergyCutSettings& energyCutSettings) const noexcept;

    double GetCut(double energy) const;
    double GetCut(const crosssection::KinematicLimits&, double energy) const;
    size_t GetHash() const noexcept;
    double GetEcut() const noexcept { return ecut_; }
    double GetVcut() const noexcept { return vcut_; }
    bool GetContRand() const noexcept { return continuous_randomization_; }
};

std::ostream& operator<<(
    std::ostream& os, PROPOSAL::EnergyCutSettings const& cut_settings);

} // namespace PROPOSAL
