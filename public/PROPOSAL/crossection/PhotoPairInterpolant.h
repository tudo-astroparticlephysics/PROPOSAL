
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

#include "PROPOSAL/crossection/CrossSectionInterpolant.h"
using std::shared_ptr;

namespace PROPOSAL {
class PhotoPairProduction;
} // namespace PROPOSAL

namespace PROPOSAL {
class PhotoPairInterpolant : public CrossSectionInterpolant {
public:
    template <typename T,
        typename = typename enable_if<is_base_of<PhotoPairProduction,
            typename decay<T>::type>::value>::type>
    PhotoPairInterpolant(
        T&&, shared_ptr<const EnergyCutSettings>, const InterpolationDef&);
    virtual ~PhotoPairInterpolant();

    double CalculatedEdx(double energy) override;
    double CalculatedE2dx(double energy) override;
};
} // namespace PROPOSAL

namespace PROPOSAL {
template <typename T,
    typename = typename enable_if<
        is_base_of<PhotoPairProduction, typename decay<T>::type>::value>::type>
PhotoPairInterpolant::PhotoPairInterpolant(
        T&& param, shared_ptr<const EnergyCutSettings> cuts, const InterpolationDef& interpol_def)
    : CrossSectionInterpolant(param, cuts, interpol_def)
{
}
} // namespace PROPOSAL
