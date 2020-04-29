
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

#include "PROPOSAL/crossection/CrossSectionIntegral.h"

namespace PROPOSAL {
class Photonuclear;
} // namespace PROPOSAL

namespace PROPOSAL {
class PhotoIntegral : public CrossSectionIntegral {
public:
    template <typename T,
        typename = typename enable_if<
            is_base_of<Photonuclear, typename decay<T>::type>::value>::type>
    PhotoIntegral(T&&, shared_ptr<const EnergyCutSettings>);
    virtual ~PhotoIntegral();
};
} // namespace PROPOSAL

namespace PROPOSAL {
template <typename T,
    typename = typename enable_if<
        is_base_of<Photonuclear, typename decay<T>::type>::value>::type>
PhotoIntegral::PhotoIntegral(T&& param, shared_ptr<const EnergyCutSettings> cut)
    : CrossSectionIntegral(param, cut)
{
}
} // namespace PROPOSAL
