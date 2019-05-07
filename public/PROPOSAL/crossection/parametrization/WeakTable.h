
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2019 TU Dortmund University, Department of Physics,          *
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
#include <string>
#include <vector>
namespace PROPOSAL {

    const extern std::vector<double> energies;

    const extern std::vector< std::vector<double> > y_nu_p;
    const extern std::vector< std::vector<double> > y_nubar_p;
    const extern std::vector< std::vector<double> > y_nu_n;
    const extern std::vector< std::vector<double> > y_nubar_n;

    const extern std::vector< std::vector<double> > sigma_nu_p;
    const extern std::vector< std::vector<double> > sigma_nubar_p;
    const extern std::vector< std::vector<double> > sigma_nu_n;
    const extern std::vector< std::vector<double> > sigma_nubar_n;

} // namespace PROPOSAL
