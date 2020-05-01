
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

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include <array>
#include <functional>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

using std::array;
using std::decay;
using std::enable_if;
using std::function;
using std::is_base_of;
using std::pair;
using std::remove_reference;
using std::shared_ptr;
using std::unique_ptr;
using std::forward;
using std::vector;
using std::placeholders::_1;
using std::placeholders::_2;

namespace PROPOSAL {
class CrossSection {
protected:
    unique_ptr<Parametrization> parametrization_;
    shared_ptr<const EnergyCutSettings> cuts_;

public:
    /* CrossSection() { std::cout <<"default cstr. " << std::endl; }; */
    CrossSection(unique_ptr<Parametrization>&&, shared_ptr<const EnergyCutSettings>);
    virtual ~CrossSection() = default;

    virtual double CalculatedEdx(double) = 0;
    virtual double CalculatedE2dx(double) = 0;
    virtual double CalculatedNdx(double, double = 1) = 0;

    double GetEnergyCut(double energy) const;
    double GetLowerEnergyLimit() const;

    virtual size_t GetHash() const;
};

using CrossSectionList = vector<shared_ptr<CrossSection>>;
} // namespace PROPOSAL
