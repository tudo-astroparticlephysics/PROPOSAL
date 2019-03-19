
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

#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL {

namespace Components {
class Component;
}

class Parametrization;

class CrossSection
{
public:
    CrossSection(const DynamicData::Type&, const Parametrization&);
    CrossSection(const CrossSection&);
    virtual ~CrossSection();

    bool operator==(const CrossSection& cross_section) const;
    bool operator!=(const CrossSection& cross_section) const;

    virtual CrossSection* clone() const = 0;

    friend std::ostream& operator<<(std::ostream&, CrossSection const&);

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double CalculatedEdx(double energy)                                     = 0;
    virtual double CalculatedE2dx(double energy)                                    = 0;
    virtual double CalculatedNdx(double energy)                                     = 0;
    virtual double CalculatedNdx(double energy, double rnd)                         = 0;
    virtual double CalculateStochasticLoss(double energy, double rnd1, double rnd2) = 0;
    virtual std::vector<Particle*> CalculateProducedParticles(double energy, double energy_loss, double rnd1, double rnd2){
        (void)energy; (void)energy_loss; (void)rnd1; (void)rnd2; return std::vector<Particle*>();}

    // ----------------------------------------------------------------- //
    // Getter
    // ----------------------------------------------------------------- //

    DynamicData::Type GetTypeId() const { return type_id_; }
    Parametrization& GetParametrization() const { return *parametrization_; }

protected:
    typedef std::vector<Components::Component*> ComponentVec;

    virtual bool compare(const CrossSection&) const = 0;

    // ----------------------------------------------------------------- //
    // Protected methods
    // ----------------------------------------------------------------- //

    virtual double CalculateStochasticLoss(double energy, double rnd1) = 0;

    // ----------------------------------------------------------------- //
    // Protected member
    // ----------------------------------------------------------------- //

    const DynamicData::Type type_id_;

    Parametrization* parametrization_;

    std::vector<double> prob_for_component_; //!< probability for each medium component to
                                             //!< interact with the particle (formerly h_)
    double sum_of_rates_;

    const ComponentVec& components_;

    double rnd_; //!< This random number will be stored in CalculateDNdx to avoid calculate dNdx a second time in
                 //! ClaculateSochasticLoss when it is already done
};

std::ostream& operator<<(std::ostream&, PROPOSAL::CrossSection const&);

} // namespace PROPOSAL
