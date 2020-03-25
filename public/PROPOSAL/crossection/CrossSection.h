
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
#include <memory>
#include <utility>

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/EnergyCutSettings.h"

namespace PROPOSAL {

using Function = std::function<double(double)>;

namespace Components {
class Component;
}

class Parametrization;


class CrossSection
{
public:
    CrossSection(const Parametrization&, std::shared_ptr<const EnergyCutSettings>);
    CrossSection(const CrossSection&);
    virtual ~CrossSection();

    bool operator==(const CrossSection& cross_section) const;
    bool operator!=(const CrossSection& cross_section) const;

    //virtual CrossSection* clone() const = 0;

    friend std::ostream& operator<<(std::ostream&, CrossSection const&);

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double CalculatedEdx(double energy)                                     = 0;
    virtual double CalculatedE2dx(double energy)                                    = 0;
    virtual double CalculatedNdx(double energy)                                     = 0;
    virtual double CalculatedNdx(double energy, double rnd)                         = 0;
    virtual double CalculateStochasticLoss(double energy, double rnd1, double rnd2) = 0;

    virtual double GetEnergyCut(double energy);

    // CalculateProducedParticles Return values:
    // First Parameter: List of produced particles by stochastic interaction (default: no particles, e.g. empty list)
    // Second parameter: Is the interaction a fatal interaction (e.g. will the initial particle vanish after interaction?)
    virtual std::pair<std::vector<DynamicData>, bool> CalculateProducedParticles(
            double energy, double energy_loss, const Vector3D& initial_direction){
        (void)energy; (void)energy_loss; (void)initial_direction; return std::make_pair(std::vector<DynamicData>(), false);
    }

    virtual std::pair<double, double> StochasticDeflection(double energy, double energy_loss);

    virtual double CalculateCumulativeCrossSection(double energy, int component, double v) = 0;

    // ----------------------------------------------------------------- //
    // Getter
    // ----------------------------------------------------------------- //

    int GetTypeId() const { return static_cast<int>(type_id_); }
    Parametrization& GetParametrization() const { return *parametrization_; }

protected:

    virtual bool compare(const CrossSection&) const = 0;

    // ----------------------------------------------------------------- //
    // Protected methods
    // ----------------------------------------------------------------- //

    virtual double CalculateStochasticLoss(double energy, double rnd1) = 0;

    // ----------------------------------------------------------------- //
    // Protected member
    // ----------------------------------------------------------------- //

    InteractionType type_id_;

    Parametrization* parametrization_;

    std::vector<double> prob_for_component_; //!< probability for each medium component to
                                             //!< interact with the particle (formerly h_)
    double sum_of_rates_;

    const std::vector<Components::Component>& components_;

    double rnd_; //!< This random number will be stored in CalculateDNdx to avoid calculate dNdx a second time in
                 //! ClaculateSochasticLoss when it is already done

    std::shared_ptr<const EnergyCutSettings> cuts_;
};

class CrossSectionBuilder final : public CrossSection
{
public:
    CrossSectionBuilder() : CrossSection(*param, nullptr){
        // Set zero return functions as default
        dEdx_function = [](double x) { (void)x; return 0;};
        dE2dx_function = [](double x) { (void)x; return 0;};
        dNdx_function = [](double x) { (void)x; return 0;};
        dNdx_rnd_function = [](double x1, double x2) { (void)x1; (void)x2; return 0;};
        StochasticLoss_function = [](double x1, double x2, double x3) { (void)x1; (void)x2; (void)x3; return 0;};
        CumulativeCrossSection_function = [](double x1, double x2, double x3) { (void)x1; (void)x2; (void)x3; return 0;};
    }

    static Parametrization* param;

    double CalculatedEdx(double energy) override;
    double CalculatedE2dx(double energy) override;
    double CalculatedNdx(double energy) override;
    double CalculatedNdx(double energy, double rnd) override;
    double CalculateStochasticLoss(double energy, double rnd1, double rnd2) override;
    double CalculateCumulativeCrossSection(double energy, int component, double v) override;

    void SetdEdx_function(Function func);
    void SetdE2dx_function(Function func);
    void SetdNdx_function(Function func);
    void SetdNdx_rnd_function(std::function<double(double, double)> func);
    void SetStochasticLoss_function(std::function<double(double, double, double)> func);
    void SetCumulativeCrossSection_function(std::function<double(double, double, double)> func);

protected:
    virtual bool compare(const CrossSection& cross) const override;
    double CalculateStochasticLoss(double energy, double rnd1) override;

private:
    Function dEdx_function;
    Function dE2dx_function;
    Function dNdx_function;
    std::function<double(double, double)> dNdx_rnd_function;
    std::function<double(double, double, double)> StochasticLoss_function;
    std::function<double(double, double, double)> CumulativeCrossSection_function;
};

std::ostream& operator<<(std::ostream&, PROPOSAL::CrossSection const&);
using CrossSectionList = std::vector<std::shared_ptr<CrossSection>>;

} // namespace PROPOSAL
