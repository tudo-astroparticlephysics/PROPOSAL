
#pragma once

#include <vector>

#include "PROPOSAL/particle/PROPOSALParticle.h"

namespace PROPOSAL
{

namespace Components
{
    class Component;
}

class Parametrization;

class CrossSection
{
    public:
    CrossSection(const DynamicData::Type&, const Parametrization&);
    CrossSection(const CrossSection&);
    virtual ~CrossSection();

    virtual CrossSection* clone() const = 0;

    friend std::ostream& operator<<(std::ostream&, CrossSection const&);

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double CalculatedEdx(double energy) = 0;
    virtual double CalculatedE2dx(double energy) = 0;
    virtual double CalculatedNdx(double energy) = 0;
    virtual double CalculatedNdx(double energy, double rnd) = 0;
    virtual double CalculateStochasticLoss(double energy, double rnd1, double rnd2) = 0;

    // ----------------------------------------------------------------- //
    // Getter
    // ----------------------------------------------------------------- //

    DynamicData::Type GetTypeId() const { return type_id_; }
    Parametrization& GetParametrization() const { return *parametrization_; }

    protected:
    typedef std::vector<Components::Component*> ComponentVec;

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
                 //!ClaculateSochasticLoss when it is already done
};

} /* PROPOSAL */
