
#pragma once

#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/EnergyCutSettings.h"

namespace PROPOSAL
{

class Parametrization
{
    public:
    struct Definition
    {
        double multiplier; //!< multiplier to in- or decrease the Epairproduction cross-sections

        bool lpm_effect_enabled;
        double order_of_interpolation;
        std::string path_to_tables; /// Path to interpolation tables
        bool raw;                   /// Determine if output format of interpolation tables is binary or txt.

        Definition();
        ~Definition();
    };

    public:
        Parametrization(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition = Definition());
        Parametrization(const Parametrization&);
        virtual ~Parametrization();

        virtual Parametrization* clone() const = 0;

        //bounds of integration
        struct IntegralLimits
        {
            double vMax; //!< upper bound of integration
            double vUp;  //!< lower bound of integration
            double vMin; //!< lowest physical possible bound of integration
        };

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        virtual double DifferentialCrossSection(double energy, double v) = 0;
        virtual double CalculateParametrization(double energy, double v) = 0;
        virtual double FunctionToDEdxIntegral(double energy, double v) = 0;
        virtual double FunctionToDNdxIntegral(double energy, double v) = 0;
        virtual IntegralLimits GetIntegralLimits(double energy) = 0;

        // ----------------------------------------------------------------- //
        // Getter
        // ----------------------------------------------------------------- //

        const ParticleDef& GetParticleDef() const { return particle_def_; }
        const Medium& GetMedium() const { return *medium_; }
        const EnergyCutSettings& GetEnergyCuts() const { return cut_settings_; }
        double GetMultiplier() const { return param_def_.multiplier; }

        // ----------------------------------------------------------------- //
        // Setter
        // ----------------------------------------------------------------- //

        void SetCurrentComponent(Components::Component* component) {current_component_ = component;}

    protected:
        typedef std::vector<Components::Component*> ComponentVec;

        const ParticleDef particle_def_;
        const Medium* medium_;
        const EnergyCutSettings cut_settings_;

        const Components::Component* current_component_;

        Definition param_def_;
        bool init_lpm_effect_;
        // double multiplier_;
        // bool lpm_effect_enabled_;
};

} /* PROPOSAL */
