
#pragma once

// #include <string>
// #include <vector>

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/particle/ParticleDef.h"
// #include "PROPOSAL/methods.h"
// #include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/crossection/factories/PhotonuclearFactory.h"
#include "PROPOSAL/crossection/factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crossection/factories/EpairProductionFactory.h"
#include "PROPOSAL/crossection/factories/IonizationFactory.h"

namespace PROPOSAL {

// class ContinuousRandomization;
class CrossSection;
class Photonuclear;
class Bremsstrahlung;
class Ionization;
class EpairProduction;

class Medium;
class EnergyCutSettings;
class Geometry;

class Utility
{
    public:
    struct Definition
    {
        BremsstrahlungFactory::Definition brems_def;
        PhotonuclearFactory::Definition photo_def;
        EpairProductionFactory::Definition epair_def;
        IonizationFactory::Definition ioniz_def;

        Definition();
        ~Definition();
    };

    public:
    // Utility(const ParticleDef&);
    Utility(const ParticleDef&,
            const Medium&,
            const EnergyCutSettings&,
            Definition);
    Utility(const ParticleDef&,
            const Medium&,
            const EnergyCutSettings&,
            Definition,
            const InterpolationDef&);
    Utility(const std::vector<CrossSection*>&);
    Utility(const Utility&);
    virtual ~Utility();

    const Definition& GetDefinition() const { return utility_def_; }
    const ParticleDef& GetParticleDef() const { return particle_def_; }
    const Medium& GetMedium() const { return *medium_; }
    const std::vector<CrossSection*>& GetCrosssections() const { return crosssections_; }

    protected:
    Utility& operator=(const Utility&); // Undefined & not allowed

    // --------------------------------------------------------------------- //
    // Protected members
    // --------------------------------------------------------------------- //

    ParticleDef particle_def_;
    Medium* medium_;
    EnergyCutSettings cut_settings_;

    std::vector<CrossSection*> crosssections_;

    Definition utility_def_;
};

class UtilityDecorator
{
    public:
    UtilityDecorator(const Utility&);
    UtilityDecorator(const UtilityDecorator&);
    virtual ~UtilityDecorator();

    virtual UtilityDecorator* clone() const = 0;

    virtual double FunctionToIntegral(double energy);
    virtual double Calculate(double ei, double ef, double rnd) = 0;
    virtual double GetUpperLimit(double ei, double rnd) = 0;

    const Utility& GetUtility() const { return utility_; }

    protected:
    UtilityDecorator& operator=(const UtilityDecorator&); // Undefined & not allowed

    const Utility& utility_;
};


// #<{(|! \class ProcessSector ProcessSector.h "CrossSections.h"
//     \brief initializes all cross sections and keeps references to them
//  |)}>#
// class PropagationUtility
// {
//     public:
//
//     struct Definition
//     {
//         double brems_multiplier; //!< multiplier to in- or decrease the Bremsstrahlung cross-sections
//         double photo_multiplier; //!< multiplier to in- or decrease the Photonucler cross-sections
//         double ioniz_multiplier; //!< multiplier to in- or decrease the Ionization cross-sections
//         double epair_multiplier; //!< multiplier to in- or decrease the Epairproduction cross-sections
//
//         PhotonuclearFactory::Enum photo_parametrization;
//         PhotonuclearFactory::Shadow photo_shadow;
//         bool hardbb_enabled;
//
//         BremsstrahlungFactory::Enum brems_parametrization;
//
//         bool lpm_effect_enabled;
//
//         bool do_exact_time_calculation;   //!< exact local time calculation enabled if true
//         bool do_continuous_randomization; //!< exact local time calculation enabled if true
//
//         double order_of_interpolation;
//         bool raw;                   /// Determine if output format of interpolation tables is binary or txt.
//         std::string path_to_tables; /// Path to interpolation tables
//
//         Definition();
//         ~Definition();
//     };
//
//     public:
//     PropagationUtility(const ParticleDef&);
//     PropagationUtility(const ParticleDef&, const Medium&, const EnergyCutSettings&, const Definition& def = Definition());
//     PropagationUtility(const PropagationUtility&);
//     virtual ~PropagationUtility();
//
//     virtual PropagationUtility* clone() const = 0; // virtual constructor idiom (used for deep copies)
//
//     // bool operator==(const PropagationUtility& collection) const;
//     // bool operator!=(const PropagationUtility& collection) const;
//     // friend std::ostream& operator<<(std::ostream& os, PropagationUtility const& collection);
//
//     // --------------------------------------------------------------------- //
//     // Member functions
//     // --------------------------------------------------------------------- //
//
//
//     #<{(|*
//      *  Makes Decay
//      *
//      *  \return pair of energy loss [MeV] and kind of interaction
//      |)}>#
//     double MakeDecay(double energy);
//
//
//     #<{(|*
//      *  Randomize the continous energy losses
//      *
//      *  \return randomized energy
//      |)}>#
//     double Randomize(double initial_energy, double final_energy, double rnd);
//
//     #<{(|!
//     * returns the value of the distance integral from ei to ef;
//     * \f[ dx = \int_{e_i}^{e_f} - \frac{1}{  \frac{dE}{dx}\big|_{Ioniz}
//     * +\frac{dE}{dx}\big|_{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}
//     * \big|_{epair}  } dE \f]
//     *
//     * \param ei lower integration limit (initial energy)
//     * \param ef upper integration limit (final energy)
//     * \param dist ???
//     *
//     * \return value of the distance integral from ei to ef
//     |)}>#
//     virtual double CalculateDisplacement(double ei, double ef, double dist) = 0;
//
//     #<{(|!
//     * final energy, corresponding to the given value of displacement dist;
//     * returns \f$e_f\f$ which fullfills
//     * \f[\int_{e_i}^{\hat{e}_f} - \frac{1}{  \frac{dE}{dx}\big|_{Ioniz} +
//     * \frac{dE}{dx}\big|_{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}
//     * \big|_{epair}  } dE = dist\f]
//     *
//     * \param ei initial energy [MeV]
//     * \param dist value of displacement
//     *
//     * \return final energy [MeV]
//     |)}>#
//
//     virtual double CalculateFinalEnergy(double ei, double dist) = 0;
//
//     #<{(|*
//      * final energy, corresponding to the given value rnd of the
//      * tracking integral
//      *
//      *  \param  ei                      initial energy
//      *  \param  rnd                     random number which is used for the calculation
//      *  \param  particle_interaction    particle interaction? (false = decay)
//      *  \return final energy due to continous energy losses [MeV]
//      |)}>#
//
//     virtual double CalculateFinalEnergy(double ei, double rnd, bool particle_interaction) = 0;
//
//     #<{(|*
//      * Calculates the value of the tracking integral
//      *
//      *  \param  initial_energy          initial energy
//      *  \param  rnd                     random number which is used for the calculation
//      *  \param  particle_interaction    particle interaction? (false = decay)
//      *  \return value of the tracking integral [ 1 ]
//      |)}>#
//
//     virtual double CalculateTrackingIntegal(double initial_energy, double rnd, bool particle_interaction) = 0;
//
//     #<{(|!
//     * time delta, corresponding to the given propagation distance
//     *
//     * \param    ei  initial energy
//     * \param    ef  final energy
//     * \return   time delta
//     |)}>#
//
//     virtual double CalculateParticleTime(double ei, double ef) = 0;
//
//     #<{(|!
//     * Quantitiy need for randomization
//     *
//     * \param    ei  initial energy
//     * \param    ef  final energy
//     * \return   \int_{E_i}^{E_f} \frac{d^2E}{d\epsilon^2} dE
//     |)}>#
//     virtual double CalculateDE2de(double ei, double ef) = 0;
//
//     #<{(|!
//     * function for range calculation for given energy - interface to Integral;
//     * \f[f(E) =- \frac{1}{ \frac{dE}{dx}\big|_{Ioniz} +\frac{dE}{dx}\big|
//     * _{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}\big|_{epair}  }\f]
//     * \return range calculation for given energy [cm/E]
//     * \param E energy [MeV]
//     |)}>#
//
//     virtual double FunctionToIntegral(double energy);
//
//     #<{(|*
//      * function for energy range calculation - interface to Integral
//      * For decay
//      *
//      *  \param  energy particle energy
//      *  \return Returns the probability [1/MeV]
//      |)}>#
//
//     virtual double FunctionToPropIntegralDecay(double energy);
//
//     //----------------------------------------------------------------------------//
//
//     #<{(|*
//      * function for energy range calculation - interface to Integral
//      * For Interaction
//      *
//      *  \param  energy particle energy
//      *  \return Returns the probability [1/MeV]
//      |)}>#
//
//     virtual double FunctionToPropIntegralInteraction(double energy);
//
//
//     #<{(|!
//     * function for time delta calculation - interface to Integral
//     *
//     |)}>#
//     virtual double FunctionToTimeIntegral(double energy);
//
//
//     #<{(|*
//      * function for continous randomization calculation - interface to Integral
//      *
//      *  \param  energy particle energy
//      *  \return //TODO(mario):  Sat 2017/09/16
//      |)}>#
//     double FunctionToDE2deIntegral(double energy);
//
//
//     // --------------------------------------------------------------------- //
//     // Getter
//     // --------------------------------------------------------------------- //
//
//     double GetIni() const { return ini_; }
//
//     const Medium& GetMedium() const { return *medium_; }
//     const ParticleDef& GetParticleDef() const { return particle_def_; }
//     const EnergyCutSettings& GetCutSettings() const { return cut_settings_; }
//     const std::vector<CrossSection*>& GetCrosssections() const { return crosssections_; }
//
//     protected:
//     PropagationUtility& operator=(const PropagationUtility&); // Undefined & not allowed
//
//     // --------------------------------------------------------------------- //
//     // Protected members
//     // --------------------------------------------------------------------- //
//
//     // Just a temporary to store -> bad design
//     double ini_;
//
//     Definition utility_def_;
//
//     ParticleDef particle_def_;
//     Medium* medium_;
//     EnergyCutSettings cut_settings_;
//
//     std::vector<CrossSection*> crosssections_;
// };
}
