
#pragma once

// #include <string>
// #include <vector>

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/crossection/PhotonuclearFactory.h"
#include "PROPOSAL/crossection/BremsstrahlungFactory.h"

namespace PROPOSAL {

// class ContinuousRandomization;
class CrossSection;
class Medium;
class EnergyCutSettings;
class Geometry;

/*! \class ProcessSector ProcessSector.h "CrossSections.h"
    \brief initializes all cross sections and keeps references to them
 */
class Sector
{
    public:
    struct ParticleLocation
    {
        enum Enum
        {
            InfrontDetector = 0,
            InsideDetector,
            BehindDetector
        };
    };

    struct Definition
    {
        bool do_weighting;      //!< Do weigthing? Set to false in constructor
        double weighting_order; //!< Re-weighting order. Set to 0 in constructor

        double brems_multiplier; //!< multiplier to in- or decrease the Bremsstrahlung cross-sections
        double photo_multiplier; //!< multiplier to in- or decrease the Photonucler cross-sections
        double ioniz_multiplier; //!< multiplier to in- or decrease the Ionization cross-sections
        double epair_multiplier; //!< multiplier to in- or decrease the Epairproduction cross-sections

        PhotonuclearFactory::Enum photo_parametrization;
        PhotonuclearFactory::Shadow photo_shadow;
        bool hardbb_enabled;

        BremsstrahlungFactory::Enum brems_parametrization;

        // TODO(mario): must be removed Fri 2017/08/25
        bool do_scattering;                       //!< if true moliere scattering is enabled
        ScatteringFactory::Enum scattering_model; //!< if true moliere scattering is enabled

        bool do_continuous_randomization; //!< exact local time calculation enabled if true
        bool do_exact_time_calculation;   //!< exact local time calculation enabled if true
        bool lpm_effect_enabled;

        // int location;              //!< 0 = infront of the detector, 1 = inside the detector, 2 = behind the detector
        Sector::ParticleLocation::Enum
            location; //!< 0 = infront of the detector, 1 = inside the detector, 2 = behind the detector

        double order_of_interpolation;
        bool raw;                   /// Determine if output format of interpolation tables is binary or txt.
        std::string path_to_tables; /// Path to interpolation tables

        Definition();
        ~Definition();
    };

    public:
    Sector(PROPOSALParticle&);
    Sector(PROPOSALParticle&, const Medium&, const Geometry&, const EnergyCutSettings&, const Definition& def = Definition());
    Sector(const Sector&);
    virtual ~Sector();

    virtual Sector* clone() const = 0; // virtual constructor idiom (used for deep copies)

    // Sector& operator=(const Sector& collection);
    // bool operator==(const Sector& collection) const;
    // bool operator!=(const Sector& collection) const;
    // friend std::ostream& operator<<(std::ostream& os, Sector const& collection);

    // --------------------------------------------------------------------- //
    // Member functions
    // --------------------------------------------------------------------- //

    /**
     * Propagates the particle of initial energy e to the distance r.
     * Returns the final energy if the
     * particle has survived or the track length to the
     * point of disappearance with a minus sign otherwise.
     *
     *  \param  distance   maximum track length
     *  \param  energy   initial energy
     *  \return energy at distance OR -(track length)
     */

    virtual double Propagate(double distance);

    /**
     * Calculates the contiuous loss till the first stochastic loss happend
     * and subtract it from initial energy
     * Also caluclate the energy at which the particle decay
     * These to energys can be compared to decide if a decay or particle interaction
     * happens
     *
     *  \param  initial_energy   initial energy
     *  \return pair.first final energy befor first interaction pair.second decay energy at which the
     *          particle decay
     */
    virtual std::pair<double, double> CalculateEnergyTillStochastic(double initial_energy);

    /*!
    * advances the particle by the given distance
    * Sets the x,y and z coordinates of particle_
    * and its time and propagation distance
    *
    * \param    dr  flight distance
    * \param    ei  initial energy
    * \param    ef  final energy
    */

    void AdvanceParticle(double dr, double ei, double ef);

    /**
     *  Makes Stochastic Energyloss
     *
     *  \return pair of energy loss [MeV] and kind of interaction
     */
    virtual std::pair<double, DynamicData::Type> MakeStochasticLoss();

    /**
     *  Makes Decay
     *
     *  \return pair of energy loss [MeV] and kind of interaction
     */
    double MakeDecay(double energy);

    /*!
    returns the value of the distance integral from ei to ef;
    \f[ dx = \int_{e_i}^{e_f} - \frac{1}{  \frac{dE}{dx}\big|_{Ioniz}
    +\frac{dE}{dx}\big|_{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}
    \big|_{epair}  } dE \f]
    \return value of the distance integral from ei to ef
    \param ei lower integration limit (initial energy)
    \param ef upper integration limit (final energy)
    \param dist ???
    */
    virtual double CalculateDisplacement(double ei, double ef, double dist) = 0;

    /*!
    final energy, corresponding to the given value of displacement dist;
    returns \f$e_f\f$ which fullfills
    \f[\int_{e_i}^{\hat{e}_f} - \frac{1}{  \frac{dE}{dx}\big|_{Ioniz} +
    \frac{dE}{dx}\big|_{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}
    \big|_{epair}  } dE = dist\f]
    \return final energy [MeV]
    \param ei initial energy [MeV]
    \param dist value of displacement
    */

    virtual double CalculateFinalEnergy(double ei, double dist) = 0;

    /**
     * final energy, corresponding to the given value rnd of the
     * tracking integral
     *
     *  \param  ei                      initial energy
     *  \param  rnd                     random number which is used for the calculation
     *  \param  particle_interaction    particle interaction? (false = decay)
     *  \return final energy due to continous energy losses [MeV]
     */

    virtual double CalculateFinalEnergy(double ei, double rnd, bool particle_interaction) = 0;

    /**
     * Calculates the value of the tracking integral
     *
     *  \param  initial_energy          initial energy
     *  \param  rnd                     random number which is used for the calculation
     *  \param  particle_interaction    particle interaction? (false = decay)
     *  \return value of the tracking integral [ 1 ]
     */

    virtual double CalculateTrackingIntegal(double initial_energy, double rnd, bool particle_interaction) = 0;

    /*!
    * time delta, corresponding to the given propagation distance
    *
    * \param    ei  initial energy
    * \param    ef  final energy
    * \return   time delta
    */

    virtual double CalculateParticleTime(double ei, double ef) = 0;

    /*!
    * function for range calculation for given energy - interface to Integral;
    * \f[f(E) =- \frac{1}{ \frac{dE}{dx}\big|_{Ioniz} +\frac{dE}{dx}\big|
    * _{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}\big|_{epair}  }\f]
    * \return range calculation for given energy [cm/E]
    * \param E energy [MeV]
    */

    virtual double FunctionToIntegral(double energy);

    /**
     * randomize the continuous energy loss
     *
     *  \param  initial_energy          initial energy
     *  \param  final_energy            final energy
     */

    // virtual double Randomize(double initial_energy, double final_energy);

    // --------------------------------------------------------------------- //
    // Enable options & Setter
    // --------------------------------------------------------------------- //

    void EnableLpmEffect();
    void DisableLpmEffect();
    // void EnableContinuousRandomization();
    // void DisableContinuousRandomization();
    void EnableExactTimeCalculation();
    void DisableExactTimeCalculation();
    void SetLocation(ParticleLocation::Enum location) { collection_def_.location = location; }

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    ParticleLocation::Enum GetLocation() const { return collection_def_.location; }
    double GetIni() const { return ini_; }

    bool GetDoScattering() const { return collection_def_.do_scattering; }
    bool GetLpmEffectEnabled() const { return collection_def_.lpm_effect_enabled; }
    // bool GetDoRandomization() const { return do_continuous_randomization_; }
    // bool GetEnableRandomization() const { return enable_randomization_; }

    Scattering* GetScattering() const { return scattering_; }
    Medium* GetMedium() const { return medium_; }
    PROPOSALParticle& GetParticle() const { return particle_; }
    Geometry* GetGeometry() const { return geometry_; }
    const EnergyCutSettings& GetCutSettings() const { return cut_settings_; }
    std::vector<CrossSection*> GetCrosssections() const { return crosssections_; }
    // ContinuousRandomization* GetContinuousRandomization() const { return randomizer_; }

    protected:
    Sector& operator=(const Sector&); // Undefined & not allowed

    // --------------------------------------------------------------------- //
    // Protected members
    // --------------------------------------------------------------------- //

    // Just a temporary to store -> bad design
    double ini_;

    Definition collection_def_;

    // TODO(mario): Do better weight enabling Fri 2017/08/25
    // bool do_weighting_;          //!< Do weigthing? Set to false in constructor
    // double weighting_order_;     //!< Re-weighting order. Set to 0 in constructor
    double weighting_starts_at_; //!< Distance at which re-weighting starts. Set to 0 in constructor
    //
    // // bool enable_randomization_; //!< if true continuous randomization will be enabled (to remember if
    // randomization
    //                             //!should be enable when cross sections are initalized)
    //
    // // //TODO(mario): must be removed Fri 2017/08/25
    // bool do_scattering_;               //!< if true moliere scattering is enabled
    // bool lpm_effect_enabled_;
    // bool do_exact_time_calculation_; //!< exact local time calculation enabled if true
    //
    // int location_; //!< 0 = infront of the detector, 1 = inside the detector, 2 = behind the detector
    // double density_correction_; //!< density correction factor

    PROPOSALParticle& particle_;
    Geometry* geometry_;
    Medium* medium_;
    EnergyCutSettings cut_settings_;

    // ContinuousRandomization* randomizer_;
    Scattering* scattering_;

    std::vector<CrossSection*> crosssections_;

    // --------------------------------------------------------------------- //
    // Protected member functions
    // --------------------------------------------------------------------- //

    /**
     * function for energy range calculation - interface to Integral
     * For decay
     *
     *  \param  energy particle energy
     *  \return Returns the probability [1/MeV]
     */

    virtual double FunctionToPropIntegralDecay(double energy);

    //----------------------------------------------------------------------------//

    /**
     * function for energy range calculation - interface to Integral
     * For Interaction
     *
     *  \param  energy particle energy
     *  \return Returns the probability [1/MeV]
     */

    virtual double FunctionToPropIntegralInteraction(double energy);

    //----------------------------------------------------------------------------//

    /*!
    * function for time delta calculation - interface to Integral
    *
    */

    virtual double FunctionToTimeIntegral(double energy);
};
}
