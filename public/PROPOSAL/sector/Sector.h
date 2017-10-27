
#pragma once

// #include <string>
// #include <vector>

#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

namespace PROPOSAL {

class ContinuousRandomizer;
// class CrossSection;
// class Medium;
// class EnergyCutSettings;
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
        bool stopping_decay; //!< Let particle decay if elow is reached but no decay was sampled

        bool do_continuous_randomization;
        bool do_exact_time_calculation;

        ScatteringFactory::Enum scattering_model;

        Sector::ParticleLocation::Enum location;

        Utility::Definition utility_def;

        Definition();
        ~Definition();
    };

    public:
    // Sector(Particle&);
    Sector(Particle&, const Medium&, const EnergyCutSettings&, const Geometry&, const Definition&);
    Sector(Particle&, const Medium&, const EnergyCutSettings&, const Geometry&, const Definition&, const InterpolationDef&);
    Sector(Particle&, const Sector&);
    // Sector(Particle&, const Geometry&, const Utility&, const Scattering&, bool do_interpolation, const Definition& def = Definition());
    Sector(const Sector&);
    virtual ~Sector();

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

    // --------------------------------------------------------------------- //
    // Enable options & Setter
    // --------------------------------------------------------------------- //

    void SetLocation(ParticleLocation::Enum location) { sector_def_.location = location; }

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    ParticleLocation::Enum GetLocation() const { return sector_def_.location; }
    // bool GetDoRandomization() const { return do_continuous_randomization_; }
    // bool GetEnableRandomization() const { return enable_randomization_; }

    Scattering* GetScattering() const { return scattering_; }
    Particle& GetParticle() const { return particle_; }
    Geometry* GetGeometry() const { return geometry_; }
    const Medium* GetMedium() const { return &utility_.GetMedium(); }
    // ContinuousRandomization* GetContinuousRandomization() const { return randomizer_; }

    protected:
    Sector& operator=(const Sector&); // Undefined & not allowed

    // --------------------------------------------------------------------- //
    // Protected members
    // --------------------------------------------------------------------- //

    Definition sector_def_;

    // TODO(mario): Do better weight enabling Fri 2017/08/25
    double weighting_starts_at_; //!< Distance at which re-weighting starts. Set to 0 in constructor

    Particle& particle_;
    Geometry* geometry_;

    Utility utility_;
    UtilityDecorator* displacement_calculator_;
    UtilityDecorator* interaction_calculator_;
    UtilityDecorator* decay_calculator_;
    UtilityDecorator* exact_time_calculator_;

    ContinuousRandomizer* cont_rand_;
    Scattering* scattering_;
};
}
