
#pragma once

#include "PROPOSAL/sector/Collection.h"
#include "PROPOSAL/math/Interpolant.h"

namespace PROPOSAL {

class CollectionInterpolant : public Collection
{
    public:
    CollectionInterpolant();
    CollectionInterpolant(const Medium&,
                          const Geometry&,
                          const EnergyCutSettings&,
                          const CollectionDef& def = CollectionDef());
    CollectionInterpolant(const CollectionInterpolant&);
    virtual ~CollectionInterpolant();

    virtual Collection* clone() const { return new CollectionInterpolant(*this); }


    double Propagate(PROPOSALParticle&, double distance);
    double CalculateDisplacement(const PROPOSALParticle&, double ei, double ef, double dist);
    double CalculateFinalEnergy(const PROPOSALParticle&, double ei, double dist);
    double CalculateFinalEnergy(const PROPOSALParticle&, double ei, double rnd, bool particle_interaction);
    double CalculateTrackingIntegal(const PROPOSALParticle&,
                                    double initial_energy,
                                    double rnd,
                                    bool particle_interaction);
    double CalculateParticleTime(const PROPOSALParticle&, double ei, double ef);

    private:

    CollectionInterpolant& operator=(const CollectionInterpolant&); // Undefined & not allowed

    // --------------------------------------------------------------------- //
    // Private methods
    // --------------------------------------------------------------------- //

    void InitInterpolation(const PROPOSALParticle&, std::string filename, bool raw);
    void InitTimeInterpolation(const PROPOSALParticle&, std::string filename, bool raw);

    double FunctionToBuildInterpolant(const PROPOSALParticle&, double);
    double InterpolPropDecay(const PROPOSALParticle&, double);
    double InterpolPropInteraction(const PROPOSALParticle& particle, double);
    double InterpolTimeParticleDiff(const PROPOSALParticle&, double);

    // --------------------------------------------------------------------- //
    // Private members
    // --------------------------------------------------------------------- //

    bool initialized_interpolation_; /// Stores if Interplation tables are already created or loaded.

    // ----------------------------------------------------------------------------
    /// @brief indicates if the interpolated function is increasing or decreasing.
    ///
    /// True  =   Interpolated function is increasing. \n
    /// False =   Interpolated function is decreasing. BigLow is set to f(x_min)
    ///           and the interpolation points that are saved are BigLow-f(x).
    // ----------------------------------------------------------------------------
    bool up_;

    double big_low_interatction_; //!< See the describtion of the variable up.
    double big_low_decay_;
    double store_dif_interaction_; //!< Stores the interpolated values for farther calculations
    double store_dif_decay_;

    Interpolant* interpolant_;
    Interpolant* interpolant_diff_;

    Interpolant* interpol_time_particle_;
    Interpolant* interpol_time_particle_diff_;

    Interpolant* interpol_prop_decay_;           //!< Interpolant object of the Integral of the function FunctionToPropIntegralDecay
    Interpolant* interpol_prop_decay_diff_;      //!< Interpolant object of the function FunctionToPropIntegralDecay
    Interpolant* interpol_prop_interaction_;     //!< Interpolant object of the Integral of the function FunctionToPropIntegralInteraction
    Interpolant* interpol_prop_interaction_diff_;//!< Interpolant object of the function FunctionToPropIntegralInteraction

};

} /* PROPOSAL */
