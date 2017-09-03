
#pragma once

#include "PROPOSAL/sector/Sector.h"
#include "PROPOSAL/math/Interpolant.h"

namespace PROPOSAL {

class SectorInterpolant : public Sector
{
    public:
    SectorInterpolant(PROPOSALParticle&);
    SectorInterpolant(PROPOSALParticle&, const Medium&,
                          const Geometry&,
                          const EnergyCutSettings&,
                          const Definition& def = Definition());
    SectorInterpolant(const SectorInterpolant&);
    virtual ~SectorInterpolant();

    virtual Sector* clone() const { return new SectorInterpolant(*this); }


    double CalculateDisplacement( double ei, double ef, double dist);
    double CalculateFinalEnergy( double ei, double dist);
    double CalculateFinalEnergy( double ei, double rnd, bool particle_interaction);
    double CalculateTrackingIntegal(
                                    double initial_energy,
                                    double rnd,
                                    bool particle_interaction);
    double CalculateParticleTime( double ei, double ef);

    private:

    SectorInterpolant& operator=(const SectorInterpolant&); // Undefined & not allowed

    // --------------------------------------------------------------------- //
    // Private methods
    // --------------------------------------------------------------------- //

    void InitInterpolation( std::string filename, bool raw);
    void InitTimeInterpolation( std::string filename, bool raw);

    double FunctionToBuildInterpolant( double);
    double InterpolPropDecay( double);
    double InterpolPropInteraction( double);
    double InterpolTimeParticleDiff( double);

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
