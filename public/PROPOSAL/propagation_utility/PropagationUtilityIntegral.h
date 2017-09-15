
#pragma once

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/math/Integral.h"

namespace PROPOSAL {

class PropagationUtilityIntegral : public PropagationUtility
{
    public:
    PropagationUtilityIntegral(const ParticleDef&);
    PropagationUtilityIntegral(const ParticleDef&, const Medium&,
                       const EnergyCutSettings&,
                       const Definition& def = Definition());
    PropagationUtilityIntegral(const PropagationUtilityIntegral&);
    virtual ~PropagationUtilityIntegral();

    virtual PropagationUtility* clone() const { return new PropagationUtilityIntegral(*this); }

    double CalculateDisplacement(double ei, double ef, double dist);
    double CalculateFinalEnergy(double ei, double dist);
    double CalculateFinalEnergy(double ei, double rnd, bool particle_interaction);
    double CalculateTrackingIntegal(double initial_energy,
                                    double rnd,
                                    bool particle_interaction);
    double CalculateParticleTime(double ei, double ef);

    private:

    PropagationUtilityIntegral& operator=(const PropagationUtilityIntegral&); // Undefined & not allowed

    Integral integral_;
    Integral prop_interaction_;
    Integral prop_decay_;
    Integral time_particle_;
};

} /* PROPOSAL */
