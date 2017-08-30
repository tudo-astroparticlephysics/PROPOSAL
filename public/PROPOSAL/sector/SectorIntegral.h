
#pragma once

#include "PROPOSAL/sector/Sector.h"
#include "PROPOSAL/math/Integral.h"

namespace PROPOSAL {

class SectorIntegral : public Sector
{
    public:
    SectorIntegral();
    SectorIntegral(const Medium&,
                       const Geometry&,
                       const EnergyCutSettings&,
                       const SectorDef& def = SectorDef());
    SectorIntegral(const SectorIntegral&);
    virtual ~SectorIntegral();

    virtual Sector* clone() const { return new SectorIntegral(*this); }

    double CalculateDisplacement(const PROPOSALParticle&, double ei, double ef, double dist);
    double CalculateFinalEnergy(const PROPOSALParticle&, double ei, double dist);
    double CalculateFinalEnergy(const PROPOSALParticle&, double ei, double rnd, bool particle_interaction);
    double CalculateTrackingIntegal(const PROPOSALParticle&,
                                    double initial_energy,
                                    double rnd,
                                    bool particle_interaction);
    double CalculateParticleTime(const PROPOSALParticle&, double ei, double ef);

    private:

    SectorIntegral& operator=(const SectorIntegral&); // Undefined & not allowed

    Integral integral_;
    Integral prop_interaction_;
    Integral prop_decay_;
    Integral time_particle_;
};

} /* PROPOSAL */
