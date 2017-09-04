
#pragma once

#include "PROPOSAL/sector/Sector.h"
#include "PROPOSAL/math/Integral.h"

namespace PROPOSAL {

class SectorIntegral : public Sector
{
    public:
    SectorIntegral(PROPOSALParticle&);
    SectorIntegral(PROPOSALParticle&, const Medium&,
                       const Geometry&,
                       const EnergyCutSettings&,
                       const Definition& def = Definition());
    SectorIntegral(const SectorIntegral&);
    virtual ~SectorIntegral();

    virtual Sector* clone() const { return new SectorIntegral(*this); }

    double CalculateDisplacement(double ei, double ef, double dist);
    double CalculateFinalEnergy(double ei, double dist);
    double CalculateFinalEnergy(double ei, double rnd, bool particle_interaction);
    double CalculateTrackingIntegal(double initial_energy,
                                    double rnd,
                                    bool particle_interaction);
    double CalculateParticleTime(double ei, double ef);

    private:

    SectorIntegral& operator=(const SectorIntegral&); // Undefined & not allowed

    Integral integral_;
    Integral prop_interaction_;
    Integral prop_decay_;
    Integral time_particle_;
};

} /* PROPOSAL */
