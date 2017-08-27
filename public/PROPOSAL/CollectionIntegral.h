
#pragma once

#include "PROPOSAL/Collection.h"
#include "PROPOSAL/Integral.h"

namespace PROPOSAL {

class CollectionIntegral : public Collection
{
    public:
    CollectionIntegral();
    CollectionIntegral(const Medium&,
                       const Geometry&,
                       const EnergyCutSettings&,
                       const CollectionDef& def = CollectionDef());
    CollectionIntegral(const CollectionIntegral&);
    virtual ~CollectionIntegral();

    virtual Collection* clone() const { return new CollectionIntegral(*this); }

    double CalculateDisplacement(const PROPOSALParticle&, double ei, double ef, double dist);
    double CalculateFinalEnergy(const PROPOSALParticle&, double ei, double dist);
    double CalculateFinalEnergy(const PROPOSALParticle&, double ei, double rnd, bool particle_interaction);
    double CalculateTrackingIntegal(const PROPOSALParticle&,
                                    double initial_energy,
                                    double rnd,
                                    bool particle_interaction);
    double CalculateParticleTime(const PROPOSALParticle&, double ei, double ef);

    private:

    CollectionIntegral& operator=(const CollectionIntegral&); // Undefined & not allowed

    Integral integral_;
    Integral prop_interaction_;
    Integral prop_decay_;
    Integral time_particle_;
};

} /* PROPOSAL */
