
#pragma once

// #include <vector>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL
{

class PROPOSALParticle;

class ManyBodyPhaseSpace : public DecayChannel
{
    public:
    ManyBodyPhaseSpace(std::vector<ParticleDef> daughters);
    ManyBodyPhaseSpace(const ManyBodyPhaseSpace& mode);
    virtual ~ManyBodyPhaseSpace();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() { return new ManyBodyPhaseSpace(*this); }

    // ----------------------------------------------------------------------------
    /// @brief Many body phase space decay
    ///
    /// Calculate decay products with the help of the Raubold Lynch algorithm.
    ///
    /// @param PROPOSALParticle
    ///
    /// @return Vector of particles, the decay products
    // ----------------------------------------------------------------------------
    DecayProducts Decay(PROPOSALParticle*);

    const std::string& GetName() const { return name_; }

    private:
    ManyBodyPhaseSpace& operator=(const ManyBodyPhaseSpace&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;
    void print(std::ostream&) const;

    double Momentum(double m1, double m2, double m3);
    Vector3D GenerateRandomDirection();

    std::vector<ParticleDef> daughters_;
    std::vector<double> daughter_masses_;
    double sum_daughter_masses_;

    static const std::string name_;
};

} /* PROPOSAL */

