
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
        class Builder;

    public:
    ManyBodyPhaseSpace(std::vector<ParticleDef> daughters);
    ManyBodyPhaseSpace(const ManyBodyPhaseSpace& mode);
    virtual ~ManyBodyPhaseSpace();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() const { return new ManyBodyPhaseSpace(*this); }

    // ----------------------------------------------------------------------------
    /// @brief Many body phase space decay
    ///
    /// Calculate decay products with the help of the Raubold Lynch algorithm.
    ///
    /// @param PROPOSALParticle
    ///
    /// @return Vector of particles, the decay products
    // ----------------------------------------------------------------------------
    DecayProducts Decay(const Particle&);

    const std::string& GetName() const { return name_; }

    private:
    ManyBodyPhaseSpace& operator=(const ManyBodyPhaseSpace&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;
    void print(std::ostream&) const;

    std::vector<ParticleDef> daughters_;
    std::vector<double> daughter_masses_;
    double sum_daughter_masses_;

    static const std::string name_;
};

class ManyBodyPhaseSpace::Builder
{
    public:

    Builder();
    Builder(const Builder&);
    ~Builder();

    // --------------------------------------------------------------------- //
    // Setter
    // --------------------------------------------------------------------- //

    Builder& addDaughter(const ParticleDef& daughter);
    ManyBodyPhaseSpace build();

    private:
    std::vector<ParticleDef> daughters_;

};

} /* PROPOSAL */

