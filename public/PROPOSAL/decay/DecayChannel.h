
#pragma once

#include <vector>

namespace PROPOSAL
{

class Particle;

class DecayChannel
{

    public:
    typedef std::vector<Particle*> DecayProducts;

    DecayChannel() {}
    virtual ~DecayChannel() {}

    bool operator==(const DecayChannel&) const;
    bool operator!=(const DecayChannel&) const;

    virtual DecayChannel* clone() = 0;

    virtual DecayProducts Decay(Particle*) = 0;

    protected:
    virtual bool compare(const DecayChannel&) const = 0;
    // DecayChannel(const DecayChannel&); // Undefined
    // DecayChannel& operator=(const DecayChannel&); // Undefined

    // DecayProducts decay_products;
};

} /* PROPOSAL */
