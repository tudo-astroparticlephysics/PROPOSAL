
#pragma once

#include <vector>

namespace PROPOSAL
{

class PROPOSALParticle;

class DecayChannel
{

    public:
    typedef std::vector<PROPOSALParticle*> DecayProducts;

    DecayChannel() {}
    virtual ~DecayChannel() {}

    bool operator==(const DecayChannel&) const;
    bool operator!=(const DecayChannel&) const;

    virtual DecayChannel* clone() = 0;

    virtual DecayProducts Decay(PROPOSALParticle*) = 0;

    protected:
    virtual bool compare(const DecayChannel&) const = 0;
    // DecayChannel(const DecayChannel&); // Undefined
    // DecayChannel& operator=(const DecayChannel&); // Undefined

    // DecayProducts decay_products;
};

} /* PROPOSAL */
