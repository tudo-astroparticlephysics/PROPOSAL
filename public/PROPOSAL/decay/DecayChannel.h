
#pragma once

#include <string>
#include <vector>

namespace PROPOSAL
{

class Vector3D;
class DynamicData;
class PROPOSALParticle;

class DecayChannel
{

    public:
    typedef std::vector<PROPOSALParticle*> DecayProducts;

    DecayChannel() {}
    virtual ~DecayChannel() {}

    bool operator==(const DecayChannel&) const;
    bool operator!=(const DecayChannel&) const;

    friend std::ostream& operator<<(std::ostream&, DecayChannel const&);

    virtual DecayChannel* clone() = 0;

    virtual DecayProducts Decay(PROPOSALParticle*) = 0;
    static void Boost(PROPOSALParticle*, const Vector3D& direction, double beta);

    virtual const std::string& GetName() const = 0;

    protected:
    virtual bool compare(const DecayChannel&) const = 0;
    virtual void print(std::ostream&) const {};
    // DecayChannel(const DecayChannel&); // Undefined
    // DecayChannel& operator=(const DecayChannel&); // Undefined

    // DecayProducts decay_products;
};

} /* PROPOSAL */
