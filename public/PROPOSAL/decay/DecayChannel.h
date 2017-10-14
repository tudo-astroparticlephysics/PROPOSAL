
#pragma once

#include <string>
#include <vector>

namespace PROPOSAL
{

class Vector3D;
class DynamicData;
class Particle;

class DecayChannel
{

    public:
    typedef std::vector<Particle*> DecayProducts;

    DecayChannel() {}
    virtual ~DecayChannel() {}

    bool operator==(const DecayChannel&) const;
    bool operator!=(const DecayChannel&) const;

    friend std::ostream& operator<<(std::ostream&, DecayChannel const&);

    virtual DecayChannel* clone() = 0;

    virtual DecayProducts Decay(Particle*) = 0;
    static void Boost(Particle*, const Vector3D& direction, double beta);

    virtual const std::string& GetName() const = 0;

    protected:
    virtual bool compare(const DecayChannel&) const = 0;
    virtual void print(std::ostream&) const {};
    // DecayChannel(const DecayChannel&); // Undefined
    // DecayChannel& operator=(const DecayChannel&); // Undefined

    // DecayProducts decay_products;
};

} /* PROPOSAL */
