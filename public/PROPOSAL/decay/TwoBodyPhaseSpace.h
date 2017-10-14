
#pragma once

// #include <vector>

#include "PROPOSAL/decay/DecayChannel.h"

namespace PROPOSAL
{

class PROPOSALParticle;

class TwoBodyPhaseSpace : public DecayChannel
{
    public:
    TwoBodyPhaseSpace(double m1, double m2);
    TwoBodyPhaseSpace(const TwoBodyPhaseSpace& mode);
    virtual ~TwoBodyPhaseSpace();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() { return new TwoBodyPhaseSpace(*this); }

    DecayProducts Decay(PROPOSALParticle*);

    const std::string& GetName() const { return name_; }

    private:
    TwoBodyPhaseSpace& operator=(const TwoBodyPhaseSpace&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;
    void print(std::ostream&) const;

    double first_daughter_mass_;
    double second_daughter_mass_;

    static const std::string name_;
};

} /* PROPOSAL */

