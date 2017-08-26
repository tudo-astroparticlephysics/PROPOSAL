/*
 * Decay.h
 *
 *  Created on: 07.05.2013
 *      Author: koehne
 */
#pragma once

#include <vector>

#include "PROPOSAL/RootFinder.h"

namespace PROPOSAL {

class PROPOSALParticle;

/*! \class  Decay  Decay.h " Decay.h"
   \brief class contains functions necessary for the calculation of decay
 */
class Decay
{

    protected:
    // ParticleType::Enum out_;

    bool store_neutrinos_;
    double multiplier_;

    public:
    Decay();

    Decay(const Decay&);

    Decay& operator=(const Decay& decay);
    bool operator==(const Decay& decay) const;
    bool operator!=(const Decay& decay) const;

    // Memberfunctions
    void swap(Decay& decay);

    /**
    * this cross section describes decay
    */
    double MakeDecay(const PROPOSALParticle&);

    // double CalculateProductEnergy( double ernd, double arnd, double srnd );

    // Getter
    bool GetStoreNeutrinos() const { return store_neutrinos_; }
    double GetMultiplier() const { return multiplier_; }

    // Setter
    void SetStoreNeutrinos(bool store_neutrinos);
    void SetMultiplier(double multiplier);
};

// ----------------------------------------------------------------------------
/// @brief Abstract class of decay channels
// ----------------------------------------------------------------------------
class DecayChannel
{

    public:
    typedef std::vector<PROPOSALParticle*> DecayProducts;

    // DecayChannel();
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

class StableChannel : public DecayChannel
{
    public:
    StableChannel();
    StableChannel(const StableChannel& mode);
    virtual ~StableChannel();
    // No copy and assignemnt -> done by clone
    StableChannel* clone() { return new StableChannel(*this); }

    DecayProducts Decay(PROPOSALParticle*);

    private:
    StableChannel& operator=(const StableChannel&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;

    private:
    /* data */
};

// ----------------------------------------------------------------------------
/// @brief Leptonic decay channel
// ----------------------------------------------------------------------------
class LeptonicDecayChannel : public DecayChannel
{
    public:
    LeptonicDecayChannel();
    LeptonicDecayChannel(const LeptonicDecayChannel& mode);
    virtual ~LeptonicDecayChannel();
    // No copy and assignemnt -> done by clone
    LeptonicDecayChannel* clone() { return new LeptonicDecayChannel(*this); }

    DecayProducts Decay(PROPOSALParticle*);

    private:
    LeptonicDecayChannel& operator=(const LeptonicDecayChannel&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;

    // ----------------------------------------------------------------------------
    /// @brief Function for electron energy calculation - interface to FindRoot
    // ----------------------------------------------------------------------------
    double DecayRate(double);

    // ----------------------------------------------------------------------------
    /// @brief Function for electron energy calculation - interface to FindRoot
    // ----------------------------------------------------------------------------
    double DifferentialDecayRate(double);

    RootFinder root_finder_;
};

// ----------------------------------------------------------------------------
/// @brief Two body phase space
// ----------------------------------------------------------------------------
class TwoBodyPhaseSpace : public DecayChannel
{
    public:
    TwoBodyPhaseSpace(double m1, double m2);
    TwoBodyPhaseSpace(const TwoBodyPhaseSpace& mode);
    virtual ~TwoBodyPhaseSpace();
    // No copy and assignemnt -> done by clone
    virtual TwoBodyPhaseSpace* clone() { return new TwoBodyPhaseSpace(*this); }

    DecayProducts Decay(PROPOSALParticle*);

    private:
    TwoBodyPhaseSpace& operator=(const TwoBodyPhaseSpace&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;

    double first_daughter_mass_;
    double second_daughter_mass_;
};

}
