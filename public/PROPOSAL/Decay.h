/*
 * Decay.h
 *
 *  Created on: 07.05.2013
 *      Author: koehne
 */

#ifndef DECAY_H_
#define DECAY_H_

#include "PROPOSAL/RootFinder.h"
#include "PROPOSAL/Particle.h"


/*! \class  Decay  Decay.h " Decay.h"
   \brief class contains functions necessary for the calculation of decay
 */


class Decay
{

protected:

    RootFinder* root_finder_;
    Particle*   particle_;
    std::string out_;

    bool        store_neutrinos_;
    double      multiplier_;

//----------------------------------------------------------------------------//

    /**
    * function for electron energy calculation - interface to FindRoot
    */

    double Function(double x);

//----------------------------------------------------------------------------//

    /**
    * function for electron energy calculation - interface to FindRoot
    */

    double DifferentiatedFunction(double x);

//----------------------------------------------------------------------------//
public:


    //Constructor
    /**
    * creates internal references to p and m
    */

    Decay();
    Decay(Particle* particle);
    Decay(const Decay&);
    Decay& operator=(const Decay& decay);
    bool operator==(const Decay &decay) const;
    bool operator!=(const Decay &decay) const;
//----------------------------------------------------------------------------//
    //Memberfunctions

    void swap(Decay &decay);
//----------------------------------------------------------------------------//

    /**
    * this cross section describes decay
    */

    double MakeDecay();

//----------------------------------------------------------------------------//

    /**
    * energy of the electron that results from the muon decay
    */

    double CalculateProductEnergy( double ernd, double arnd, double srnd );


//----------------------------------------------------------------------------//

    // Getter
    std::string GetOut() const
    {
        return out_;
    }

    RootFinder* GetRootFinder() const
    {
        return root_finder_;
    }

    Particle* GetParticle() const
    {
        return particle_;
    }

    bool GetStoreNeutrinos() const
    {
        return store_neutrinos_;
    }

    double GetMultiplier() const
    {
        return multiplier_;
    }
//----------------------------------------------------------------------------//
    //Setter
    void SetOut(std::string out);

    void SetRootFinder(RootFinder *root_finder);

    void SetParticle(Particle* particle);

    void SetStoreNeutrinos(bool store_neutrinos);

    void SetMultiplier(double multiplier);
};


#endif /* DECAY_H_ */
