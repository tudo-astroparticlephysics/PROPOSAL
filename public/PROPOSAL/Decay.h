
/******************************************************************************
 *																			  *
 * This file is part of the simulation tool PROPOSAL.						  *
 *																			  *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,		  *
 * 				      Chair Experimental Physics 5b							  *
 *																			  *
 * This software may be modified and distributed under the terms of a		  *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE". 									  *
 *																			  *
 * Modifcations to the LGPL License:										  *
 *																			  *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the		  *
 *         following reference:												  *
 *																			  *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001										  *
 *																			  *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *		   GitHub webpage													  *
 *																			  *
 *		   "https://github.com/tudo-astroparticlephysics/PROPOSAL"			  *
 *																			  *
 ******************************************************************************/

#pragma once

#ifndef DECAY_H_
#define DECAY_H_

#include "PROPOSAL/RootFinder.h"
#include "PROPOSAL/PROPOSALParticle.h"

namespace PROPOSAL{

/*! \class  Decay  Decay.h " Decay.h"
   \brief class contains functions necessary for the calculation of decay
 */


class Decay
{

protected:

    RootFinder* root_finder_;
    PROPOSALParticle*   particle_;
    PROPOSALParticle*   backup_particle_;
    // std::string out_;
    ParticleType::Enum out_;

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
    Decay(PROPOSALParticle* particle);
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
    ParticleType::Enum GetOut() const
    {
        return out_;
    }

    RootFinder* GetRootFinder() const
    {
        return root_finder_;
    }

    PROPOSALParticle* GetParticle() const
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
    void SetOut(ParticleType::Enum out);

    void SetRootFinder(RootFinder *root_finder);

    void SetParticle(PROPOSALParticle* particle);

    void SetStoreNeutrinos(bool store_neutrinos);

    void SetMultiplier(double multiplier);

    PROPOSALParticle *GetBackup_particle() const;
    void SetBackup_particle(PROPOSALParticle *backup_particle);
    void RestoreBackup_particle();
};

}

#endif /* DECAY_H_ */
