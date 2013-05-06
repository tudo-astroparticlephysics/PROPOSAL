/*
 * Propagator.h
 *
 *  Created on: 2013.05.06
 *      Author: Jan-Hendrik Köhne
 */
#ifndef PROPAGATOR_H
#define PROPAGATOR_H


/*! \mainpage PROPOSAL:
 * <b>PR</b>opagator with <b>O</b>ptimal <b>P</b>recision and <b>O</b>ptimized <b>S</b>peed for <b>A</b>ll <b>L</b>eptons.
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection requirements
 *
 * \section HowTo
 *
 */

#include "PROPOSAL/Particle.h"
#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/ProcessCollection.h"


class Propagator :public MathModel
{
private:
    bool debug_;
    bool particle_interaction_;
    double rho_;                    //!< multiplicative medium density correction factor ????????????????????????????????????
    bool do_weighting_;             //!< Do weigthing? Set as false in constructor
    double weighting_order_;        //!< Re-weighting order. Set to 0 in constructor
    double weighting_starts_at_;    //!< Distance at which re-weighting starts. Set to 0 in constructor
    std::vector<double> rates_;
    double total_rate_;


    Particle* particle_;
    ProcessCollection *collection_;




public:

    //Constructors
    Propagator();
    Propagator(const Propagator&);
    Propagator& operator=(const Propagator& propagator);
    bool operator==(const Propagator &propagator) const;
    bool operator!=(const Propagator &propagator) const;

//----------------------------------------------------------------------------//
    //Destructor
    ~Propagator();
//----------------------------------------------------------------------------//
    //Memberfunctions
    /**
     * Propagates the particle of initial energy e to the distance r.
     * Returns the final energy if the
     * particle has survived or the track length to the
     * point of disappearance with a minus sign otherwise.
     *
     *  \param  distance   maximum track length
     *  \param  energy   initial energy
     *  \return energy at distance r OR -(track length)
     */

    double Propagate(double distance, double energy);
//----------------------------------------------------------------------------//

    void swap(Propagator &propagator);
//----------------------------------------------------------------------------//

};

#endif // _PROPAGATOR_H
