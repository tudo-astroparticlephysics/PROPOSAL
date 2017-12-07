
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

// #include <utility>

#include "PROPOSAL/ProcessCollection.h"
#include "PROPOSAL/ScatteringFirstOrder.h"
#include "PROPOSAL/ScatteringMoliere.h"
// #include "PROPOSAL/MathModel.h"
// #include "PROPOSAL/PROPOSALParticle.h"
// #include "PROPOSAL/Geometry.h"
// #include "PROPOSAL/Scattering.h"

namespace PROPOSAL{

class Propagator :public MathModel
{
private:
    int  order_of_interpolation_;
    bool debug_;
    bool particle_interaction_;     //!< particle interaction? (false = decay)

    int     seed_;                      //!< seed of the random number generator
    ParametrizationType::Enum  brems_;                     //!< Bremsstrahlungs parametrization
    ParametrizationType::Enum  photo_;                     //!< Photonuclear parametrization
    bool    lpm_;                       //!< Landau-Pomeranchuk-Migdal supression of EM cross-sections enabled if true
    bool    moliere_;                   //!< Moliere scattering enabled if true
    bool    stopping_decay_;            //!< Do decay of particles. formarly sdec
    bool    do_exact_time_calculation_; //!< exact local time calculation enabled if true
    bool    integrate_;                 //!< if true nothing will be interpolated
    double  brems_multiplier_;          //!< multiplier to in- or decrease the Bremsstrahlung cross-sections
    double  photo_multiplier_;          //!< multiplier to in- or decrease the Photonucler cross-sections
    double  ioniz_multiplier_;          //!< multiplier to in- or decrease the Ionization cross-sections
    double  epair_multiplier_;          //!< multiplier to in- or decrease the Epairproduction cross-sections
    double  global_ecut_inside_;        //!< ecut for inside the detector (it's used when not specified explicit for a sector in congiguration file)
    double  global_ecut_infront_;       //!< ecut for infront of the detector (it's used when not specified explicit for a sector in congiguration file)
    double  global_ecut_behind_;        //!< ecut for behind the detector (it's used when not specified explicit for a sector in congiguration file)
    double  global_vcut_inside_;        //!< vcut for inside the detector (it's used when not specified explicit for a sector in congiguration file)
    double  global_vcut_infront_;       //!< ecut for infront of the detector (it's used when not specified explicit for a sector in congiguration file)
    double  global_vcut_behind_;        //!< ecut for behind the detector (it's used when not specified explicit for a sector in congiguration file)
    double  global_cont_inside_;        //!< continuous randominzation flag for inside the detector (it's used when not specified explicit for a sector in congiguration file)
    double  global_cont_infront_;       //!< continuous randominzation flag for infront of the detector (it's used when not specified explicit for a sector in congiguration file)
    double  global_cont_behind_;        //!< continuous randominzation flag for behind the detector (it's used when not specified explicit for a sector in congiguration file)

    std::string path_to_tables_;        //!< path to interpolation tables (if not empty tables are stored)
    bool    raw_;                       //!< if true interpolation tables will  be written binary if path is not empty

    std::vector<ProcessCollection*> collections_;

    PROPOSALParticle* particle_;
    //TODO(mario): decide to hold backup particle, because particle could be deleted and cause crashes.
    // So this particle could be reinitialized at the end of propagate(). Di 2017/04/04
    PROPOSALParticle* backup_particle_;
    //FirstOrderScattering
    ScatteringFirstOrder* scatteringFirstOrder_;
    //FirstOrderMoliere
    ScatteringMoliere* scatteringFirstOrderMoliere_;
    int scattering_model_;
    ProcessCollection *current_collection_;

    Geometry*    detector_;


//----------------------------------------------------------------------------//

    /*!
    * Initalize a geomtry. Used when reading the values from config file
    *,
    */
    void InitGeometry(Geometry* geometry, std::deque<std::string>* token , std::string first_token);
//----------------------------------------------------------------------------//
    /*!
    * Init ProcessCollection from configuration file. When keyword sector is found in configuration
    * file this function is called and inits ProcessCollections:
    * 3 for muons inside/behind/infront
    * 3 for taus inside/behind/infront
    * 3 for electrons inside/behind/infront
    */
    void InitProcessCollections(std::ifstream &file);


    void MoveParticle(double distance);

public:

    //Constructors
    Propagator();
    Propagator(ParticleType::Enum particle_type,
               std::string path_to_tables,
               bool exact_time = true,
               bool lpm = true,
               bool integrate = false,
               int scattering_model = 0);
    Propagator(std::string config_file, bool DoApplyOptions=true);
    Propagator(std::string config_file, PROPOSALParticle* particle, bool DoApplyOptions=true);
    Propagator(Medium* medium,
               EnergyCutSettings* cuts,
               ParticleType::Enum particle_type,
               std::string path_to_tables,
               bool moliere = true,
               bool continuous_rand = true,
               bool exact_time = true,
               bool lpm = true,
               ParametrizationType::Enum brems = ParametrizationType::BremsKelnerKokoulinPetrukhin,
               ParametrizationType::Enum photo = ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich,
               double brems_multiplier = 1,
               double photo_multiplier = 1,
               double ioniz_multiplier = 1,
               double epair_multiplier = 1,
               bool integrate = false,
               int scattering_model = 0);
    Propagator(const Propagator&);
    Propagator& operator=(const Propagator& propagator);
    bool operator==(const Propagator &propagator) const;
    bool operator!=(const Propagator &propagator) const;

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
     *  \return energy at distance OR -(track length)
     */

    double Propagate(double distance);

//----------------------------------------------------------------------------//
    /**
     * Propagates the particle through the current set of ProcessCollections
     *  \return vector of secondarys
     */

    std::vector<PROPOSALParticle*> Propagate( PROPOSALParticle *particle, double MaxDistance_cm = 1e20 );


//----------------------------------------------------------------------------//
    std::vector<PROPOSALParticle*> propagate(double MaxDistance_cm = 1e20 ); //TODO(mario): Find new name Fr 2017/03/10

//----------------------------------------------------------------------------//
        /*!
         *  Apply options which are read from configuration file to ProcessCollections
        */
        void ApplyOptions();

//----------------------------------------------------------------------------//
    /*!
    * advances the particle by the given distance
    * Sets the x,y and z coordinates of particle_
    * and its time and propagation distance
    *
    * \param    dr  flight distance
    * \param    ei  initial energy
    * \param    ef  final energy
    */

    void AdvanceParticle(double dr, double ei, double ef);

//----------------------------------------------------------------------------//

    void swap(Propagator &propagator);

//----------------------------------------------------------------------------//

    void InitDefaultCollection(Geometry* geom);

//----------------------------------------------------------------------------//

    void EnableInterpolation(std::string path ="",bool raw=false);

//----------------------------------------------------------------------------//

    void DisableInterpolation();

//----------------------------------------------------------------------------//

    void ReadConfigFile(std::string config_file, bool DoApplyOptions=true);

//----------------------------------------------------------------------------//
    /**
     * Choose the current collection by particle type and location.
     */

    void ChooseCurrentCollection(PROPOSALParticle* particle);

//----------------------------------------------------------------------------//
    /**
     * Calculates the contiuous loss till the first stochastic loss happend
     * and subtract it from initial energy
     * Also caluclate the energy at which the particle decay
     * These to energys can be compared to decide if a decay or particle interaction
     * happens
     *
     *  \param  initial_energy   initial energy
     *  \return pair.first final energy befor first interaction pair.second decay energy at which the
     *          particle decay
     */
    std::pair<double,double> CalculateEnergyTillStochastic( double initial_energy );

//----------------------------------------------------------------------------//
    //Getter

    ProcessCollection* GetCurrentCollection() const
    {
        return current_collection_;
    }
//----------------------------------------------------------------------------//

    std::vector<ProcessCollection*> GetCollections() const
    {
        return collections_;
    }

//----------------------------------------------------------------------------//
    PROPOSALParticle* GetParticle() const
    {
        return particle_;
    }

//----------------------------------------------------------------------------//
    //Setter

    /**
     *  Sets the ProcessCollections. Need to execute AplyOptions() afterward.
     */
    void SetCollections(std::vector<ProcessCollection*>);
    /**
     *  Sets the particle for the Propagator and its current ProcessCollection
     */
    void SetParticle(PROPOSALParticle* particle);

//----------------------------------------------------------------------------//
    //Destructor
    ~Propagator();

    int GetSeed() const;
    void SetSeed(int seed);
    ParametrizationType::Enum GetBrems() const;
    void SetBrems(ParametrizationType::Enum brems);
    ParametrizationType::Enum GetPhoto() const;
    void SetPhoto(ParametrizationType::Enum photo);
    std::string GetPath_to_tables() const;
    void SetPath_to_tables(const std::string &path_to_tables);
    Geometry *GetDetector() const;
    void SetDetector(Geometry *detector);
    bool GetStopping_decay() const;
    void SetStopping_decay(bool stopping_decay);
    void RestoreBackup_particle();
    void ResetParticle();
};

}

#endif // _PROPAGATOR_H
