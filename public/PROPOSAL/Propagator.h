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
#include <utility>
#include <boost/program_options.hpp>
#include "PROPOSAL/Geometry.h"


class Propagator :public MathModel
{
private:
    int  order_of_interpolation_;
    bool debug_;
    bool particle_interaction_;     //!< particle interaction? (false = decay)
    bool do_time_interpolation_;    //!< If true, CalculateParticleTime is interpolated


    int     seed_;                      //!< seed of the random number generator
    int     brems_;                     //!< Bremsstrahlungs parametrization
    int     photo_;                     //!< Photonuclear parametrization
    bool    lpm_;                       //!< Landau-Pomeranchuk-Migdal supression of EM cross-sections enabled if true
    bool    moliere_;                   //!< Moliere scattering enabled if true
    bool    do_exact_time_calulation_;  //!< exact local time calculation enabled if true
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


    std::vector<ProcessCollection*> collections_;

    Particle* particle_;
    ProcessCollection *collection_;

    Interpolant* interpol_time_particle_;
    Interpolant* interpol_time_particle_diff_;

    Integral*    time_particle_;

    Geometry*    detector_;


//----------------------------------------------------------------------------//
    /*!
    * function for time delta calculation - interface to Integral
    *
    */

    double FunctionToTimeIntegral(double E);

//----------------------------------------------------------------------------//

    double InterpolTimeParticle(double energy);

//----------------------------------------------------------------------------//
    double InterpolTimeParticleDiff(double energy);

//----------------------------------------------------------------------------//
    /*!
    * Initalize a geomtry. Used when reading the values from config file
    *
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

//----------------------------------------------------------------------------//
    /*!
     *  Apply options which are read from configuration file to ProcessCollections
    */
    void ApplyOptions();

public:

    //Constructors
    Propagator();
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
     *  \return energy at distance r OR -(track length)
     */

    double Propagate(double distance);

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
    /*!
    * time delta, corresponding to the given propagation distance
    *
    * \param    ei  initial energy
    * \param    ef  final energy
    * \return   time delta
    */

    double CalculateParticleTime(double ei, double ef);

//----------------------------------------------------------------------------//
    void swap(Propagator &propagator);

//----------------------------------------------------------------------------//

    void InitDefaultCollection();

//----------------------------------------------------------------------------//

    void EnableInterpolation(std::string path ="");

//----------------------------------------------------------------------------//

    void DisableInterpolation();

//----------------------------------------------------------------------------//

    void ReadConfigFile(std::string config_file);

//----------------------------------------------------------------------------//
    /**
     * Disables the Interpolation for calculation the exact particle time
     */
    void DisableParticleTimeInterpolation();
//----------------------------------------------------------------------------//

    /**
     * Enables the Interpolation for calculation the exact particle time
     */
    void EnableParticleTimeInterpolation(std::string path);
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

    void Setup(int argc, char** argv);

//----------------------------------------------------------------------------//

    boost::program_options::options_description CreateOptions();

//----------------------------------------------------------------------------//
    //Getter

    ProcessCollection* GetCollection() const
    {
        return collection_;
    }
//----------------------------------------------------------------------------//

    std::vector<ProcessCollection*> GetCollections() const
    {
        return collections_;
    }

//----------------------------------------------------------------------------//
    Particle* GetParticle() const
    {
        return particle_;
    }

 //----------------------------------------------------------------------------//
    //Destructor
    ~Propagator();

};

#endif // _PROPAGATOR_H
