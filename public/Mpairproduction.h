#ifndef MPAIRPRODUCTION_H
#define MPAIRPRODUCTION_H

/*
 * Epairproduction.h
 *
 *  Created on: 15.03.2011
 *  author: mschmitz
 */


#include "Integral.h"
#include "Interpolate.h"

class CrossSections;

/*! \class Mpairproduction Mpairproduction.h "Mpairproduction.h"
    \brief class contains functions necessary for calculation of \mu+/\mu- pair production losses. All parametrizations are from Kelner, S., Kokoulin, R., & Petrukhin, A., Direct Production of Muon Pairs by High-Energy Muons,Phys. of Atomic Nuclei, Vol. 63, No 9 (2000) 1603. All variable names are consinstent with this paper.

 */


class MpairContinuous;

class MpairStochastic;


class Mpairproduction: public CrossSections{

    protected:

        Integral *integral_;                        ///< Object for integration
        std::vector<Interpolate*> interpolateJ_;    ///< Object for Interpolation

        Mpairproduction *Mpair_;                    ///< Pointer to remember everything
        MpairContinuous *continuous_;               ///< Pointer to remember everything
        MpairStochastic *stochastic_;               ///< Pointer to remember everything

        double v;                                   ///< energy fraction transferred to the particles of a pair
        int NumberOfComponent;                      ///< Which Component of media is used for mpairproduction

        double vMin_;                               ///< integrationlimit for v integration
        double vUp_;                                ///< integrationlimit for v integration
        double vMax_;                               ///< integrationlimit for v integration
        double rhomax;                              ///< integrationlimit for \rho integration
        int form_;                                  ///< parametrizationswitcher
        int cm_;                                    ///< turn mpair off with cm = 0;



        /**
        * Sets the paramter for integration/interpolation
        */
        void setEnergy(int i); 

        /**
        * function which will be integrated.
        */
        double function(double x);

        /**
          * function which will be interpolated
          */
        double functionInt(double e, double v);

        /**
          * This is what U is equal to. See eq. (22) in Kelner, S., Kokoulin, R., & Petrukhin, A., Direct Production of Muon Pairs by High-Energy Muons,Phys. of Atomic Nuclei, Vol. 63, No 9 (2000) 1603 for more information
        */
        double calculateU(double rho);


    public:

        /**
          * usual Constructor
        */
        Mpairproduction(CrossSections *cros);

        Mpairproduction(Mpairproduction *cros_);

        /**
         * this is what \simga dv is equal to
        */

        double mpair(double v, int i);


        void activate_interpolation();

        // Getter

        MpairContinuous* get_Continuous() { return continuous_; }

        int get_cm(){ return cm_;}

        MpairStochastic* get_Stochastic() { return stochastic_; }

        // Setter

        void set_cm(int newcm){cm_ = newcm; }
};

#endif // MPAIRPRODUCTION_H
