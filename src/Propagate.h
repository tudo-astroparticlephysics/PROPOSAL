/*
 * Propagate.h
 *
 *  Created on: 05.07.2010
 *      Author: koehne
 */

#ifndef PROPAGATE_H_
#define PROPAGATE_H_

#include "PhysicsModel.h"
#include <string>
#include"methods.h"
#include "Integral.h"
#include "StandardNormal.h"
#include "PROPOSALParticle.h"

class PROPOSALParticle;
class Medium;
class CrossSections;
class Energy2Loss;
class Output;

class Interpolate;




/**
 * main mmc class - must be constructed before anything else
 */

class Propagate: public PhysicsModel{

	protected:

                PROPOSALParticle *particle_;
                Medium *medium_;
                CrossSections *cros;
                Energy2Loss *E2Loss_;
		Output *o;

		double rho;		// Preset to 1 in constrcutor


                StandardNormal *StandardN;
                std::vector<Integral*> integral_;
		bool pint;
		double decayS, ionizS, bremsS, epairS, photoS, totalS;


		double ec, tc;


		bool up;		// Preset to true in constructor
		double bigLow[2];
		double storeDif[2];

		bool df;		// Preset as false in constructor
                Interpolate* interpolateJ_;
                Interpolate* interpolateJdf_;


	public:

                bool jt;		// Preset as false in constructor
                bool recc;		// Preset as false in constructor
                bool exactTime; 	// Preset as false in constructor
                bool contiCorr;		// Preset as false in constructor
                bool molieScat;		// Preset as false in constructor
                bool sdec;		// Preset as false in constructor
                static int g;
                bool dw; 		// Preset as false in constructor
                double rw, hw;		// Preset to 0 in constructor
		//----------------------------------------------------------------------------------------------------//

		// Default Constructor

		Propagate();

		
		//----------------------------------------------------------------------------------------------------//

		/**
		 * initialize all classes necessary for propagation of a muon.
		 */




		
		Propagate(std::string w, double ecut, double vcut);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * initialize all classes necessary for propagation of a muon or tau.
		 */

		Propagate(std::string w, double ecut, double vcut, std::string type);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * initialize all classes necessary for propagation. std::string w contains the name of the medium. For the
		 * definition of ecut and vcut see help for the class Medium constructor. Set type to "mu" or "tau".
		 * The value of rho sets the multiplicative medium density correction factor.
		 * To enable parametrization routines, call this.interpolate("all") after this class is created.
		 */

		Propagate(std::string w, double ecut, double vcut, std::string type, double rho);

                //----------------------------------------------------------------------------------------------------//

                /**
                        }

                        if(DEBUG){            * class initializer
                 */

                void init(std::string w, double ecut, double vcut, std::string type, double rho);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * Propagates the particle of initial energy e to the distance r. Returns the final energy if the
		 * particle has survived or the track length to the point of disappearance with a minus sign otherwise.
		 */

		double propagateTo(double r, double e);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * Propagates the particle of initial energy e to the distance r. Returns the final energy if the
		 * particle has survived or the track length to the point of disappearance with a minus sign otherwise.
		 * Also calculates particle energy at point rc. Call getPropEc() to get this energy.
		 */

		double propagateTo(double r, double e, double rc);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * returns the particle energy at rc if the particle has survived or the
		 * distance from the point of decay to rc with a minus sign otherwise.
                 */

		double getPropEc();

		//----------------------------------------------------------------------------------------------------//

		/**
		 * returns the particle time at rc if the particle has survived or at the point of decay otherwise.
		 */

		double getPropTc();

		//----------------------------------------------------------------------------------------------------//

		/**
		 * function for energy range calculation - interface to Integral
		 */

		double function(double E);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * call this routine to enable interpolations. To enable everything, set w="all"
		 */

		void interpolate(std::string w);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * call this routine to enable interpolations and their save and reread. To enable everything, set w="all"
		 */

		void interpolate(std::string w, std::string filename);

		//----------------------------------------------------------------------------------------------------//


		/**
		 * returns the value of the tracking integral
		 */

		double getpr(double ei, double rnd, bool pint);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * final energy, corresponding to the given value rnd of the tracking integral
		 */

		double getef(double ei, double rnd, bool pint);

		//----------------------------------------------------------------------------------------------------//

		/**
		 * 1d parametrization - interface to Interpolate
		 */

		double functionInt(double e);

                // Getter

                bool get_exactTime()  { return exactTime; }
                bool get_molieScat()  { return molieScat;}
                double get_rho()        {return rho;}
                StandardNormal* get_Standard();
                CrossSections* get_cros();
                Medium* get_Medium();
                PROPOSALParticle* get_particle(){return particle_;}
                Output* get_output(){return o;}

                //Setter

                void set_exactTime(bool eTime);
                void set_molieScat(bool mScat);
                void set_contiCorr(bool cCorr);
                void set_pint(bool newPint);
                void set_rho(double newRho);
};


#endif /* PROPAGATE_H_ */
