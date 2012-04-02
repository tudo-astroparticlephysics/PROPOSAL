/*
 * PhysicsModel.cxx
 *
 *  Created on: 21.06.2010
 *      Author: koehne
 */

#include "PhysicsModel.h"

double PhysicsModel::ebig_;
double PhysicsModel::elow_;
double PhysicsModel::nlow_;
using namespace std;

// Some defaultvariables


       PhysicsModel::PhysicsModel(){

			
//		Log10_=std::log(10);	// log(10)
//		sqrt2_=std::sqrt(2);	// sqrt(2)
//		sqrt3_=std::sqrt(3);	// sqrt(3)
//		sqrtE_=std::sqrt(exp(1));	// sqrt(e)
                //Alpha_=1/137.03599976;	// fine structure constant

		// computerPrecision_=MathModel.computerPrecision; // Not needed cuz of inheriance
                //HALF_PRECISION=std::sqrt(COMPUTER_PRECISION);
		
//                Cmon_=1/(2*ALPHA);		// monopole charge (in units of e)
//
                elow_=0;
                nlow_=ME;				// bounds of parameterizations
                ebig_=BIGENERGY;			// bounds of parameterizations
                //Rm_ = (1.602176487*pow(10,-19))/(4*PI*8.854187817*pow(10,-12)*Mmu_)*pow(10,-4);



		};

	PhysicsModel::~PhysicsModel(){}


