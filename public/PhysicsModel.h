/*
 * PhysicsModel.h
 *
 *  Created on: 21.06.2010
 *      Author: koehne
 */

#ifndef PhysicsModel_H
#define PhysicsModel_H

#include <cmath>
#include "FunctionInt.h"
#include "FunctionOfx.h"
#include "FunctionInt2.h"
#include "MathModel.h"



//#ifndef
//    #define
//#endif

/*! \class PhysicsModel PhysicsModel.h "PhysicsModel.h"
 *  \brief This is an entry class of the physical model we are building. Contains constant definitions.
 *
 *
 */

class PhysicsModel: public FunctionInt, public FunctionOfx, public FunctionInt2, public MathModel {


	protected:

//                const static double Pi_=3.141592653589793;	// 3.14159...
//		double Log10_;	// log(10)
//		double sqrt2_;	// sqrt(2)
//		double sqrt3_;	// sqrt(3)
//		double sqrtE_;	// sqrt(e)

//		const static int iromb_=5;			// romb # for integration
//		const static int imaxs_=40;			// max number of int. steps
//		const static double iprec_=1.e-6;		// integration precision
//		const static double iprec2_=1.e-6*10;		// integration precision

//                const static int num1_=100;			// number of interpolation
//		const static int num2_=200;			// points specified mainly
//		const static int num3_=1000;			// in Propagate.java

//		double computerPrecision_;
//		double halfPrecision_;

//                double Alpha_;                              // fine structure constant
//                const static double Me_=0.510998902;		// electron mass (MeV)
//                const static double Ry_=13.60569172;		// Rydberg energy (eV)
//                const static double K_=0.307075;		// in the ionization formula (MeV*cm2/g)

//		const static double C_=2.99792458e10;		// speed of light (cm/s)
//		const static double Re_=2.817940285e-13;	// classical electron radius (cm)
//                double Rm_;                              // classical Muon radius(cm)
//		const static double Na_=6.02214199e23;		// Avogadro's number (1/mol)

//		const static double Mmu_=105.658389;		// muon mass (MeV)
//		const static double Lmu_=2.19703e-6;		// muon lifetime (sec)
//		const static double Mtau_=1777.03;		// tau mass (MeV)
//		const static double Ltau_=290.6e-15;		// tau lifetime (sec)

//		const static double Mpi_=139.57018;		// pion mass (MeV)
//                const static double Mp_=938.271998;		// proton mass (MeV)
//                const static double Mn_=939.56533;		// neutron mass (MeV)
//
//		const static double Mrh_=769.3;			// rho-770 mass (MeV)
//		const static double Ma1_=1230;			// a1-1260 mass (MeV)
//		const static double Mrs_=1465;			// rho-1450 mass (MeV)

//		const static double Gf_=1.16639e-11;		// Fermi coupling const staticant (MeV^-2)
//		const static double Mw_=80419.;			// W+- boson mass (MeV)
//		const static double Mz_=91188.2;		// Z0 boson mass (MeV)
//		const static double Xw_=0.23117;		// sin^2(mixing angle at Mz)
//		const static double Gw_=2120.;			// W+- boson width (MeV)
//		const static double Gz_=2495.2;			// Z0 boson width (MeV)

//		const static double Ds2_=6.9e-5;		// Sun neutrino mass difference (eV^2)
//		const static double Tt2_=0.43;			// tan^2(th) of the mixing angle
//		const static double De2_=2.6e-3;		// Earth neutrino mass difference (eV^2)
//		const static double St2_=1.0;			// sin^2(2 th) of the mixing angle


//                const static double Mmon_=1.e5;			// monopole mass (MeV)
//                double Cmon_;					// monopole charge (in units of e)
//                const static double Mstau_=1.e5;		// stau mass (MeV)
//                const static double Lstau_=-1;			// stau lifetime (sec)
//
//                const static double xres_=1.e-3;		// resolution of particle position (cm)
//                const static double bigEnergy_=1.e14;		// used for radiation length (MeV)
                  static double elow_;			// bounds of parameterizations
                  static double nlow_;				// bounds of parameterizations
                  static double ebig_;			// bounds of parameterizations




	public:


		// constructors

		/// @brief  Create a PhysicsModel.
		PhysicsModel();
		// Some more values which cant be made static 



		// Memberfunctions

		/// @brief  for integration - interface to Integral - to be overwritten by subclasses
		/// @return 0
		double function(double e){return 0;}



		/// @brief  1d parametrization - interface to Interpolate - to be overwritten by subclasses
		/// @return 0
		double functionInt(double e){return 0;}

		/// @brief  2d parametrization - interface to Interpolate - to be overwritten by subclasses
		/// @return 0
		double functionInt(double e, double v){return 0;}

			// Getter

//			double get_Pi() const {return Pi_;}
//			double get_Log10() const {return Log10_;}
//			double get_sqrt2() const {return sqrt2_;}
//			double get_sqrt3() const {return sqrt3_;}
//			double get_sqrtE() const {return sqrtE_;}

//			int get_iromb() const {return iromb_;}
//			int get_imaxs() const {return imaxs_;}
//			double get_iprec() const {return iprec_;}
//			double get_iprec2() const {return iprec2_;}

//			int get_num1() const {return num1_;}
//			int get_num2() const {return num2_;}
//			int get_num3() const {return num3_;}

//			double get_computerPrecision() const {return computerPrecision_;}
//			double get_halfPrecision() const {return halfPrecision_;}

//			double get_Alpha() const {return Alpha_;}
//			double get_Me() const {return Me_;}
//			double get_Ry() const {return Ry_;}
//			double get_K() const {return K_;}
//
//                        double get_C() const {return C_;}
//                        double get_Re() const {return Re_;}
//                        double get_Rm() const {return Rm_;}
//                        double get_Na() const {return Na_;}

//			double get_Mmu() const {return Mmu_;}
//			double get_Lmu() const {return Lmu_;}
//			double get_Mtau() const {return Mtau_;}
//			double get_Ltau() const {return Ltau_;}

//			double get_Mpi() const {return Mpi_;}
//			double get_Mp() const {return Mp_;}
//			double get_Mn() const {return Mn_;}
//
//			double get_Mrh() const {return Mrh_;}
//			double get_Ma1() const {return Ma1_;}
//			double get_Mrs() const {return Mrs_;}

//			double get_Gf() const {return Gf_;}
//			double get_Mw() const {return Mw_;}
//			double get_Mz() const {return Mz_;}
//			double get_Xw() const {return Xw_;}
//			double get_Gw() const {return Gw_;}
//			double get_Gz() const {return Gz_;}

//			double get_Ds2() const {return Ds2_;}
//			double get_Tt2() const {return Tt2_;}
//			double get_De2() const {return De2_;}
//			double get_St2() const {return St2_;}
//
//
//                        double get_Mmon() const {return Mmon_;}
//                        double get_Cmon() const {return Cmon_;}
//                        double get_Mstau() const {return Mstau_;}
//                        double get_Lstau() const {return Lstau_;}
//
//                        double get_xres() const {return xres_;}
//                        double get_bigEnergy() const {return bigEnergy_;}
//                        double get_elow() const {return elow_;}
//                        double get_nlow() const {return nlow_;}
                        static double get_ebig() {return ebig_;}
                        static double get_elow() {return elow_;}

			// Setter
                        static void set_ebig( double newebig ) { ebig_ = newebig; }
                        static void set_elow( double newelow ) { elow_ = newelow; }
			// to be done if necessary

		// destructors

		///@brief Crush this PhysicsModel.
                virtual ~PhysicsModel();


};



#endif //PhysicsModel_H
