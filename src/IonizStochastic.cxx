/*
 * IonizStochastic.cxx
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */

#include "IonizStochastic.h"
#include "algorithm"
#include <cmath>
#include "Medium.h"

using namespace std;


    IonizStochastic::IonizStochastic(Ionizationloss *cros) : Ionizationloss(*cros){

		
                jt_ = false;
		// ?????? super(cros);

                integral_ = new Integral(IROMB, IMAXS, IPREC);
    }

    //----------------------------------------------------------------------------------------------------//


	double IonizStochastic::function(double v){
                return cros->get_ci()*ioniz->d2Ndvdx(v)*(1+ioniz->inelCorrection(v));
    }

    //----------------------------------------------------------------------------------------------------//


	double IonizStochastic::dNdx(){
		if(cros->get_ci()<=0){
			return 0;
		}
                if(jt_){
			return max(interpolateJo_->interpolate(particle_->e), (double) 0);
		}
		else{
			setEnergy();
			return integral_->integrateWithLog(vUp, vMax, this);
		}
    }

    //----------------------------------------------------------------------------------------------------//


	double IonizStochastic::dNdx(double rnd){

                if(cros->get_ci()<=0){
                        return 0.;
		}
                if(jt_){
                    //cout<<"ICH LAUFE"<<endl;
			this->rnd=rnd;
			return sum=max(interpolateJo_->interpolate(particle_->e), (double) 0);
		}
		else{
			setEnergy();
			return integral_->integrateWithLog(vUp, vMax, this, rnd);
		}
    }


    //----------------------------------------------------------------------------------------------------//


	double IonizStochastic::e(double rnd){

		double rand, rsum;
		rand=medium_->get_sumCharge()*rnd;
		rsum=0;
		for(int i=0; i<medium_->get_numCompontents(); i++){
			rsum+=medium_->get_atomInMolecule().at(i)*medium_->get_NucCharge().at(i);

			if(rsum>rand){
				cros->set_component(i);

                                if(jt_){
					setEnergy();

					if(vUp==vMax){
						return particle_->e*vUp;
					}

					return particle_->e*(vUp*exp(interpolateJ_->findLimit(particle_->e, this->rnd*sum)*log(vMax/vUp)));
				}
				else{
					return particle_->e*integral_->getUpperLimit();
				}
			}
		}

		cerr<<"Error (in IonizStochastic/e): m.totZ was not initialized correctly"<<endl;
		return 0;
    }

    //----------------------------------------------------------------------------------------------------//



	double IonizStochastic::functionInt(double e, double v){

		particle_->setEnergy(e);
		setEnergy();

		if(vUp==vMax){
			return 0;
		}
//cout<<"In IonizStochastic::functionInt(double e, double v): \t"<<"vUp = "<<vUp<<"\t"<< "v = "<<v<<endl;
		v=vUp*exp(v*log(vMax/vUp));

                //cout<<integral_->integrateWithLog(vUp, v, this)<<endl;
		return integral_->integrateWithLog(vUp, v, this);
    }

    //----------------------------------------------------------------------------------------------------//


	double IonizStochastic::functionInt(double e){
           // cout<<"Hier passiert was\t"<<e<<"\t"<<interpolateJ_->interpolate(e, 1)<<endl;
           // cout<<interpolateJ_->get_interpolate().at(1)->get_iY().at(1)<<endl;
            //cout<<"Stoch Zeile 130:\t"<<e<<endl;
        return interpolateJ_->interpolate(e, 1.0); //<--- Hier stand 1 statt 1.0 07.12.10 jhk

    }

    // Setter

        void IonizStochastic::set_jt(bool newJT) { jt_ = newJT; }

    // GEtter

        Interpolate* IonizStochastic::get_interpolateJ() { return interpolateJ_;};
        Interpolate* IonizStochastic::get_interpolateJo() { return interpolateJo_;};


