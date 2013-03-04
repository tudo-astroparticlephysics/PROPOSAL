/*
 * IonizStochastic.cxx
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */

#include "PROPOSAL/IonizStochastic.h"
#include "algorithm"
#include <cmath>
#include "PROPOSAL/Medium.h"

using namespace std;


    IonizStochastic::IonizStochastic(Ionizationloss *cros) : Ionizationloss(*cros){

		
                jt_ = false;

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
                        return integral_->integrateWithSubstitution(vUp,vMax,this,1);
		}
    }

    //----------------------------------------------------------------------------------------------------//


	double IonizStochastic::dNdx(double rnd){

                if(cros->get_ci()<=0){
                        return 0.;
		}
                if(jt_){
			this->rnd=rnd;
			return sum=max(interpolateJo_->interpolate(particle_->e), (double) 0);
		}
		else{
			setEnergy();
                        return integral_->integrateWithSubstitution(vUp,vMax,this,1,rnd);
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
		v=vUp*exp(v*log(vMax/vUp));

                return integral_->integrateWithSubstitution(vUp,v,this,1);
    }

    //----------------------------------------------------------------------------------------------------//


	double IonizStochastic::functionInt(double e){

        return interpolateJ_->interpolate(e, 1.0);
    }

    // Setter

        void IonizStochastic::set_jt(bool newJT) { jt_ = newJT; }

    // GEtter

        Interpolate* IonizStochastic::get_interpolateJ() { return interpolateJ_;}
        Interpolate* IonizStochastic::get_interpolateJo() { return interpolateJo_;}


