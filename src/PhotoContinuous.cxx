/*
 * PhotoContinuous.cxx
 *
 *  Created on: 28.06.2010
 *      Author: koehne
 */


#include "PhotoContinuous.h"

#include "algorithm"
#include "Medium.h"

using namespace std;

	PhotoContinuous::PhotoContinuous(Photonuclear *cros) : Photonuclear(*cros){

                jt_ = false;
		// ??????? super(cros);
                integral_ = new Integral(IROMB, IMAXS, IPREC);
    }

    //----------------------------------------------------------------------------------------------------//

	double PhotoContinuous::function(double v){
		return v*photo->photoN(v, cros->get_component());
    }

    //----------------------------------------------------------------------------------------------------//

 	double PhotoContinuous::dEdx(){
//cout<<"In PhotoContinuous::dEdx(): cros->get_cp() = "<<cros->get_cp()<<endl;

		if(cros->get_cp()<=0){
			return 0;
		}
                if(jt_){
                        return max(interpolateJ_->interpolate(particle_->e), 0.0);
		}

		double sum = 0;

		for(int i=0; i < medium_->get_numCompontents(); i++){
			setEnergy(i);
			sum+=integral_->integrateWithLog(vMin, vUp, this);
//cout<<"In PhotoContinuous::dEdx(): sum = "<<sum<<endl;
		}
		return cros->get_cp()*particle_->e*sum;
    }

    //----------------------------------------------------------------------------------------------------//

	double PhotoContinuous::functionInt(double e){
		particle_->setEnergy(e);
		return dEdx();
                }



