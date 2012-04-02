/*
 * Scattering.cxx
 *
 *  Created on: 29.06.2010
 *      Author: koehne
 */

#include "Scattering.h"
#include <cmath>
#include "algorithm"

#include "Propagate.h"
#include "PROPOSALParticle.h"
#include "Integral.h"
#include "Interpolate.h"

#include "Bremsstrahlung.h"

using namespace std;

double Scattering::cutoff=1;


	Scattering::Scattering(PROPOSALParticle *p){
		
		df = false;
		jt = false;
                this->particle_=p;
                propagate_=p->propagate_;
                integral_ = new Integral(IROMB, IMAXS, IPREC2);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * function for angle distribution spread calculation - interface to Integral
     */

	double Scattering::function(double E){
//cout<<"Hallo"<<endl;
		const bool DEBUG=false;
		double aux, aux2;
                aux=propagate_->get_cros()->function(E);
                aux2=RY*particle_->e/(particle_->p*particle_->p);
		aux*=aux2*aux2;

		if(DEBUG){
                        cout <<" $" << o->f(particle_->e) <<" \t "<< o->f(aux);
		}
//cout<<aux<<endl;
		return aux;
    }

    //----------------------------------------------------------------------------------------------------//
	double Scattering::gettho(double dr, double ei, double ef){
//            if(dr == 0)
//                return 0;
//cout << "im in gettho" <<endl;
//cout <<"In Scattering:gettho:"<< dr << "\t"<<ei<< "\t"<<ef<<endl;
                double aux;

                if(jt){
                        if(fabs(ei-ef)>fabs(ei)*HALF_PRECISION){
                                aux=interpolateJ_->interpolate(ei);
                                double aux2=aux-interpolateJ_->interpolate(ef);

                                if(fabs(aux2)>fabs(aux)*HALF_PRECISION){
                                        aux=aux2;
                                }
                                else{
                                        aux=interpolateJdf_->interpolate((ei+ef)/2)*(ef-ei);
                                }
                        }
                        else{
                                aux=interpolateJdf_->interpolate((ei+ef)/2)*(ef-ei);
                        }
                }
                else{
//cout<<"Scattering:gettho:79 "<<ei<<"\t"<<ef<<"\t"<<endl;
                        aux=integral_->integrateWithLog(ei, ef, this);
//cout<<"In scattering::gettho: aux ="<<aux<<endl;
                }

		// setLpm sets xo_ as well
                propagate_->get_cros()->get_bremsstrahlung()->setLpm();
               // cout << aux/(Ry_*particle_->e/(particle_->p*particle_->p)*Ry_*particle_->e/(particle_->p*particle_->p)) <<"\t"<< dr<<endl;
                //hotfix aux = dr * ( Ry_*particle_->e/(particle_->p*particle_->p)*Ry_*particle_->e/(particle_->p*particle_->p )
                //double constant = RY*particle_->e/(particle_->p*particle_->p)*RY*particle_->e/(particle_->p*particle_->p);
//                if(particle_->p == 0)
//                        return 0;
                   // cout << aux << "\t"<<dr*constant <<endl;
//                double test;
//                test = min(aux, cutoff);
                //cout<<dr<<"\t"<<ei<<"\t"<<ef<<"\t"<<aux<<"\t"<<fabs(ei-ef)<<"\t"<<fabs(ei)*HALF_PRECISION<<"\t"<<test<<"\t";
                aux=sqrt(max(aux, 0.0)/propagate_->get_cros()->get_bremsstrahlung()->get_xo())*particle_->c*max(1+0.038*log(dr/propagate_->get_cros()->get_bremsstrahlung()->get_xo()), 0.0);


		return min(aux, cutoff);
    }

    //----------------------------------------------------------------------------------------------------//
	double Scattering::functionInt(double e){
		if(df){
			return function(e);
		}
		else{
                        return integral_->integrateWithLog(e, PhysicsModel::ebig_, this);
		}
    }

// Getter


// Setter

        void Scattering::set_jt(bool newJT){ jt = newJT;}

        void Scattering::set_df(bool newDF){ df = newDF;}
