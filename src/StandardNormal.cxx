/*
 * StandardNormal.cxx
 *
 *  Created on: 02.08.2010
 *      Author: koehne
 */


/**
 * This class provides routines for evaluation of random numbers distributed normally.
 * @author Dmitry Chirkin
 */

#include "StandardNormal.h"
#include <algorithm>
#include "Interpolate.h"


using namespace std;


    /**
     * initializes class with default settings
     */

StandardNormal::StandardNormal(){
		jt = false;
                norm=1/sqrt(2*PI);
                init(5,20,1.e-6 );

    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes the class
     */

	StandardNormal::StandardNormal(int romberg, int maxSteps, double precision){
		jt = false;
                norm=1/sqrt(2*PI);
		
                init(romberg, maxSteps, precision);

    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * class initializer
     */

        void StandardNormal::init(int romberg, int maxSteps, double precision){

            integral_ = new Integral(romberg, maxSteps, precision);
            val1=sndpr(-1);
            val2=sndpr(1);

        }

    //----------------------------------------------------------------------------------------------------//

    /**
     * evaluates the integrated probability
     */

    double StandardNormal::sndpr(double x){
		if(x<-5){
			return 0;
		}
		else if(x>5){
			return 1;
		}
		if(jt){
                        return min(max(interpolateJ_->interpolate(x), 0.0), 1.0);
		}
		if(x<=-1){
                        return integral_->integrateWithSubstitution(1 , x, this, -2);
		}
		else if(x<=1){
                        return val1+integral_->integrateOpened(-1 , x, this);
		}
		else{
                        return val2+integral_->integrateWithSubstitution(1 , x, this, 2);
		}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * evaluates the standard normal random number
     */

    double StandardNormal::sndrn(double x){
		if(jt){
                        return interpolateJ_->findLimit(x);
		}
                cout<<"StandardNormal::sndrn "<<x<<"\t"<<val1<<"\t"<<val2<<endl;
		if(x<=val1){
                        integral_->integrateWithSubstitution(1 , -1, this, -2, -x);
                        return integral_->getUpperLimit();
		}
		else if(x<=val2){
                        integral_->integrateOpened(-1 , 1, this, -x+val1);
                        return integral_->getUpperLimit();
		}
		else{
                        integral_->integrateWithSubstitution(1 , -1, this, 2, -x+val2);
                        return integral_->getUpperLimit();
		}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * evaluates the standard normal random number
     */

    double StandardNormal::sndrn(double rnd, double average, double sigma, double xmin, double xmax, bool cutoff){

    	double x, xl, xh;

		if(xmax<xmin){
			x=xmin;
			xmin=xmax;
			xmax=x;
		}
		if(sigma==0){
			x=average;
		}
		else{
			if(cutoff){
				x=rnd;
			}
			else{
				xl=sndpr((xmin-average)/sigma);
				xh=sndpr((xmax-average)/sigma);
				x=xl+(xh-xl)*rnd;
			}

			x=average+sigma*sndrn(x);
		}
		if(x<xmin){
			x=xmin;
		}
		if(x>xmax){
			x=xmax;
		}

		return x;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * function describes standard normal distribution - interface to Integral
     */

    double StandardNormal::function(double x){
		double aux=norm*exp(-x*x/2);
                //cout<<"In standardNormal function:\t"<<norm<<"\t"<<aux<<endl;
		return aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * 1d parametrization - interface to Interpolate
     */

    double StandardNormal::functionInt(double x){
    	return sndpr(x);
    }

    // Setter

    void StandardNormal::set_jt(bool newJT){jt=newJT;}

    // Getter

    Interpolate StandardNormal::get_interpolateJ(){return *interpolateJ_;}


