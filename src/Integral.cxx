/*
* Integral.cxx
*
*  Created on: 02.08.2010
*      Author: koehne
*/

// package mmc;

/**
* This class provides routines for function integration using Romberg method.
* Include the function to be integrated in a class that implements the interface FunctionOfx (defined below).
* Methods contained here are based on the Numerical Recipes (W. H. Press et al.) with the following modifications:
* <ul><li>Algorithm of route choice in the interpolate method is different: it keeps the partial approximations
* centered on the target x by keeping two of the grid points closest to x in the middle as long as possible.</li>
* <li>Power of Substitution (POS) is any real number. If it is zero, no substitution is made. Otherwise the
* substitution is 1/x^POS if POS&gt;0 or 1/(-x)^(-POS), if POS&lt;0. If one of the limits of integration is zero or
* has opposite sign than that of POS, it is replaced with infinity. If both limits of integration are zero or have
* opposite sign than that of POS, integral is considered to have equal to each other (and to infinity) limits, and
* the returned value is zero.</li></ul>
* It is possible to evaluate x(rand) such that the integral from xmin to x(rand) is a fraction rand of the original
* full integral from xmin to xmax. Set the ratio 0&lt;rand&lt;1 as the last argument of one of the open integration
* methods (integrateOpened or integrateWithSubstitution), and get x(rand) by the subsequent call to getUpperLimit.
* If rand&lt;0, then -rand is assumed to be the absolute value of the portion of the original integral such that the
* integral from xmin to x(rand) is equal to this portion. Its sign is determined as the sign of the whole integral
* (or, rather, the N-1st approximation to its value). If rand is given, it is generally assumed that the function
* does not change sign on the integration interval. Otherwise, the resulting x(rand) is less predictable. In any
* case, an approximation to x(rand) is found during evaluation of the original integral, and then refined by the
* combination of the Newton-Raphson method and bisection.
* <pre>
* interface FunctionOfx{
*     double function(double x);
* }
* </pre>
* @author Dmitry Chirkin
*/
#include "Integral.h"
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;


//----------------------------------------------------------------------------------------------------//

/**
* initializes class with default settings
*/

        Integral::Integral()
        {
            // Integral is used with presets defined below

            maxSteps = 20;
            romberg = 5;
            precision = 1.e-6;

            max=1;
            min=0;
            romberg4refine=2;

            powerOfSubstitution=0;

            randomDo=false;
            useLog=false;
            reverse=false;

            int aux;
            if(romberg<=0)
            {
                printf("Warning (in Integral/Integral/0): romberg = %i must be > 0, setting to 1",romberg);
                romberg=1;
            }
            if(maxSteps<=0)
            {
                printf("Warning (in Integral/Integral/1): maxSteps = %i must be > 0, setting to 1", maxSteps);
                maxSteps=1;
            }
            if(precision<=0)
            {
                printf("Warning (in Integral/Integral/2): precision = %f must be > 0, setting to 1.e-6",precision);
                precision=1.e-6;
            }
            this->romberg=romberg;
            this->maxSteps=maxSteps;
            this->precision=precision;
            iX.resize(maxSteps);
            iY.resize(maxSteps);
            aux=std::max(romberg, romberg4refine);
            c.resize(aux);
            d.resize(aux);

        }

        //----------------------------------------------------------------------------------------------------//

        /**
        * initializes class - this is the main constructor
        */

        Integral::Integral(int  romberg, int maxSteps, double precision)
        {
            max=1;
            min=0;
            romberg4refine=2;

            powerOfSubstitution=0;

            randomDo=false;
            useLog=false;
            reverse=false;



            int aux;
            if(romberg<=0)
            {
                printf("Warning (in Integral/Integral/0): romberg = %i must be > 0, setting to 1",romberg);
                romberg=1;
            }
            if(maxSteps<=0)
            {
                printf("Warning (in Integral/Integral/1): maxSteps = %i must be > 0, setting to 1", maxSteps);
                maxSteps=1;
            }
            if(precision<=0)
            {
                printf("Warning (in Integral/Integral/2): precision = %f must be > 0, setting to 1.e-6",precision);
                precision=1.e-6;
            }
            this->romberg=romberg;
            this->maxSteps=maxSteps;
            this->precision=precision;
            iX.resize(maxSteps);
            iY.resize(maxSteps);
            aux=std::max(romberg, romberg4refine);
            c.resize(aux);
            d.resize(aux);
        }

        //----------------------------------------------------------------------------------------------------//

        /*
        * table of substitutions
        */

        double Integral::function(double x)
        {
            double result, t;
            if(reverse)
            {
                x=reverseX-x;
            }
            if(powerOfSubstitution==0)
            {
                t=x;
                result=1;
            }
            else if(powerOfSubstitution>0)
            {
                t=pow(x, -powerOfSubstitution);
                result=powerOfSubstitution*(t/x);
            }
            else
            {
                t=-pow(-x, powerOfSubstitution);
                result=-powerOfSubstitution*(t/x);
            }
            if(useLog)
            {
                t=exp(t);
                result*=t;
            }

            result*=function2use->function(t);
            return result;
        }

        //----------------------------------------------------------------------------------------------------//

        /*
        * returns corrected integral by trapezoid rule, n=1, 2, 4, 8, ...
        */

        double Integral::trapezoid(int n, double oldSum)
        {
            double xStep, stepSize, resultSum;

            if(n==1)
            {
                return (function(max)+function(min))*(max-min)/2;
            }

            n/=2;
            stepSize=(max-min)/n;
            resultSum=0;
            for(xStep=min+stepSize/2;xStep<max;xStep+=stepSize)
            {
                resultSum+=function(xStep);
            }
            return (oldSum+resultSum*stepSize)/2;
        }

        //----------------------------------------------------------------------------------------------------//

        /*
        * returns corrected integral by opened trapezoid rule, n=1, 3, 9, ...
        */

        double Integral::trapezoid3(int n, double oldSum){
            double xStep, stepSize, resultSum;

            if(n==1)
            {

                return (max-min)*function((max+min)/2);
            }
            stepSize=(max-min)/n;
            resultSum=0;

            for(xStep=min+stepSize/2;xStep<max;xStep+=stepSize)
            {
                resultSum+=function(xStep);
                xStep+=2*stepSize;
                resultSum+=function(xStep);

            }
//cout<<oldSum/3+resultSum*stepSize<<endl;
            return oldSum/3+resultSum*stepSize;
        }

        //----------------------------------------------------------------------------------------------------//

        /*
        * returns corrected integral by opened trapezoid rule, n=1, 3, 9, ...
        * and computes the approximation to the value of the x(rand)
        */

        double Integral::trapezoid3S(int n, double oldSum, int stepNumber)
        {
            double xStep, stepSize, resultSum;
            double smallSum=0, approX=0, functionValue1, functionValue2, sumDifference;
            double functionDifference, functionSum, aEq, bEq, bEq2, determinant;
            bool flag;
            if(n==1)
            {

                return (max-min)*function((max+min)/2);
            }

            stepSize=(max-min)/n;
            if(stepNumber>=romberg-1)
            {
                if(randomNumber>=0)
                {
                    smallSum=randomNumber*oldSum/(1.5*stepSize);
                }
                else
                {
                    smallSum=-randomNumber/(1.5*stepSize);
                    if(oldSum<0) smallSum*=-1;
                }
            }
            resultSum=0;
            flag=false;
            for(xStep=min+stepSize/2;xStep<max;xStep+=stepSize)
            {
                resultSum+=functionValue1=function(xStep);
                xStep+=2*stepSize;
                resultSum+=functionValue2=function(xStep);
                if(!flag) if(stepNumber>=romberg-1) if((resultSum>=smallSum && smallSum>0) || (resultSum<=smallSum && smallSum<0)){
                        functionSum=functionValue1+functionValue2;
                        sumDifference=(smallSum-(resultSum-functionSum))*1.5*stepSize;
                        functionDifference=functionValue2-functionValue1;
                        aEq=functionDifference/stepSize;
                        bEq=(functionValue2-5*functionValue1)/2;
                        bEq2=bEq*bEq;
                        if(fabs(aEq*sumDifference)<precision*bEq2) approX=sumDifference*2/functionSum;
                        else{
                                determinant=bEq2+4*aEq*sumDifference;
                                if(determinant>=0){
                                        determinant=sqrt(determinant);
                                        approX=(bEq+determinant)/aEq;
                                        if(approX<0 || approX>3*stepSize) approX=(bEq-determinant)/aEq;
                                        else if(approX<0 || approX>3*stepSize) approX=sumDifference*2/functionSum;
                                }
                                else{
                                        printf("Warning (in Integral/trapezoid3S): Determinant is negative, proceeding with linear approximation");
                                        approX=sumDifference*2/functionSum;
                                }
                        }
                        approX+=xStep-2.5*stepSize;
                        flag=true;
                }
            }
            if(stepNumber>=romberg-1){
                    if(!flag) approX=max;
                    randomX=approX;
            }

            return oldSum/3+resultSum*stepSize;
        }

//----------------------------------------------------------------------------------------------------//

/*
* finds f(x) for f: iY[i]=f(iX[i]), i=start, ..., start+romberg
*/

void Integral::interpolate(int start, double x){
	int num, i, k;
	bool dd;
	double error=0, result=0;
	double aux, aux2, dx1, dx2;
	num=0; aux=fabs(x-iX[start+0]);

	for(i=0;i<romberg;i++){
		aux2=fabs(x-iX[start+i]);
		if(aux2<aux){ num=i; aux=aux2; }
		c[i]=d[i]=iY[start+i];

	}
	if(num==0) dd=true;
	else if(num==romberg-1) dd=false;
	else{
		if(fabs(x-iX[start+num-1]) > fabs(x-iX[start+num+1])) dd=true;
		else dd=false;
	}
        //cout<<"iY: "<<iY[start+num]<<endl;
	result=iY[start+num];
	for(k=1;k<romberg;k++){
		for(i=0;i<romberg-k;i++){
			dx1=iX[start+i]-x;
			dx2=iX[start+i+k]-x;

			aux=c[i+1]-d[i];

			aux2=dx1-dx2;
			if(aux2!=0){

				aux=aux/aux2;
				c[i]=dx1*aux;
				d[i]=dx2*aux;


			}
			else{
				c[i]=0;
				d[i]=0;
			}
		}
		if(num==0) dd=true;
		if(num==romberg-k) dd=false;

		if(dd) error=c[num];
		else{ num--; error=d[num]; }
		dd=!dd;
                //cout <<"result: "<<result<<"\t"<<"error: "<<error<<endl;
		result+=error;
	}

	integralError=error;
	integralValue=result;
       // cout<<"integralvalue in interpolate\t"<<integralValue<<endl;
}

//----------------------------------------------------------------------------------------------------//

/*
* finds integral for closed intervals
*/

double Integral::rombergIntegrateClosed(){
	double k = 1;
	double n = 1;
	double error, result, value=0;
	result=0;

	for(int i=0;i<maxSteps;i++){
		result=trapezoid(k, result);
		iX[i]=n;
		iY[i]=result;
		if(i>=romberg-1){
			interpolate(i-(romberg-1), 0);
			error=integralError;
			value=integralValue;
			if(value!=0) error/=value;
			if(fabs(error)<precision) return value;
		}	

		k=k*2; 
		n=n/4;
	}
        printf("Warning (in Integral/rombergIntegrateClosed): Precision %f has not been reached after %i steps \n", precision, maxSteps);
	return value;
}

//----------------------------------------------------------------------------------------------------//

/*
* finds integral for opened intervals
*/

double Integral::rombergIntegrateOpened(){
                   // cout<<"rombergIntegrateOpened läuft an"<<endl;
	int i, k;
	double n;
	double error, result, value=0;
	k=1; n=1;
	result=0;
//        if(min<3)
//   cout<<"\t \t Maxsteps: "<<maxSteps<<endl;
	for(i=0;i<maxSteps;i++){

		if(randomNumber==0 || randomNumber==1){

			result=trapezoid3(k, result);
                     //   cout << "\t \t result="<<result<<endl;
                        if(randomNumber==0){
                            randomX=min;
                        }
                        else{
                            randomX=max;
                        }
		}
                else{
                    result=trapezoid3S(k, result, i);
                }
//cout<<"Integral::rombergIntegralOpened():430 "<<i<<endl;
		iX[i]=n;
		iY[i]=result;
                //cout<<iY.size()<<"\t"<<result<<endl;
		if(i>=romberg-1){
			interpolate(i-(romberg-1), 0);
			error=integralError;
			value=integralValue;
                        //cout<<"ERROR: \t"<<error<<endl;
			if(value!=0) error/=value;
                            //cout << "FABS UND SO\t "<<fabs(error)<<"\t"<<precision<<endl;
                           // cout << "VALUE: "<<value<<endl;
                        if(fabs(error)<precision)

                            return value;
               }
		k*=3; n/=9;
	}
        printf("Warning (in Integral/rombergIntegrateOpened/0): Precision %f has not been reached after %i steps \n", precision, maxSteps);
	return value;
}

//----------------------------------------------------------------------------------------------------//

/*
* finds integral for opened intervals; precision of the result
* is evaluated with respect to the value provided in the argument
*/

double Integral::rombergIntegrateOpened(double bigValue){
	int i, k;
	double n;
	double error, result, value=0;
	k=1; n=1;
	result=0;
	for(i=0;i<maxSteps;i++){
		result=trapezoid3(k, result);
		iX[i]=n;
		iY[i]=result;
		if(i>=romberg-1){
			interpolate(i-(romberg-1), 0);
			error=integralError;
			value=integralValue;
			error/=bigValue;
                        if(fabs(error)<precision)
                            return value;
		}
		k*=3; n/=9;
	}
        printf("Warning (in Integral/rombergIntegrateOpened/1): Precision %f has not been reached after %i steps \n", precision, maxSteps);
	return value;
}

//----------------------------------------------------------------------------------------------------//

/**
* finds integral for closed intervals
*/

double Integral::integrateClosed(double min, double max, FunctionOfx *function2use){
	double aux;
	reverse=false;
	useLog=false;
        if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION) return 0;
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; }
	else aux=1;
	this->min=min;
	this->max=max;
	this->function2use=function2use;
	powerOfSubstitution=0;
	randomDo=false;
	return aux*rombergIntegrateClosed();
}

//----------------------------------------------------------------------------------------------------//

/**
* finds integral for opened intervals
*/

double Integral::integrateOpened(double min, double max, FunctionOfx *function2use){
	double aux;
	reverse=false;
	useLog=false;
        if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION) return 0;
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; }
	else aux=1;
	this->min=min;
	this->max=max;
	this->function2use=function2use;
	powerOfSubstitution=0;
	randomNumber=0;
	randomDo=false;
	return aux*rombergIntegrateOpened();
}

//----------------------------------------------------------------------------------------------------//

/**
* finds integral for opened intervals
* and computes the value of the x(rand)
*/

double Integral::integrateOpened(double min, double max, FunctionOfx *function2use, double randomRatio){
	double aux, result;
	reverse=false;
	useLog=false;
        if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION){ this->min=min; this->max=max; randomDo=true; return 0; }
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; reverse=!reverse; }
	else aux=1;
	this->min=min;
	this->max=max;
	if(reverse) reverseX=this->min+this->max;
	this->function2use=function2use;
	powerOfSubstitution=0;
	if(randomRatio>1) randomRatio=1;
	randomNumber=randomRatio;
	result=rombergIntegrateOpened();
	if(randomNumber<0){
		randomNumber/=-fabs(result);
		if(randomNumber>1) randomNumber=1;
		if(randomNumber<0) randomNumber=0;
	}
	savedResult=result;
	randomDo=true;
	return aux*result;
}

//----------------------------------------------------------------------------------------------------//

/**
* finds integral for opened intervals using substitution x -&gt; 1/x^(powerOfSubstitution)
*/

double Integral::integrateWithSubstitution(double min, double max, FunctionOfx *function2use, double powerOfSubstitution){
	double aux;
	reverse=false;
	useLog=false;
        if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION) return 0;
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; }
	else aux=1;
	if(powerOfSubstitution>0){
		if(max>0 && min>0){
			this->min= pow(max, -1/powerOfSubstitution);
			this->max= pow(min, -1/powerOfSubstitution);
		}
		else if(max>0){
			this->min=0;
			this->max= pow(max, -1/powerOfSubstitution);
			aux=-aux;
		}
		else return 0;
	}
	else if(powerOfSubstitution<0){
		if(max<0 && min<0){
			this->min=- pow(-max, 1/powerOfSubstitution);
			this->max=- pow(-min, 1/powerOfSubstitution);
		}
		else if(min<0){
			this->min=- pow(-min, 1/powerOfSubstitution);
			this->max=0;
			aux=-aux;
		}
		else return 0;
	}
	else{
		this->min=min;
		this->max=max;
	}
	this->function2use = function2use;

	this->powerOfSubstitution=powerOfSubstitution;
	randomNumber=0;
	randomDo=false;
	return aux*rombergIntegrateOpened();
}

//----------------------------------------------------------------------------------------------------//

/**
* finds integral for opened intervals using substitution x -&gt; 1/x^(powerOfSubstitution)
* and computes the value of the x(rand)
*/

double Integral::integrateWithSubstitution(double min, double max, FunctionOfx *function2use, double powerOfSubstitution, double randomRatio){
	double aux, result;
	reverse=false;
	useLog=false;
        if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION){ this->min=min; this->max=max; randomDo=true; return 0; }
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; reverse=!reverse; }
	else aux=1;
	if(powerOfSubstitution>0){
		if(max>0 && min>0){
			this->min= pow(max, -1/powerOfSubstitution);
			this->max= pow(min, -1/powerOfSubstitution);
			reverse=!reverse;
		}
		else if(max>0){
			this->min=0;
			this->max= pow(max, -1/powerOfSubstitution);
			aux=-aux;
		}
		else return 0;
	}
	else if(powerOfSubstitution<0){
		if(max<0 && min<0){
			this->min=- pow(-max, 1/powerOfSubstitution);
			this->max=- pow(-min, 1/powerOfSubstitution);
			reverse=!reverse;
		}
		else if(min<0){
			this->min=- pow(-min, 1/powerOfSubstitution);
			this->max=0;
			aux=-aux;
		}
		else return 0;
	}
	else{
		this->min=min;
		this->max=max;
	}
	if(reverse) reverseX=this->min+this->max;
	this->function2use=function2use;
	this->powerOfSubstitution=powerOfSubstitution;
	if(randomRatio>1) randomRatio=1;
	randomNumber=randomRatio;
	result=rombergIntegrateOpened();
	if(randomNumber<0){
		randomNumber/=-fabs(result);
		if(randomNumber>1) randomNumber=1;
		if(randomNumber<0) randomNumber=0;
	}
	savedResult=result;
	randomDo=true;
	return aux*result;
}

//----------------------------------------------------------------------------------------------------//

/*
* using Newton's method refines the value of the upper limit
* that results in the ratio of integrals equal to randomNumber
*/

void Integral::refineUpperLimit(double result){
	int i, rombergStore;
	double maxStore, minStore, functionValue, f, df;
	double deltaX, deltaOld, currentX, aux;
	double xlow, xhi, flow, fhi;
	if(randomNumber==0 || randomNumber==1){ return; }
        xlow=min;
        xhi=max;
	flow=-randomNumber*result;
	fhi=(1-randomNumber)*result;
	if(flow*fhi>0){ printf("Error (in Integral/refineUpperLimit): Root must be bracketed"); return; }
	if(flow>0){
		aux=xlow; xlow=xhi; xhi=aux;
		aux=flow; flow=fhi; fhi=aux;
	}
	deltaX=deltaOld=max-min;
	if(randomX<min || randomX>max) randomX=(min+max)/2;
	currentX=randomX;
	rombergStore=romberg;
	minStore=min;
	maxStore=max;
	max=randomX;
	functionValue=rombergIntegrateOpened(result)-randomNumber*result;
	romberg=romberg4refine;
	f=functionValue;
	df=function(currentX);
	for(i=0;i<maxSteps;i++){
		if(f<0){ xlow=currentX; flow=f; }
		else{ xhi=currentX; fhi=f; }
		if(((currentX-xhi)*df-f)*((currentX-xlow)*df-f)>0 || fabs(2*f)>fabs(deltaOld*df)){
			deltaOld=deltaX;
			deltaX=(xhi-xlow)/2;
			currentX=xlow+deltaX;
			if(xlow==currentX) break;
		}
		else{
			deltaOld=deltaX;
			deltaX=f/df;
			aux=currentX;
			currentX-=deltaX;
			if(aux==currentX) break;
		}
		if(df==0){ if(fabs(deltaX)<precision*(maxStore-minStore)) break; }
		else{ if(fabs(df*deltaX)<precision*fabs(result)) break; }
		min=randomX;
		max=currentX;
		if(min>max){ aux=min; min=max; max=aux; aux=-1; }
                else aux=1;
		f=functionValue+aux*rombergIntegrateOpened(result);
		df=function(currentX);
	}
        if(i==maxSteps) printf("Warning (in Integral/refineUpperLimit): Precision %f has not been reached after %i steps \n", precision, maxSteps);
	randomX=currentX;
	romberg=rombergStore;
	min=minStore;
	max=maxStore;
}

//----------------------------------------------------------------------------------------------------//

/**
* refines and returns the value of the upper limit x(rand)
*/

double Integral::getUpperLimit(){
	if(randomDo){
                if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION)
                    return min;
		refineUpperLimit(savedResult);
                if(reverse)
                    randomX=reverseX-randomX;
                if(powerOfSubstitution>0)
                    randomX= pow(randomX, -powerOfSubstitution);
                else if(powerOfSubstitution<0)
                    randomX=- pow(-randomX, powerOfSubstitution);
                if(useLog)
                    randomX= exp(randomX);
		return randomX;
	}
	else{
		printf("Error (in Integral/getUpperLimit): no previous call to upper limit functions was made");
		return 0;
	}
}

//----------------------------------------------------------------------------------------------------//

/**
* finds integral for opened intervals using log substitution
*/

//

double Integral::integrateWithLog(double min, double max, FunctionOfx *function2use){

   //cout<<"\t Starting to calc integratewith log with min="<<min<<"   max="<<max<<endl;

	double aux;

	reverse=false;

	useLog=true;

        //cout<<fabs(max-min)<<"\t"<<fabs(min)*COMPUTER_PRECISION<<endl;
        if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION)
		return 0;
	else if(min>max)
	{
		aux=min;
		min=max;
		max=aux;
		aux=-1;

	}
	else
		aux=1;

	this->min= log(min);
	this->max= log(max);
      //  if(min<200&&max>39900)cout<<"min = "<<min<<"\t log(min) = "<<log(min)<<endl;
        this->function2use=function2use;
	powerOfSubstitution=0;
	randomNumber=0;
	randomDo=false;
    //    cout <<"\t aux="<<aux<<endl;
//cout <<"In Integral::IBntegrateWithLog:798\t rombergIntegrapteOpened: "<<endl;
        return aux*rombergIntegrateOpened();
}

//----------------------------------------------------------------------------------------------------//

/**
* finds integral for opened intervals
* and computes the value of the x(rand)
*/

double Integral::integrateWithLog(double min, double max, FunctionOfx *function2use, double randomRatio){
	double aux, result;
	reverse=false;
	useLog=true;
        if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION)
        {
            this->min=min;
            this->max=max;
            randomDo=true;
            return 0;
        }
        else if(min>max)
        {
            aux=min;
            min=max;
            max=aux;
            aux=-1;
            reverse=!reverse;
        }
	else aux=1;
	this->min= log(min);
	this->max= log(max);
        if(reverse)
        {
            reverseX=this->min+this->max;
        }
	this->function2use=function2use;
	powerOfSubstitution=0;
	if(randomRatio>1) randomRatio=1;
	randomNumber=randomRatio;
	result=rombergIntegrateOpened();
        //cout<<"In integrateWithLog after rombergIntegrateOpened() was called "<<endl;
	if(randomNumber<0){

		randomNumber/=-fabs(result);
                if(randomNumber>1)
                {
                    randomNumber=1;
                }
                if(randomNumber<0)
                {
                    randomNumber=0;
                }
	}
	savedResult=result;
	randomDo=true;
        //cout<<"In integrateWithLog befor return"<<endl;
	return aux*result;	
}

//----------------------------------------------------------------------------------------------------//

/**
* finds integral for opened intervals using substitution x -&gt; 1/log(x)^(powerOfSubstitution)
*/

double Integral::integrateWithLogSubstitution(double min, double max, FunctionOfx *function2use, double powerOfSubstitution){
	double aux;
	reverse=false;
	useLog=true;
        if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION) return 0;
	if(min<0 || max<0) return 0;
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; }
	else aux=1;
	if(powerOfSubstitution>0){
		if(max>1 && min>1){
			this->min= pow( log(max), -1/powerOfSubstitution);
			this->max= pow( log(min), -1/powerOfSubstitution);
		}
		else if(max>1){
			this->min=0;
			this->max= pow( log(max), -1/powerOfSubstitution);
			aux=-aux;
		}
		else return 0;
	}
	else if(powerOfSubstitution<0){
		if(max<1 && min<1){
			this->min=- pow(- log(max), 1/powerOfSubstitution);
			this->max=- pow(- log(min), 1/powerOfSubstitution);
		}
		else if(min<1){
			this->min=- pow(- log(min), 1/powerOfSubstitution);
			this->max=0;
			aux=-aux;
		}
		else return 0;
	}
	else{
		this->min= log(min);
		this->max= log(max);
	}
	this->function2use=function2use;
	this->powerOfSubstitution=powerOfSubstitution;
	randomNumber=0;
	randomDo=false;
	return aux*rombergIntegrateOpened();
}


