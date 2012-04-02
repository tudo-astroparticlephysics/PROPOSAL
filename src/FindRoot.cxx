/*
 * FindRoot.cxx
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */

/**
 * Original comment of the Java class by Dmitry Chirkin:
 * This class provides routines for calculating roots by the combination of the Newton-Raphson method and bisection.
 * Include the function to be integrated in a class that implements the interface FunctionOfx (defined below).
 * Methods contained here are based on the Numerical Recipes (W. H. Press et al.).
 * <pre>
 * interface DFunctionOfx extends FunctionOfx{
 *     double dFunction(double x);
 * }
 * </pre>
 * For the definition of interface FunctionOfx see the manual page for class Integral.
 * @author Dmitry Chirkin
 * Ported to C++ by Jens Dreyer
 */
#include<iostream>
#include "FindRoot.h"
#include <stdio.h>


using namespace std;



        FindRoot::FindRoot()
        {
            _maxSteps=20;
            _precision=1e-6;
        }


        FindRoot::FindRoot(int maxSteps,double precision)
        {
            if(maxSteps<=0)
            {
                printf("Warning (in Integral/Integral/1): maxSteps = %d, must be >0, setting to 1",maxSteps);
                maxSteps=1;
            }
            if(precision<=0)
            {
                printf("Warning (in Integral/Integral/2): precision = %f must be > 0, setting to 1.e-6",precision);
                precision=1.e-6;
            }
            this->_maxSteps=maxSteps;
            this->_precision=precision;

        }


        FindRoot::~FindRoot(){}


        double FindRoot::FindRoot::function(double x)
        {
            return function2use->function(x);
        }


        double FindRoot::dFunction(double x)
        {
            return function2use->dFunction(x);
        }


        //----------------------------------------------------------------------------------------------------//

        /**
         * returns the value of the root bracketed between min and max. Starting value of x is determined by 0&lt;=startX&lt;=1
         */

        double FindRoot::findRoot(double min, double max, double startX, DFunctionOfx *function2use, double rightSide)
        {
            int i;
            double deltaX, deltaOld, currentX, aux;
            double f, df, fmin, fmax, result, xdiff;

            this->function2use=function2use;

            fmin=function(min)-rightSide;
            fmax=function(max)-rightSide;

            if(fmin==0)
            {
                return min;
            }
            if(fmax==0)
            {
                return max;
            }
            if(fmin*fmax>0)
            {
                printf("Error (in FindRoot/FindRoot): Root must be bracketed");
                return min;
            }

            if(fmin>0)
            {
                aux=min;
                min=max;
                max=aux;
                aux=fmin;
                fmin=fmax;
                fmax=aux;
            }
            result=fmax-fmin;
            xdiff=fabs(max-min);
            deltaX=deltaOld=xdiff;
            if(startX>1 || startX<0)
            {
                startX=0.5;
            }
            currentX=min*(1-startX)+max*startX;
            f=function(currentX)-rightSide;
            df=dFunction(currentX);

            for(i=0;i<_maxSteps;i++)
            {

                // Output.err.println("x = "+currentX+" f = "+f+" dx = "+(max-min)+" df = "+(f/df));
                if(f<0)
                {
                    min=currentX;
                    fmin=f;
                }
                else
                {
                    max=currentX;
                    fmax=f;
                }
                if(((currentX-max)*df-f)*((currentX-min)*df-f)>0 || fabs(2*f)>fabs(deltaOld*df))
                {
                    deltaOld=deltaX;
                    deltaX=(max-min)/2;
                    currentX=min+deltaX;
                    if(min==currentX)
                    {
                        break;
                    }
                }
                else{
                    deltaOld=deltaX;
                    deltaX=f/df;
                    aux=currentX;
                    currentX-=deltaX;
                    if(aux==currentX)
                    {
                        break;
                    }
                }
                if(fabs(f)<_precision*result)
                {
                    break;
                }
                // if(Math.abs(df)*xdiff<precision*result){ if(Math.abs(deltaX)<precision*xdiff) break; }
                // else{ if(Math.abs(df*deltaX)<precision*result) break; }
                f=function(currentX)-rightSide;
                df=dFunction(currentX);
            }

            // Output.err.println("Number of steps in FindRoot was "+i);
            // Combine this with "your program" | " awk '/Number of steps in FindRoot/ {a[$(NF)]++} END {for(i in a) print i, a[i]}' |
            // sort -n | awk '{a+=$1*$2; b+=$2; print} END {print "Average is ", a/b}' " to optimize the settings

            if(i==_maxSteps)
            {
                cerr<<"Warning (in FindRoot/findRoot): Precision"<<_precision<<" has not been reached after " << _maxSteps<<" steps"<<endl;
            }
            return currentX;
        }



