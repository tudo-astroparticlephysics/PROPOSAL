/*
 * RootFinder.cxx
 *
 *  Created on: 23.04.2013
 *      Author: koehne
 */

// #include <stdio.h>
#include <cmath>
// #include <iostream>

#include "PROPOSAL/math/RootFinder.h"
#include "PROPOSAL/Output.h"

using namespace std;
using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double RootFinder::FindRoot(double min,
                            double max,
                            double startX,
                            boost::function<double (double)> function,
                            boost::function<double (double)> differentiated_function) const
{

    int i;
    double deltaX, deltaOld, currentX, aux;
    double f, df, fmin, fmax, result, xdiff;


    fmin    =   function(min);
    fmax    =   function(max);

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
        log_warn("Root must be bracketed");
        //return min;
    }

    if(fmin>0)
    {
        aux     =   min;
        min     =   max;
        max     =   aux;
        aux     =   fmin;
        fmin    =   fmax;
        fmax    =   aux;
    }

    result  =   fmax - fmin;
    xdiff   =   fabs(max-min);
    deltaX  =   deltaOld    =   xdiff;

    if(startX>1 || startX<0)
    {
        startX  =   0.5;
    }

    currentX    =   min*(1 - startX) + max*startX;
    f           =   function(currentX);
    df          =   differentiated_function(currentX);

    for(i=0; i<maxSteps_; i++)
    {
        if(f<0)
        {
            min     =   currentX;
            fmin    =   f;
        }

        else
        {
            max     =   currentX;
            fmax    =   f;
        }

        if(((currentX - max)*df - f)*((currentX - min)*df -f) > 0 || fabs(2*f) > fabs(deltaOld*df))
        {
            deltaOld    =   deltaX;
            deltaX      =   (max - min)/2;
            currentX    =   min + deltaX;

            if(min==currentX)
            {
                break;
            }
        }
        else
        {
            deltaOld    =   deltaX;
            deltaX      =   f/df;
            aux         =   currentX;
            currentX    -=  deltaX;

            if(aux==currentX)
            {
                break;
            }
        }

        if(fabs(f) < precision_*result)
        {
            break;
        }

        f   =   function(currentX);
        df  =   differentiated_function(currentX);
    }

    // cerr<<"Number of steps in RootFinder was "<<i<<endl;
    // Combine this with "your program" | " awk '/Number of steps in RootFinder/ {a[$(NF)]++} END {for(i in a) print i, a[i]}' |
    // sort -n | awk '{a+=$1*$2; b+=$2; print} END {print "Average is ", a/b}' " to optimize the settings

    if(i==maxSteps_)
    {
        log_warn("Precision %e has not been reached after %i steps",precision_, maxSteps_);
    }

    return currentX;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


RootFinder::RootFinder()
    :maxSteps_  ( 20 )
    ,precision_ ( 1e-6 )
{

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


RootFinder::RootFinder(const RootFinder &finder)
    :maxSteps_  ( finder.maxSteps_ )
    ,precision_ ( finder.precision_ )
{

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


RootFinder::RootFinder(int maxSteps, double precision)
{
    if(maxSteps<=0)
    {
        printf("Warning (in Integral/Integral/1): maxSteps = %d, must be >0, setting to 1",maxSteps);
        maxSteps    =   1;
    }

    if(precision<=0)
    {
        printf("Warning (in Integral/Integral/2): precision = %f must be > 0, setting to 1.e-6",precision);
        precision   =   1.e-6;
    }

    maxSteps_     =   maxSteps;
    precision_    =   precision;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


RootFinder& RootFinder::operator=(const RootFinder &finder)
{
    if (this != &finder)
    {
      RootFinder tmp(finder);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool RootFinder::operator==(const RootFinder &finder) const
{
    if (maxSteps_  != finder.maxSteps_ )    return false;
    if (precision_ != finder.precision_ )   return false;

    //else
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool RootFinder::operator!=(const RootFinder &finder) const
{
    return !(*this == finder);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void RootFinder::swap(RootFinder &finder)
{
    using std::swap;

    swap( maxSteps_ , finder.maxSteps_ );
    swap( precision_, finder.precision_ );

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void RootFinder::SetMaxSteps(int maxSteps) {
    maxSteps_ = maxSteps;
}

void RootFinder::SetPrecision(double precision) {
    precision_ = precision;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


RootFinder::~RootFinder(){}
