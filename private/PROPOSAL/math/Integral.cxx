/*! \file   Integral.cxx
*   \brief  Source file for the integration routines.
*
*   For more details see the class documentation.
*
*   \date   02.08.2010
*   \author Jan-Hendrik Koehne
*/

#include <cmath>
// #include <math.h>
#include <algorithm>
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/Output.h"
// #include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"

// using namespace std;
using namespace PROPOSAL;


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


//   finds integral: choose the in integration method with the last parameter
//   method = 1: IntegrateClosed
//   method = 2: IntegrateOpened
//   method = 3: IntegrateWithSubstitution
//   method = 4: IntegrateWithLog
//   method = 5: IntegrateWithLogSubstitution

double Integral::Integrate(double min, double max, boost::function<double (double)> integrand, int method, double powerOfSubstitution)
{
    if (min == 0. && max == 0. )
    {
        return 0.;
    }

    switch (method)
    {
        case 1: return IntegrateClosed(min,max,integrand);
        case 2: return IntegrateOpened(min,max,integrand);
        case 3: return IntegrateWithSubstitution(min,max,integrand,powerOfSubstitution);
        case 4:
            if (min <= 0. || max <= 0.)
            {
                return 0;
            }
            return IntegrateWithLog(min,max,integrand);
        case 5:
            if (min <= 0. || max <= 0.)
            {
                return 0;
            }
            return IntegrateWithLogSubstitution(min,max,integrand,powerOfSubstitution);
        default:
            log_fatal("Unknown integration method! 0 is returned!");
            return 0;

    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::IntegrateWithRandomRatio(double min, double max, boost::function<double (double)> integrand, int method, double randomRatio, double powerOfSubstitution)
{
    if (min == 0. && max == 0. )
    {
        return 0.;
    }

    switch (method)
    {
        case 3: return IntegrateWithSubstitution(min,max,integrand,powerOfSubstitution, randomRatio);
        case 4:
            if (min <= 0. || max <= 0.)
            {
                return 0;
            }
            return IntegrateWithLog(min,max,integrand, randomRatio);
        default:
            log_fatal("Unknown integration method! 0 is returned!");
            return 0;

    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::GetUpperLimit()
{

    if(randomDo_)
    {
        if(fabs(max_-min_)<=fabs(min_)*COMPUTER_PRECISION)
        {
            return min_;
        }

        RefineUpperLimit(savedResult_);

        if(reverse_)
        {
            randomX_ =   reverseX_ - randomX_;
        }

        if(powerOfSubstitution_>0)
        {
            randomX_ =   pow(randomX_, -powerOfSubstitution_);
        }
        else if(powerOfSubstitution_<0)
        {
            randomX_ =   -pow(-randomX_, powerOfSubstitution_);
        }

        if(useLog_)
        {
            randomX_ =   exp(randomX_);
        }

        return randomX_;
    }
    else{
        log_error("No previous call to upper limit functions was made!");
        return 0;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::IntegrateWithSubstitution(double min, double max, boost::function<double (double)> integrand, double powerOfSubstitution, double randomRatio)
{
    double aux, result;

    reverse_ =   false;
    useLog_  =   false;

    if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION)
    {
        this->min_   =   min;
        this->max_   =   max;
        randomDo_    =   true;
        return 0;
    }
    else if(min>max)
    {
        std::swap(min,max);
        aux     =   -1;
        reverse_ =   !reverse_;
    }
    else
    {
        aux=1;
    }

    if(powerOfSubstitution>0)
    {
        if(max>0 && min>0)
        {
            this->min_   =   pow(max, -1/powerOfSubstitution);
            this->max_   =   pow(min, -1/powerOfSubstitution);
            reverse_     =   !reverse_;
        }
        else if(max>0)
        {
            this->min_   =   0;
            this->max_   =   pow(max, -1/powerOfSubstitution);
            aux         =   -aux;
        }
        else
        {
            return 0;
        }
    }
    else if(powerOfSubstitution<0)
    {
        if(max<0 && min<0)
        {
            this->min_   =   -pow(-max, 1/powerOfSubstitution);
            this->max_   =   -pow(-min, 1/powerOfSubstitution);
            reverse_     =   !reverse_;
        }
        else if(min<0)
        {
            this->min_   =   -pow(-min, 1/powerOfSubstitution);
            this->max_   =   0;
            aux         =   -aux;
        }
        else
        {
            return 0;
        }
    }
    else
    {
            this->min_   =   min;
            this->max_   =   max;
    }

    if(reverse_)
    {
        reverseX_    =   this->min_ + this->max_;
    }

    this->integrand_            =   integrand;
    this->powerOfSubstitution_   =   powerOfSubstitution;

    if(randomRatio>1)
    {
        randomRatio =   1;
    }

    randomNumber_    =   randomRatio;
    result          =   RombergIntegrateOpened();

    if(randomNumber_<0)
    {
            randomNumber_    /=  -fabs(result);

            if(randomNumber_>1)
            {
                randomNumber_    =   1;
            }

            if(randomNumber_<0)
            {
                randomNumber_    =   0;
            }
    }

    savedResult_ =   result;
    randomDo_    =   true;

    return aux*result;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::IntegrateWithLog(double min, double max, boost::function<double (double)> integrand, double randomRatio)
{
    double aux, result;

    reverse_ =   false;
    useLog_  =   true;

    if(fabs(max-min)<=fabs(min)*COMPUTER_PRECISION)
    {
        this->min_   =   min;
        this->max_   =   max;
        randomDo_    =   true;

        return 0;
    }
    else if(min>max)
    {
        std::swap(min,max);
        aux     =   -1;
        reverse_ =   !reverse_;
    }
    else
    {
        aux=1;
    }

    this->min_   =   log(min);
    this->max_   =   log(max);

    if(reverse_)
    {
        reverseX_    =   this->min_ + this->max_;
    }

    this->integrand_    =   integrand;
    powerOfSubstitution_ =   0;

    if(randomRatio>1)
    {
        randomRatio=1;
    }

    randomNumber_        =   randomRatio;
    result              =   RombergIntegrateOpened();

    if(randomNumber_<0)
    {
        randomNumber_        /=  -fabs(result);

        if(randomNumber_>1)
        {
            randomNumber_    =   1;
        }

        if(randomNumber_<0)
        {
            randomNumber_    =   0;
        }
    }

    savedResult_ =   result;
    randomDo_    =   true;
    return aux*result;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Integral::Integral()
    :maxSteps_   (20)
    ,romberg_    (5)
    ,precision_  (1.e-6)
    ,max_        (1)
    ,min_        (0)
    ,iX_()
    ,iY_()
    ,c_()
    ,d_()
    ,romberg4refine_(2)
    ,powerOfSubstitution_(0)
    ,randomDo_(false)
    ,useLog_(false)
    ,randomNumber_(0)
    ,randomX_(0)
    ,reverse_(false)
    ,reverseX_(0)
    ,savedResult_(0)
{
    int aux;
    if(romberg_<=0)
    {
        log_warn("Warning (in Integral/Integral/0): romberg = %i must be > 0, setting to 1",romberg_);
        romberg_     =   1;
    }

    if(maxSteps_<=0)
    {
        log_warn("Warning (in Integral/Integral/1): maxSteps = %i must be > 0, setting to 1", maxSteps_);
        maxSteps_    =   1;
    }

    if(precision_<=0)
    {
        log_warn("Warning (in Integral/Integral/2): precision = %f must be > 0, setting to 1.e-6",precision_);
        precision_   =   1.e-6;
    }

    iX_.resize(maxSteps_);
    iY_.resize(maxSteps_);

    aux=std::max(romberg_, romberg4refine_);
    c_.resize(aux);
    d_.resize(aux);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Integral::Integral(const Integral &integral)
    :maxSteps_(integral.maxSteps_)
    ,romberg_(integral.romberg_)
    ,precision_(integral.precision_)
    ,max_(integral.max_)
    ,min_(integral.min_)
    ,iX_(integral.iX_)
    ,iY_(integral.iY_)
    ,c_(integral.c_)
    ,d_(integral.d_)
    ,romberg4refine_(integral.romberg4refine_)
    ,powerOfSubstitution_(integral.powerOfSubstitution_)
    ,randomDo_(integral.randomDo_)
    ,useLog_(integral.useLog_)
    ,randomNumber_(integral.randomNumber_)
    ,randomX_(integral.randomX_)
    ,reverse_(integral.reverse_)
    ,reverseX_(integral.reverseX_)
    ,savedResult_(integral.savedResult_)
{
    integrand_ = boost::ref(integral.integrand_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Integral::Integral(int  romberg, int maxSteps, double precision)
    :max_          (1)
    ,min_          (0)
    ,iX_()
    ,iY_()
    ,c_()
    ,d_()
    ,romberg4refine_(2)
    ,powerOfSubstitution_(0)
    ,randomDo_(false)
    ,useLog_(false)
    ,randomNumber_(0)
    ,randomX_(0)
    ,reverse_(false)
    ,reverseX_(0)
    ,savedResult_(0)
{
    int aux;
    if(romberg<=0)
    {
        log_warn("Warning (in Integral/Integral/0): romberg = %i must be > 0, setting to 1",romberg);
        romberg     =   1;
    }

    if(maxSteps<=0)
    {
        log_warn("Warning (in Integral/Integral/1): maxSteps = %i must be > 0, setting to 1", maxSteps);
        maxSteps    =   1;
    }

    if(precision<=0)
    {
        log_warn("Warning (in Integral/Integral/2): precision = %f must be > 0, setting to 1.e-6",precision);
        precision   =   1.e-6;
    }

    this->romberg_   =   romberg;
    this->maxSteps_  =   maxSteps;
    this->precision_ =   precision;

    iX_.resize(maxSteps_);
    iY_.resize(maxSteps_);

    aux=std::max(romberg_, romberg4refine_);
    c_.resize(aux);
    d_.resize(aux);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Integral& Integral::operator=(const Integral &integral)
{
    if (this != &integral)
    {
      Integral tmp(integral);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Integral::operator==(const Integral &integral) const
{
   // if(integrand_ != integral.integrand_)     return false;
    if(iX_.size()   !=  integral.iX_.size())    return false;
    if(iY_.size()   !=  integral.iY_.size())    return false;
    if(c_.size()    !=  integral.c_.size())     return false;
    if(d_.size()    !=  integral.d_.size())     return false;


    for(unsigned int i = 0; i < iX_.size(); i++)
    {
        if(iX_.at(i)    !=  integral.iX_.at(i)) return false;
    }
    for(unsigned int i = 0; i < iY_.size(); i++)
    {
        if(iY_.at(i)    !=  integral.iY_.at(i)) return false;
    }
    for(unsigned int i = 0; i < c_.size(); i++)
    {
        if(c_.at(i)     !=  integral.c_.at(i))  return false;
    }
    for(unsigned int i =    0; i < d_.size(); i++)
    {
        if(d_.at(i)     !=  integral.d_.at(i))  return false;
    }
    if(maxSteps_        != integral.maxSteps_)      return false;
    if(romberg_         != integral.romberg_)       return false;
    if(precision_       != integral.precision_)     return false;
    if(max_             != integral.max_)           return false;
    if(min_             != integral.min_)           return false;
    if(romberg4refine_  != integral.romberg4refine_)return false;
    if(randomDo_        != integral.randomDo_)      return false;
    if(useLog_          != integral.useLog_)        return false;
    if(randomNumber_    != integral.randomNumber_)  return false;
    if(randomX_         != integral.randomX_)       return false;
    if(reverse_         != integral.reverse_)       return false;
    if(reverseX_        != integral.reverseX_)      return false;
    if(savedResult_     != integral.savedResult_)   return false;

    if(powerOfSubstitution_ != integral.powerOfSubstitution_) return false;

    else return true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Integral::operator!=(const Integral &integral) const
{
  return !(*this == integral);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------/


void Integral::swap(Integral &integral)
{
    using std::swap;

    swap(maxSteps_,integral.maxSteps_);
    swap(romberg_,integral.romberg_);
    swap(precision_,integral.precision_);
    swap(max_,integral.max_);
    swap(min_,integral.min_);

    iX_.swap(integral.iX_);
    iY_.swap(integral.iY_);

    c_.swap(integral.c_);
    d_.swap(integral.d_);

    integrand_ = boost::ref(integral.integrand_);


    swap(romberg4refine_,integral.romberg4refine_);
    swap(powerOfSubstitution_,integral.powerOfSubstitution_);
    swap(randomDo_,integral.randomDo_);
    swap(useLog_,integral.useLog_);
    swap(randomNumber_,integral.randomNumber_);
    swap(randomX_,integral.randomX_);
    swap(reverse_,integral.reverse_);
    swap(reverseX_,integral.reverseX_);
    swap(savedResult_,integral.savedResult_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::Function(double x)
{
    double result, t;
    if(reverse_)
    {
        x       =   reverseX_ - x;
    }

    if(powerOfSubstitution_==0)
    {
        t       =   x;
        result  =   1;
    }
    else if(powerOfSubstitution_>0)
    {
        t       =   pow(x, -powerOfSubstitution_);
        result  =   powerOfSubstitution_*(t/x);
    }
    else
    {
        t       =   -pow(-x, powerOfSubstitution_);
        result  =   -powerOfSubstitution_*(t/x);
    }

    if(useLog_)
    {
        t       =   exp(t);
        result  *=  t;
    }
    result  *=  integrand_(t);

    if(result!=result)
    {
        if(integrand_(t) == 0 )
        {
            log_info("substitution not suitable! returning 0!");
            return 0;
        }
        else
        {
        log_fatal("result is nan! returning 0");
        return result;
        }
    }
    return result;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::Trapezoid(int n, double oldSum)
{
    double xStep, stepSize, resultSum;

    if(n==1)
    {
        return (Function(max_)+Function(min_))*(max_-min_)/2;
    }

    n           /=  2;
    stepSize    =   (max_-min_)/n;
    resultSum   =   0;

    for(xStep=min_+stepSize/2 ; xStep<max_ ; xStep+=stepSize)
    {
        resultSum   +=  Function(xStep);
    }

    return (oldSum+resultSum*stepSize)/2;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



double Integral::Trapezoid3(int n, double oldSum)
{
    double xStep, stepSize, resultSum;

    if(n==1)
    {
        return (max_-min_)*Function((max_+min_)/2);
    }

    stepSize    =   (max_-min_)/n;
    resultSum   =   0;

    for(xStep=min_+stepSize/2 ; xStep<max_ ; xStep+=stepSize)
    {
        resultSum   +=  Function(xStep);
        xStep       +=  2*stepSize;
        resultSum   +=  Function(xStep);
    }

    return oldSum/3+resultSum*stepSize;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::Trapezoid3S(int n, double oldSum, int stepNumber)
{
    double xStep, stepSize, resultSum;
    double smallSum=0, approX=0, functionValue1, functionValue2, sumDifference;
    double functionDifference, functionSum, aEq, bEq, bEq2, determinant;
    bool flag;

    if(n==1)
    {
        return (max_-min_)*Function((max_+min_)/2);
    }

    stepSize    =   (max_-min_)/n;
    if(stepNumber>=romberg_-1)
    {
        if(randomNumber_>=0)
        {
            smallSum    =   randomNumber_*oldSum/(1.5*stepSize);
        }
        else
        {
            smallSum    =-  randomNumber_/(1.5*stepSize);
            if(oldSum<0)
            {
                smallSum    *=  -1;
            }
        }
    }
    resultSum   =   0;
    flag        =   false;
    for(xStep=min_+stepSize/2 ; xStep<max_ ; xStep+=stepSize)
    {
        resultSum   +=  functionValue1=Function(xStep);
        xStep       +=  2*stepSize;
        resultSum   +=  functionValue2=Function(xStep);

        if(!flag) if(stepNumber>=romberg_-1) if((resultSum>=smallSum && smallSum>0) || (resultSum<=smallSum && smallSum<0)){
                functionSum         =   functionValue1+functionValue2;

                sumDifference       =   (smallSum-(resultSum-functionSum));
                sumDifference       *=  1.5*stepSize;

                functionDifference  =   functionValue2 - functionValue1;

                aEq     =   functionDifference/stepSize;
                bEq     =   (functionValue2 - 5*functionValue1)/2;
                bEq2    =   bEq*bEq;

                if(fabs(aEq*sumDifference)<precision_*bEq2)
                {
                    approX  =   sumDifference*2/functionSum;
                }
                else
                {
                    determinant =   bEq2 + 4*aEq*sumDifference;
                    if(determinant>=0)
                    {
                            determinant =   sqrt(determinant);
                            approX      =   (bEq + determinant)/aEq;

                            if(approX<0 || approX>3*stepSize)
                            {
                                approX  =   (bEq-determinant)/aEq;
                            }
                            else if(approX<0 || approX>3*stepSize)
                            {
                                approX  =   sumDifference*2/functionSum;
                            }
                    }
                    else
                    {
                            log_warn("Warning (in Integral/trapezoid3S): Determinant is negative, proceeding with linear approximation");
                            approX  =   sumDifference*2/functionSum;
                    }
                }

                approX  +=  xStep - 2.5*stepSize;
                flag    =   true;
        }
    }

    if(stepNumber>=romberg_-1)
    {
            if(!flag)
            {
                approX  =   max_;
            }
            randomX_ =   approX;
    }
    return oldSum/3+resultSum*stepSize;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Integral::InterpolationResults Integral::Interpolate(int start, double x)
{
    int num, i, k;
    bool dd;
    double error=0, result=0;
    double aux, aux2, dx1, dx2;
    Integral::InterpolationResults interpolation_results;

    num =   0;
    aux =   fabs(x-iX_[start+0]);

    for(i=0 ; i<romberg_ ; i++)
    {
        aux2    =   fabs(x-iX_[start+i]);
        if(aux2<aux)
        {
            num =   i;
            aux =   aux2;
        }
        c_[i]    =   iY_[start+i];
        d_[i]    =   iY_[start+i];
    }

    if(num==0)
    {
        dd=true;
    }
    else if(num==romberg_-1)
    {
        dd=false;
    }
    else
    {
        if(fabs(x-iX_[start+num-1]) > fabs(x-iX_[start+num+1]))
        {
            dd=true;
        }
        else
        {
            dd=false;
        }
    }

    result  =   iY_[start+num];
    for(k=1 ; k<romberg_ ; k++)
    {
        for(i=0 ; i<romberg_-k ; i++)
        {
            dx1     =   iX_[start+i]-x;
            dx2     =   iX_[start+i+k]-x;

            aux     =   c_[i+1] - d_[i];
            aux2    =   dx1 - dx2;

            if(aux2!=0)
            {
                aux     =   aux/aux2;
                c_[i]    =   dx1*aux;
                d_[i]    =   dx2*aux;
            }
            else
            {
                c_[i]    =   0;
                d_[i]    =   0;
            }
        }

        if(num==0)
        {
            dd  =   true;
        }

        if(num==romberg_-k)
        {
            dd  =   false;
        }

        if(dd)
        {
            error   =   c_[num];
        }
        else
        {
            num--;
            error   =   d_[num];
        }

        dd      =   !dd;
        result  +=  error;
    }

    interpolation_results.Error = error;
    interpolation_results.Value = result;
    return interpolation_results;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::RombergIntegrateClosed()
{
    double k = 1;
    double n = 1;
    double error, result, value;
    Integral::InterpolationResults interpolation_results;


    value   =   0;
    result  =   0;

    for(int i=0 ; i<maxSteps_ ; i++)
    {
        result  =   Trapezoid(k, result);
        iX_[i]   =   n;
        iY_[i]   =   result;

        if(i>=romberg_-1)
        {
            interpolation_results = Interpolate(i-(romberg_-1), 0);

            error   =   interpolation_results.Error;
            value   =   interpolation_results.Value;

            if(value!=0)
            {
                error   /=  value;
            }

            if(fabs(error)<precision_)
            {
                return value;
            }
        }

        k   =   k*2;
        n   =   n/4;
    }

    log_warn("Precision %f has not been reached after %i steps! Returning %f!", precision_, maxSteps_,value);
    return value;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::RombergIntegrateOpened()
{
    int i, k;
    double n;
    double error, result, value;
    Integral::InterpolationResults interpolation_results;

    k       =   1;
    n       =   1;
    value   =   0;
    result  =   0;
    for(i=0 ; i<maxSteps_ ; i++)
    {
        if(randomNumber_==0 || randomNumber_==1)
        {
            result  =   Trapezoid3(k, result);
            if(result != result){
                log_error("Function Value of is nan! Returning 0!");
                return 0;
            }

            if(randomNumber_==0)
            {
                randomX_=min_;
            }
            else
            {
                randomX_=max_;
            }
        }
        else
        {
            result=Trapezoid3S(k, result, i);
        }

        iX_[i]   =   n;
        iY_[i]   =   result;
        if(i>=romberg_-1)
        {
            interpolation_results = Interpolate(i-(romberg_-1), 0);
            error = interpolation_results.Error;
            value = interpolation_results.Value;

            if(value!=0)
            {
                error/=value;
            }

            if(fabs(error)<precision_)
            {
                return value;
            }
        }

        k   *=  3;
        n   /=  9;
    }

    log_warn("Precision %f has not been reached after %i steps! Returning %f!", precision_, maxSteps_,value);
    return value;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::RombergIntegrateOpened(double bigValue)
{
    int i, k;
    double n;
    double error, result, value;
    InterpolationResults interpolation_results;

    k       =   1;
    n       =   1;
    value   =   0;
    result  =   0;

    for(i=0 ; i<maxSteps_ ; i++)
    {
        result  =   Trapezoid3(k, result);
        iX_[i]   =   n;
        iY_[i]   =   result;
        if(i>=romberg_-1)
        {
            interpolation_results = Interpolate(i-(romberg_-1), 0);
            error = interpolation_results.Error;
            value = interpolation_results.Value;
            error   /=  bigValue;

            if(fabs(error)<precision_)
            {
                return value;
            }
        }
        k   *=  3;
        n   /=  9;
    }

    log_warn("Precision %f has not been reached after %i steps! Returning %f!", precision_, maxSteps_,value);
    return value;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::InitIntegralOpenedAndClosed(double min, double max, boost::function<double (double)> integrand)
{
    double aux;

    reverse_ =   false;
    useLog_  =   false;

    if(min>max)
    {
        std::swap(min,max);
        aux =   -1;
    }
    else
    {
        aux=1;
    }

    this->min_           =   min;
    this->max_           =   max;
    integrand_           =   integrand;
    powerOfSubstitution_ =   0;
    randomDo_            =   false;

    return aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::IntegrateClosed(double min, double max, boost::function<double (double)> integrand)
{
    double aux;
    aux = InitIntegralOpenedAndClosed(min, max, integrand);

    if(fabs(max_-min_)<=fabs(min_)*COMPUTER_PRECISION)
    {
        return 0;
    }
    return aux*RombergIntegrateClosed();
}


//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::IntegrateOpened(double min, double max, boost::function<double (double)> integrand)
{
    double aux;
    aux = InitIntegralOpenedAndClosed(min, max, integrand);

    if(fabs(max_-min_)<=fabs(min_)*COMPUTER_PRECISION)
    {
        return 0;
    }

    return aux*RombergIntegrateOpened();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::InitIntegralWithSubstitution(double min, double max, boost::function<double (double)> integrand, double powerOfSubstitution)
{
    double aux;

    reverse_ =   false;
    useLog_  =   false;

    if(min>max)
    {
        std::swap(min,max);
        aux =   -1;
    }
    else
    {
        aux=1;
    }

    if(powerOfSubstitution>0)
    {
        if(max>0 && min>0)
        {
            this->min_   =   pow(max, -1/powerOfSubstitution);
            this->max_   =   pow(min, -1/powerOfSubstitution);
        }
        else if(max>0)
        {
            this->min_   =   0;
            this->max_   =   pow(max, -1/powerOfSubstitution);
            aux         =   -aux;
        }
        else
        {
            return 0;
        }
    }
    else if(powerOfSubstitution<0)
    {
        if(max<0 && min<0){
            this->min_ = -pow(-max, 1/powerOfSubstitution);
            this->max_ = -pow(-min, 1/powerOfSubstitution);
        }
        else if(min<0){
            this->min_ = -pow(-min, 1/powerOfSubstitution);
            this->max_=0;
            aux=-aux;
        }
        else
        {
            return 0;
        }
    }
    else
    {
            this->min_   =   min;
            this->max_   =   max;
    }

    this->integrand_             =   integrand;
    this->powerOfSubstitution_   =   powerOfSubstitution;
    randomNumber_                =   0;
    randomDo_                    =   false;

    return aux;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::IntegrateWithSubstitution(double min, double max, boost::function<double (double)> integrand, double powerOfSubstitution)
{
    double aux;

    aux = InitIntegralWithSubstitution(min , max , integrand ,powerOfSubstitution);

    if(fabs(max_-min_)<=fabs(min_)*COMPUTER_PRECISION)
    {
        return 0;
    }


    return aux*RombergIntegrateOpened();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Integral::RefineUpperLimit(double result)
{
    int i, rombergStore;
    double maxStore, minStore, functionValue, f, df;
    double deltaX, deltaOld, currentX, aux;
    double xlow, xhi, flow, fhi;

    //This was commented out but is necessary since
    //it ensures to find a value if the searched integral value
    //is out of integral range.
    if(randomNumber_==0 || randomNumber_==1)
    {
        return;
    }

    xlow    =   min_;
    xhi     =   max_;
    flow    =   -randomNumber_*result;
    fhi     =   (1-randomNumber_)*result;

    if(flow*fhi>0)
    {
        log_error("Root must be bracketed");
        return;
    }

    if(flow>0)
    {
            std::swap(xlow,xhi);
            std::swap(flow,fhi);
    }

    deltaX      =   max_-min_;
    deltaOld    =   max_-min_;

    if(randomX_<min_ || randomX_>max_)
    {
        randomX_ =   (min_ + max_)/2;
    }
    currentX        =   randomX_;
    rombergStore    =   romberg_;
    minStore        =   min_;
    maxStore        =   max_;
    max_             =   randomX_;

    functionValue   =   RombergIntegrateOpened(result) - randomNumber_*result;

    romberg_ =   romberg4refine_;
    f       =   functionValue;
    df      =   Function(currentX);

    for(i=0 ; i<maxSteps_ ; i++)
    {
        if(f<0)
        {
            xlow    =   currentX;
            flow    =   f;
        }
        else
        {
            xhi     =   currentX;
            fhi     =   f;
        }

        if(((currentX-xhi)*df-f)*((currentX-xlow)*df-f)>0 || fabs(2*f)>fabs(deltaOld*df))
        {
            deltaOld    =   deltaX;
            deltaX      =   (xhi-xlow)/2;
            currentX    =   xlow+deltaX;

            if(xlow==currentX)
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

        if(df==0)
        {
            if(fabs(deltaX)<precision_*(maxStore-minStore))
            {
                break;
            }
        }
        else
        {
            if(fabs(df*deltaX)<precision_*fabs(result))
            {
                break;
            }
        }

        min_ =   randomX_;
        max_ =   currentX;

        if(min_>max_)
        {
            std::swap(min_,max_);
            aux=-1;
        }
        else
        {
            aux=1;
        }

        f   =   functionValue + aux*RombergIntegrateOpened(result);
        df  =   Function(currentX);
    }

    if(i==maxSteps_)
    {
        log_warn("Precision %f has not been reached after %i steps!", precision_, maxSteps_);
    }

    randomX_ =   currentX;
    romberg_ =   rombergStore;
    min_     =   minStore;
    max_     =   maxStore;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::InitIntegralWithLog(double min, double max, boost::function<double (double)> integrand)
{
    double aux;

    reverse_ =   false;
    useLog_  =   true;

    if(min>max)
    {
        std::swap(min,max);
        aux =   -1;

    }
    else
    {
        aux=1;
    }

    this->min_           =   log(min);
    this->max_           =   log(max);
    this->integrand_    =   integrand;
    powerOfSubstitution_ =   0;
    randomNumber_        =   0;
    randomDo_            =   false;

    return aux;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::IntegrateWithLog(double min, double max, boost::function<double (double)> integrand)
{
    double aux;

    aux =InitIntegralWithLog(min, max, integrand);

    if(fabs(max_-min_)<=fabs(min_)*COMPUTER_PRECISION)
    {
        return 0;
    }

    return aux* RombergIntegrateOpened();

}


//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



double Integral::InitIntegralWithLogSubstitution(double min, double max, boost::function<double (double)> integrand, double powerOfSubstitution)
{
    double aux;

    reverse_ =   false;
    useLog_  =   true;

    if(min<0 || max<0)
    {
        return 0;
    }
    else if(min>max)
    {
        std::swap(min,max);
        aux=-1;
    }
    else
    {
        aux=1;
    }

    if(powerOfSubstitution>0)
    {
        if(max>1 && min>1)
        {
            this->min_   =   pow( log(max), -1/powerOfSubstitution);
            this->max_   =   pow( log(min), -1/powerOfSubstitution);
        }
        else if(max>1){
            this->min_   =   0;
            this->max_   =   pow( log(max), -1/powerOfSubstitution);
            aux         =   -aux;
        }
        else
        {
            return 0;
        }
    }
    else if(powerOfSubstitution<0)
    {
        if(max<1 && min<1)
        {
            this->min_   =   -pow(- log(max), 1/powerOfSubstitution);
            this->max_   =   -pow(- log(min), 1/powerOfSubstitution);
        }
        else if(min<1)
        {
            this->min_   =   -pow(- log(min), 1/powerOfSubstitution);
            this->max_   =   0;
            aux         =   -aux;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        this->min_   =   log(min);
        this->max_   =   log(max);
    }

    this->integrand_             =   integrand;
    this->powerOfSubstitution_   =   powerOfSubstitution;
    randomNumber_                =   0;
    randomDo_                    =   false;

    return aux;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Integral::IntegrateWithLogSubstitution(double min, double max, boost::function<double (double)> integrand, double powerOfSubstitution)
{

    double aux;

    aux = InitIntegralWithLogSubstitution(min, max, integrand, powerOfSubstitution);

    if(fabs(max_-min_)<=fabs(min_)*COMPUTER_PRECISION)
    {
        return 0;
    }

    return aux*RombergIntegrateOpened();
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Integral::SetIntegrand(boost::function<double(double)> integrand) {
	integrand_ = integrand;
}

void Integral::SetMax(double max) {
	max_ = max;
}

void Integral::SetMaxSteps(int maxSteps) {
	maxSteps_ = maxSteps;
}

void Integral::SetMin(double min) {
	min_ = min;
}

void Integral::SetPowerOfSubstitution(double powerOfSubstitution) {
	powerOfSubstitution_ = powerOfSubstitution;
}

void Integral::SetPrecision(double precision) {
	precision_ = precision;
}

void Integral::SetRandomDo(bool randomDo) {
	randomDo_ = randomDo;
}

void Integral::SetRandomNumber(double randomNumber) {
	randomNumber_ = randomNumber;
}

void Integral::SetRandomX(double randomX) {
	randomX_ = randomX;
}

void Integral::SetReverse(bool reverse) {
	reverse_ = reverse;
}

void Integral::SetReverseX(double reverseX) {
	reverseX_ = reverseX;
}

void Integral::SetRomberg(int romberg) {
	romberg_ = romberg;
}

void Integral::SetRomberg4refine(int romberg4refine) {
	romberg4refine_ = romberg4refine;
}

void Integral::SetSavedResult(double savedResult) {
	savedResult_ = savedResult;
}

void Integral::SetUseLog(bool useLog) {
	useLog_ = useLog;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Integral::~Integral()
{

}
