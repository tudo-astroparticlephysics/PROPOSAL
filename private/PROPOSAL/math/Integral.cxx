/*! \file   Integral.cxx
*   \brief  Source file for the integration routines.
*
*   For more details see the class documentation.
*
*   \date   02.08.2010
*   \author Jan-Hendrik Koehne
*/

#include <cmath>
#include <numeric>
// #include <math.h>
#include <algorithm>
#include <limits>
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
    ,q_limit_(50)
    // ,q_epsabs_(1e-50)
    // ,q_epsrel_(1e-6)
    ,q_last_3_results_()
    ,q_rlist2_()
    ,q_iord_()
    ,q_alist_()
    ,q_blist_()
    ,q_elist_()
    ,q_rlist_()
    ,q_fv1_()
    ,q_fv2_()
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

    q_last_3_results_.resize(3);
    q_rlist2_.resize(q_limit_ + 2);
    q_alist_.resize(q_limit_);
    q_blist_.resize(q_limit_);
    q_elist_.resize(q_limit_);
    q_rlist_.resize(q_limit_);
    q_iord_.resize(q_limit_);
    q_fv1_.resize(10);
    q_fv2_.resize(10);
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
    ,q_limit_(integral.q_limit_)
    // ,q_epsabs_(integral.q_epsabs_)
    // ,q_epsrel_(integral.q_epsrel_)
    ,q_last_3_results_(integral.q_last_3_results_)
    ,q_rlist2_(integral.q_rlist2_)
    ,q_iord_(integral.q_iord_)
    ,q_alist_(integral.q_alist_)
    ,q_blist_(integral.q_blist_)
    ,q_elist_(integral.q_elist_)
    ,q_rlist_(integral.q_rlist_)
    ,q_fv1_(integral.q_fv1_)
    ,q_fv2_(integral.q_fv2_)
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
    ,q_limit_(50)
    // ,q_epsabs_(1e-50)
    // ,q_epsrel_(1e-6)
    ,q_last_3_results_()
    ,q_rlist2_()
    ,q_iord_()
    ,q_alist_()
    ,q_blist_()
    ,q_elist_()
    ,q_rlist_()
    ,q_fv1_()
    ,q_fv2_()
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

    q_last_3_results_.resize(3);
    q_rlist2_.resize(q_limit_ + 2);
    q_alist_.resize(q_limit_);
    q_blist_.resize(q_limit_);
    q_elist_.resize(q_limit_);
    q_rlist_.resize(q_limit_);
    q_iord_.resize(q_limit_);
    q_fv1_.resize(10);
    q_fv2_.resize(10);
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

    double q_value = qags();
    for(i=0 ; i<16 ; i++)
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
                std::cout << "k = " << k << ", n = " << n << std::endl;
                if(std::abs(value - q_value) > 1e-6)
                {
                    log_warn("RombergIntegrateOpened qags %f and romberg %f differs", q_value, value);
                    // std::cout << q_value << " qags" << std::endl;
                    // std::cout << value << " romberg" << std::endl;
                }
                return value;
            }
        }

        k   *=  3;
        n   /=  9;
    }
    std::cout.precision(16);
    std::cout << q_value << " qags" << std::endl;
    std::cout << value << " romberg" << std::endl;
    std::cout << "k = " << k << ", n = " << n << std::endl;
    // log_warn("RombergIntegrateOpened qags %f and romberg %f differs", q_value, value);
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
    double q_value = qags();

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
                if(std::abs(value - q_value) > 1e-6)
                {
                    log_warn("RombergIntegrateOpened_bigValue qags %f and romberg %f differs", q_value, value);
                    // std::cout << q_value << " qags" << std::endl;
                    // std::cout << value << " romberg" << std::endl;
                }
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



double Integral::qags()
{

    const double q_epmach_ = std::numeric_limits<double>::epsilon(); // machine epsilon
    const double q_uflow_ = std::numeric_limits<double>::min(); // smallest finite value
    const double q_oflow_ = std::numeric_limits<double>::max(); // largest finite value
    // output parameter
    double result = 0.0;
    double abserr = 0.0;
    int neval = 0;
    int ier = 0;
    double q_epsabs_ = 1.0e-50;
    double q_epsrel_ = 1.0e-6;

    if (q_limit_ < 1)
    {
        std::cout<< "abnormal return from dqags" << std::endl;
        std::cout<< " Error Number = " << ier << std::endl;
        return result;
    }

    if (q_epsabs_ < 0. && q_epsrel_ < 0.)
    {
        ier = 6;
        return result;
    }
    std::pair<Integral::InterpolationResults, Integral::InterpolationResults> qk21_output;
    // std::cout << "first approx" << std::endl;

    // First approximation to the integral.
    qk21_output = q_gaus_kronrod_21(min_, max_);
    result = qk21_output.first.Value;
    abserr = qk21_output.first.Error;
    double defabs = qk21_output.second.Value;
    double resabs = qk21_output.second.Error;

    // Test on accuracy.
    double dres = std::abs(result);
    double errbnd = std::max(q_epsabs_, q_epsrel_*dres);
    // std::cout << "abserr = " << abserr << ", errbnd = " << errbnd << ", resabs = " << resabs << std::endl;
    int last = 1;

    if(abserr <= 100*q_epmach_*defabs && abserr > errbnd)
        ier = 2;

    if(q_limit_ == 1)
        ier = 1;

    if(ier != 0 || (abserr <= errbnd && abserr != resabs ) || abserr == 0.0)
    {
        neval = 42 * last - 21;
        std::cout.precision(16);
        std::cout << "ier = " << ier << ", neval = " << neval << ", result = " << result << ", abserr = " << abserr << std::endl;
        return result;
    }

    int ierro = 0;
    q_alist_[0] = min_;
    q_blist_[0] = max_;
    q_elist_[0] = abserr;
    q_rlist_[0] = result;
    q_iord_[0] = 1;


    // local parameters
    double abseps,area1,area12,area2,a1,
        a2,b1,b2,correc,defab1,defab2,
        erlarg,erlast,
        error1,error2,erro12,ertest,reseps,
        small;
    int jupbnd;
    Integral::InterpolationResults extrapolation_output;
    // Initialization
    q_rlist2_[0] = result;
    double area = result;
    double errmax = abserr;
    double errsum = abserr;
    int maxerr = 1;
    int nrmax = 1;
    int nres = 0;
    int numrl2 = 2;
    int ktmin = 0;
    bool extrap = false;
    bool noext = false;
    int iroff1 = 0;
    int iroff2 = 0;
    int iroff3 = 0;

    int ksgn = -1;
    if(dres >= (1.0 - 50*q_epmach_) * defabs)
        ksgn = 1;

    abserr = q_oflow_;

    // main loop
    for (last = 2; last <= q_limit_; last++)
    {
        // Bisect the subinterval with the nrmax-th largest error estimate.
        a1 = q_alist_[maxerr-1];
        b1 = 0.5 * (q_alist_[maxerr-1] + q_blist_[maxerr-1]);
        a2 = b1;
        b2 = q_blist_[maxerr-1];
        erlast = errmax;

        qk21_output = q_gaus_kronrod_21(a1, b1);
        area1 = qk21_output.first.Value;
        error1 = qk21_output.first.Error;
        resabs = qk21_output.second.Value;
        defab1 = qk21_output.second.Error;
        qk21_output = q_gaus_kronrod_21(a2, b2);
        area2 = qk21_output.first.Value;
        error2 = qk21_output.first.Error;
        resabs = qk21_output.second.Value;
        defab2 = qk21_output.second.Error;

        // Improve previous approximations to integral and error and test for accuracy.
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - q_rlist_[maxerr-1];

        if(defab1 != error1 && defab2 != error2 )
        {
            if(std::abs(q_rlist_[maxerr-1] - area12) <= 1.0e-5 * std::abs(area12) && erro12 >= 0.99*errmax)
            {
                if(extrap)
                {
                    iroff2 = iroff2 + 1;
                }
                else
                {
                    iroff1 = iroff1 + 1;
                }
            }
            if ( last > 10 && erro12 > errmax )
                iroff3 = iroff3 + 1;
        }

        q_rlist_[maxerr-1] = area1;
        q_rlist_[last-1] = area2;
        errbnd = std::max(q_epsabs_, q_epsrel_*std::abs(area));

        // Test for roundoff error and eventually set error flag.

        if(iroff1 + iroff2 >= 10 || iroff3 >= 20)
            ier = 2;

        if(iroff2 >= 5)
            ierro = 3;

        // Set error flag in the case that the number of subintervals equals limit.
        if(last == q_limit_)
            ier = 1;

        // Set error flag in the case of bad integrand behavior at a point of the integration range.
        // if(std::max(std::abs(a1), std::abs(b2)) <= (1 + 1000*q_epmach_) * (std::abs(a2) + 1000*q_uflow_)) // diff between f90 & f77
        if(std::max(std::abs(a1), std::abs(b2)) <= (1 + 100*q_epmach_) * (std::abs(a2) + 1000*q_uflow_))
            ier = 4;

        // Append the newly-created intervals to the list.
        if ( error2 <= error1 )
        {
            q_alist_[last-1] = a2;
            q_blist_[maxerr-1] = b1;
            q_blist_[last-1] = b2;
            q_elist_[maxerr-1] = error1;
            q_elist_[last-1] = error2;
        }
        else
        {
            q_alist_[maxerr-1] = a2;
            q_alist_[last-1] = a1;
            q_blist_[last-1] = b1;
            q_rlist_[maxerr-1] = area2;
            q_rlist_[last-1] = area1;
            q_elist_[maxerr-1] = error2;
            q_elist_[last-1] = error1;
        }

        // Call QSORT to maintain the descending ordering in the list of error estimates 
        // and select the subinterval with nrmax-th largest error estimate (to be bisected next).
        nrmax = q_sort(last, maxerr, nrmax);
        maxerr = q_iord_[nrmax-1];
        errmax = q_elist_[maxerr-1];

        if (errsum <= errbnd)
        {
            // Compute global integral sum.
            result = std::accumulate(q_rlist_.begin(), q_rlist_.begin()+(last), 0.0);
            abserr = errsum;

            if(2 < ier)
                ier = ier - 1;

            neval = 42*last - 21;
            std::cout.precision(16);
            std::cout << "in loop " << last << std::endl;
            std::cout << "ier = " << ier << ", neval = " << neval << ", result = " << result << ", abserr = " << abserr << std::endl;
            return result;
        }

        if(ier != 0)
            break;

        if(last == 2)
        {
            small = std::abs(max_ - min_) * 0.375;
            erlarg = errsum;
            ertest = errbnd;
            q_rlist2_[1] = area;
            continue;
        }

        if(noext)
            continue;

        erlarg = erlarg - erlast;

        if(std::abs(b1 - a1) > small)
            erlarg = erlarg + erro12;

        // Test whether the interval to be bisected next is the smallest interval.
        if (!extrap)
        {
            if(std::abs(q_blist_[maxerr-1] - q_alist_[maxerr-1]) > small)
                continue;

            extrap = true;
            nrmax = 2;
        }

        // The smallest interval has the largest error.
        // Before bisecting decrease the sum of the errors over 
        // the larger intervals (erlarg) and perform extrapolation.
        if ( ierro != 3 && erlarg > ertest )
        {
            if(last > (2 + q_limit_/2))
            {
                jupbnd = q_limit_ + 3 - last;
            }
            else
            {
                jupbnd = last;
            }

            // int id = nrmax;
            // for (int idx = id; idx <= jupbnd; idx++)
            // {
            //     maxerr = q_iord_[nrmax-1];
            //     errmax = q_elist_[maxerr-1];
            //     if(std::abs(q_blist_[maxerr-1] - q_alist_[maxerr-1]) > small)
            //         continue;

            //     nrmax = nrmax+1;
            // }
            while(nrmax <= jupbnd)
            {
                maxerr = q_iord_[nrmax-1];
                errmax = q_elist_[maxerr-1];
                if(std::abs(q_blist_[maxerr-1] - q_alist_[maxerr-1]) > small)
                    continue;

                nrmax = nrmax+1;
            }
        }

        // Perform extrapolation.
        numrl2 = numrl2 + 1;
        q_rlist2_[numrl2-1] = area;
        nres += 1;
        extrapolation_output = q_epsilon_extrapolation(numrl2, nres); // f90 q_rlist2_ -> epstab
        reseps = extrapolation_output.Value;
        abseps = extrapolation_output.Error;
        ktmin = ktmin + 1;

        if(ktmin > 5 && abserr < 1.0e-3 * errsum)
            ier = 5;
        
        if (abseps < abserr)
        {
            ktmin = 0;
            abserr = abseps;
            result = reseps;
            correc = erlarg;
            ertest = std::max(q_epsabs_, q_epsrel_*std::abs(reseps));

            if (abserr <= ertest)
                break;
        }

        // Prepare bisection of the smallest interval.
        if (numrl2 == 1)
            noext = true;

        if (ier == 5)
            break;

        maxerr = q_iord_[0];
        errmax = q_elist_[maxerr-1];
        nrmax = 1;
        extrap = false;
        small = small * 0.5;
        erlarg = errsum;
    } // end of main loop

    // Set final result and error estimate.
    if (abserr == q_oflow_)
    {
        // Compute global integral sum.
        result = std::accumulate(q_rlist_.begin(), q_rlist_.begin()+(last), 0.0);
        abserr = errsum;
    }
    else
    {
        // Test on divergence.
        if (ier + ierro == 0)
        {
            if (ksgn != -1 || std::max(std::abs(result), std::abs(area)) >  defabs*0.01)
            {
                if(0.01 > (result/area) || (result/area) > 100 || errsum > std::abs(area))
                    ier = 6;
            }
        }

        if (ierro == 3)
            abserr = abserr + correc;

        if (ier == 0)
            ier = 3;

        if ( result != 0.0 && area != 0.0 )
        {
            if ( abserr/std::abs(result) > errsum/std::abs(area) )
            {
                // // Compute global integral sum.
                result = std::accumulate(q_rlist_.begin(), q_rlist_.begin()+(last), 0.0);
                abserr = errsum;
            }
        }
        else
        {
            if ( abserr > errsum )
            {
                // Compute global integral sum.
                result = std::accumulate(q_rlist_.begin(), q_rlist_.begin()+(last), 0.0);
                abserr = errsum;
            }
            else
            {
                if ( area != 0.0 )
                {
                    if (ksgn != -1 || std::max(std::abs(result), std::abs(area)) >  defabs*0.01)
                    {
                        if(0.01 > (result/area) || (result/area) > 100 || errsum > std::abs(area))
                            ier = 6;
                    }
                }
            }
        }
    }

    if ( 2 < ier )
        ier = ier - 1;

    neval = 42*last - 21;

    std::cout.precision(16);
    std::cout << "after loop " << std::endl;
    std::cout << "ier = " << ier << ", neval = " << neval << ", result = " << result << ", abserr = " << abserr << std::endl;
    return result;
}

std::pair<Integral::InterpolationResults, Integral::InterpolationResults> Integral::q_gaus_kronrod_21(
    double min_lim, 
    double max_lim)
{
    const double q_epmach_ = std::numeric_limits<double>::epsilon(); // machine epsilon
    const double q_uflow_ = std::numeric_limits<double>::min(); // smallest finite value
    const double q_weights_gaus_10p_ [5] =
    {
        0.066671344308688137593568809893332,
        0.149451349150580593145776339657697,
        0.219086362515982043995534934228163,
        0.269266719309996355091226921569469,
        0.295524224714752870173892994651338
    };
    const double q_weights_kronrod_21p_ [11] =
    {
        0.011694638867371874278064396062192,
        0.032558162307964727478818972459390,
        0.054755896574351996031381300244580,
        0.075039674810919952767043140916190,
        0.093125454583697605535065465083366,
        0.109387158802297641899210590325805,
        0.123491976262065851077958109831074,
        0.134709217311473325928054001771707,
        0.142775938577060080797094273138717,
        0.147739104901338491374841515972068,
        0.149445554002916905664936468389821
    };
    const double q_abscissae_kronrod_21p_[11] =
    {
        0.995657163025808080735527280689003,
        0.973906528517171720077964012084452,
        0.930157491355708226001207180059508,
        0.865063366688984510732096688423493,
        0.780817726586416897063717578345042,
        0.679409568299024406234327365114874,
        0.562757134668604683339000099272694,
        0.433395394129247190799265943165784,
        0.294392862701460198131126603103866,
        0.148874338981631210884826001129720,
        0.000000000000000000000000000000000
    };
    // output parameters
    double result, abserr, resabs, resasc;
    Integral::InterpolationResults qk21_output;
    Integral::InterpolationResults qk21_abs_output;
    // internal parameters
    double absc, fsum, fval1, fval2, reskh;
    int jtw, jtwm1;

    double centr = 0.5*(min_lim + max_lim);
    double hlgth = 0.5*(max_lim - min_lim);
    double dhlgth = std::abs(hlgth);

    // Compute the 21-point Kronrod approximation to the
    // integral, and estimate the absolute error.

    double resg = 0.0;
    double fc = Function(centr);
    double resk = q_weights_kronrod_21p_[10]*fc;
    resabs = std::abs(resk);

    for(int idx = 1; idx <= 5; idx++)
    {
        jtw = 2*idx;
        absc = hlgth*q_abscissae_kronrod_21p_[jtw-1];
        fval1 = Function(centr-absc);
        fval2 = Function(centr+absc);
        q_fv1_[jtw-1] = fval1;
        q_fv2_[jtw-1] = fval2;
        fsum = fval1 + fval2;
        resg = resg + q_weights_gaus_10p_[idx-1]*fsum;
        resk = resk + q_weights_kronrod_21p_[jtw-1]*fsum;
        resabs = resabs + q_weights_kronrod_21p_[jtw-1]*(std::abs(fval1) + std::abs(fval2));
    }

    for(int idx = 1; idx <= 5; idx++)
    {
        jtwm1 = 2*idx - 1;
        absc = hlgth*q_abscissae_kronrod_21p_[jtwm1-1];
        fval1 = Function(centr-absc);
        fval2 = Function(centr+absc);
        q_fv1_[jtwm1-1] = fval1;
        q_fv2_[jtwm1-1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + q_weights_kronrod_21p_[jtwm1-1]*fsum;
        resabs = resabs + q_weights_kronrod_21p_[jtwm1-1]*(std::abs(fval1) + std::abs(fval2));
    }

    reskh = resk*0.5;
    resasc = q_weights_kronrod_21p_[10]*std::abs(fc - reskh);

    for(int idx = 0; idx < 10; idx++)
    {
        resasc = resasc + q_weights_kronrod_21p_[idx]*(std::abs(q_fv1_[idx] - reskh) + std::abs(q_fv2_[idx] - reskh));
    }
    result = resk * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);

    if(resasc != 0.0 && abserr != 0.0)
        abserr = resasc * std::min(1., std::pow(200.0 * abserr / resasc, 1.5));
    // {
    //     double tmp_for_pow = std::sqrt(200. * abserr / resasc);
    //     abserr = resasc * std::min(1., tmp_for_pow*tmp_for_pow*tmp_for_pow);
    // }

    if(resabs > q_uflow_ / (50. * q_epmach_))
        abserr = std::max((q_epmach_ * 50.0)*resabs, abserr);

    qk21_output.Value = result;
    qk21_output.Error = abserr;
    qk21_abs_output.Value = resabs;
    qk21_abs_output.Error = resasc;
    return std::make_pair(qk21_output, qk21_abs_output);
}

int Integral::q_sort(int last, int maxerr, int nrmax)
{
    double errmax, errmin;
    int i,ibeg,isucc,jbnd,jupbn,k, nrmax_intern;

    nrmax_intern = nrmax;
    // Check whether the list contains more than two error estimates.
    if(last <= 2)
    {
        q_iord_[0] = 1;
        q_iord_[1] = 2;
        return nrmax_intern;
    }

    // This part of the routine is only executed if, due to a difficult integrand,
    // subdivision increased the error estimate. in the normal case the insert procedure should
    // start after the nrmax-th largest error estimate.

    errmax = q_elist_[maxerr-1];

    if(nrmax_intern != 1)
    {
        int max_iter = nrmax_intern - 1;
        for (int idx = 1; idx <= max_iter; idx++)
        {
            isucc = q_iord_[nrmax_intern-2];

            if(q_elist_[isucc-1] <= errmax)
                break;

            q_iord_[nrmax_intern-1] = isucc;
            nrmax_intern = nrmax_intern - 1;
        }
    }

    // Compute the number of elements in the list to be maintained in descending order.
    // This number depends on the number of subdivisions still allowed.

    jupbn = last;
    if((q_limit_/2 + 2) < last)
        jupbn = q_limit_ + 3 - last;

    errmin = q_elist_[last-1];

    // Insert errmax by traversing the list top-down, starting
    // comparison from the element elist(iord(nrmax+1)).

    jbnd = jupbn - 1;
    ibeg = nrmax_intern + 1;

    if(ibeg > jbnd)
    {
        q_iord_[jbnd-1] = maxerr;
        q_iord_[jupbn-1] = last;
        return nrmax_intern;
    }
    
    bool jump_to_60 = false;
    for(i = ibeg; i <= jbnd; i++)
    {
        isucc = q_iord_[i-1];
        if(q_elist_[isucc-1] <= errmax)
        {
            jump_to_60 = true;
            break;
        }

        q_iord_[i-2] = isucc;
    }

    if(!jump_to_60)
    {
        q_iord_[jbnd-1] = maxerr;
        q_iord_[jupbn-1] = last;
        return nrmax_intern;
    }

    // Insert errmin by traversing the list bottom-up.
    q_iord_[i-2] = maxerr;
    k = jbnd;

    for(int idx = i; idx <= jbnd; idx++)
    {
        isucc = q_iord_[k-1];
        if(errmin < q_elist_[isucc-1])
        {
            q_iord_[k] = last;
            return nrmax_intern;
        }
        q_iord_[k] = isucc;
        k = k-1;
    }
    q_iord_[i-1] = last;

    return nrmax_intern;
}

Integral::InterpolationResults Integral::q_epsilon_extrapolation(int numrl2, int nres)
{
    const double q_epmach_ = std::numeric_limits<double>::epsilon(); // machine epsilon
    const double q_oflow_ = std::numeric_limits<double>::max(); // largest finite value
    const int q_limit_epsilon_table_ = 50;
    // output
    Integral::InterpolationResults extrapolation_output;
    double result, abserr;
    abserr = q_oflow_;
    result = q_rlist2_[numrl2-1];
    if(numrl2 < 3)
    {
        abserr = std::max(abserr, 5.0 * q_epmach_ * std::abs(result));
        // abserr = std::max(abserr, 0.5 * q_epmach_ * std::abs(result)); diff between f90 & f77
        extrapolation_output.Value = result;
        extrapolation_output.Error = abserr;
        return extrapolation_output;
    }

    double delta1,delta2,delta3,
        epsinf,error,err1,err2,err3,e0,e1,e1abs,e2,e3,
        res,ss,tol1,tol2,tol3;
    int ib,ib2,indx,k1,k2,k3,newelm,num;

    q_rlist2_[numrl2+1] = q_rlist2_[numrl2-1];
    newelm = (numrl2-1) / 2;
    q_rlist2_[numrl2-1] = q_oflow_;
    num = numrl2;
    k1 = numrl2;

    for(int idx = 1; idx <= newelm; idx++)
    {
        k2 = k1 - 1;
        k3 = k1 - 2;
        res = q_rlist2_[k1+1];
        e0 = q_rlist2_[k3-1];
        e1 = q_rlist2_[k2-1];
        e2 = res;
        e1abs = std::abs(e1);
        delta2 = e2 - e1;
        err2 = std::abs(delta2);
        tol2 = std::max(std::abs(e2), e1abs) * q_epmach_;
        delta3 = e1 - e0;
        err3 = std::abs(delta3);
        tol3 = std::max(e1abs, std::abs(e0)) * q_epmach_;

        // If e0, e1 and e2 are equal to within machine accuracy, convergence is assumed.
        if(err2 <= tol2 && err3 <= tol3)
        {
            result = res;
            abserr = err2 + err3;

            abserr = std::max(abserr, 5.0 * q_epmach_ * std::abs(result));
            // abserr = std::max(abserr, 0.5 * q_epmach_ * std::abs(result)); diff between f90 & f77
            extrapolation_output.Value = result;
            extrapolation_output.Error = abserr;
            return extrapolation_output;
        }

        e3 = q_rlist2_[k1-1];
        q_rlist2_[k1-1] = e1;
        delta1 = e1 - e3;
        err1 = std::abs(delta1);
        tol1 = std::max(e1abs, std::abs(e3)) * q_epmach_;

        // If two elements are very close to each other, omit a part
        // of the table by adjusting the value of N.
        if(err1 <= tol1 || err2 <= tol2 || err3 <= tol3)
        {
            numrl2 = idx + idx - 1;
            break;
        }

        ss = 1.0 / delta1 + 1.0 / delta2 - 1.0 / delta3;
        epsinf = std::abs(ss * e1);

        if(epsinf <= 1.0e-4)
        {
            numrl2 = idx + idx - 1;
            break;
        }

        res = e1 + 1.0 / ss;
        q_rlist2_[k1-1] = res;
        k1 = k1 - 2;
        error = err2 + std::abs(res - e2) + err3;

        if(error <= abserr)
        {
            abserr = error;
            result = res;
        }
    }

    // Shift the table.
    if(numrl2 == q_limit_epsilon_table_)
        numrl2 = 2 * (q_limit_epsilon_table_ / 2) - 1;

    ib = 1;
    if((num / 2) * 2 == num )
        ib = 2;

    for(int idx = 1; idx <= newelm + 1; idx++)
    {
        ib2 = ib + 2;
        q_rlist2_[ib-1] = q_rlist2_[ib2-1];
        ib = ib2;
    }

    if(num != numrl2)
    {
        indx = num - numrl2 + 1;
        for(int idx = 1; idx <= numrl2; idx++)
        {
            q_rlist2_[idx-1] = q_rlist2_[indx-1];
            indx = indx + 1;
        }
    }

    if(nres < 4)
    {
      q_last_3_results_[nres-1] = result;
      abserr = q_oflow_;
    }
    else
    {
        abserr = std::abs(result - q_last_3_results_[3]) +
            std::abs(result - q_last_3_results_[2]) +
            std::abs(result - q_last_3_results_[1]);
        q_last_3_results_[1] = q_last_3_results_[2];
        q_last_3_results_[2] = q_last_3_results_[3];
        q_last_3_results_[3] = result;
    }

    abserr = std::max(abserr, 5.0 * q_epmach_ * std::abs(result));
    // abserr = std::max(abserr, 0.5 * q_epmach_ * std::abs(result)); diff between f90 & f77
    extrapolation_output.Value = result;
    extrapolation_output.Error = abserr;
    return extrapolation_output;
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
