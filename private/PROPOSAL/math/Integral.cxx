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
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/Integral.h"
#include <algorithm>
#include <limits>
// #include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"

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

double Integral::Integrate(double min,
                           double max,
                           std::function<double(double)> integrand,
                           int method,
                           double powerOfSubstitution)
{
    if (min == 0. && max == 0.)
    {
        return 0.;
    }

    switch (method)
    {
        case 1:
            return IntegrateClosed(min, max, integrand);
        case 2:
            return IntegrateOpened(min, max, integrand);
        case 3:
            return IntegrateWithSubstitution(min, max, integrand, powerOfSubstitution);
        case 4:
            if (min <= 0. || max <= 0.)
            {
                return 0;
            }
            return IntegrateWithLog(min, max, integrand);
        case 5:
            if (min <= 0. || max <= 0.)
            {
                return 0;
            }
            return IntegrateWithLogSubstitution(min, max, integrand, powerOfSubstitution);
        default:
            log_fatal("Unknown integration method! 0 is returned!");
            return 0;
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::IntegrateWithRandomRatio(double min,
                                          double max,
                                          std::function<double(double)> integrand,
                                          int method,
                                          double randomRatio,
                                          double powerOfSubstitution)
{
    if (min == 0. && max == 0.)
    {
        return 0.;
    }

    switch (method)
    {
        case 3:
            return IntegrateWithSubstitution(min, max, integrand, powerOfSubstitution, randomRatio);
        case 4:
            if (min <= 0. || max <= 0.)
            {
                return 0;
            }
            return IntegrateWithLog(min, max, integrand, randomRatio);
        default:
            log_fatal("Unknown integration method! 0 is returned!");
            return 0;
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::GetUpperLimit()
{

    if (randomDo_)
    {
        if (std::abs(max_ - min_) <= std::abs(min_) * COMPUTER_PRECISION)
        {
            return min_;
        }

        RefineUpperLimit(savedResult_);

        if (reverse_)
        {
            randomX_ = reverseX_ - randomX_;
        }

        if (powerOfSubstitution_ > 0)
        {
            randomX_ = std::pow(randomX_, -powerOfSubstitution_);
        } else if (powerOfSubstitution_ < 0)
        {
            randomX_ = -std::pow(-randomX_, powerOfSubstitution_);
        }

        if (useLog_)
        {
            randomX_ = std::exp(randomX_);
        }

        return randomX_;
    } else
    {
        log_error("No previous call to upper limit functions was made!");
        return 0;
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::IntegrateWithSubstitution(double min,
                                           double max,
                                           std::function<double(double)> integrand,
                                           double powerOfSubstitution,
                                           double randomRatio)
{
    double aux, result;

    reverse_ = false;
    useLog_  = false;

    if (std::abs(max - min) <= std::abs(min) * COMPUTER_PRECISION)
    {
        this->min_ = min;
        this->max_ = max;
        randomDo_  = true;
        return 0;
    } else if (min > max)
    {
        std::swap(min, max);
        aux      = -1;
        reverse_ = !reverse_;
    } else
    {
        aux = 1;
    }

    if (powerOfSubstitution > 0)
    {
        if (max > 0 && min > 0)
        {
            this->min_ = std::pow(max, -1 / powerOfSubstitution);
            this->max_ = std::pow(min, -1 / powerOfSubstitution);
            reverse_   = !reverse_;
        } else if (max > 0)
        {
            this->min_ = 0;
            this->max_ = std::pow(max, -1 / powerOfSubstitution);
            aux        = -aux;
        } else
        {
            return 0;
        }
    } else if (powerOfSubstitution < 0)
    {
        if (max < 0 && min < 0)
        {
            this->min_ = -std::pow(-max, 1 / powerOfSubstitution);
            this->max_ = -std::pow(-min, 1 / powerOfSubstitution);
            reverse_   = !reverse_;
        } else if (min < 0)
        {
            this->min_ = -std::pow(-min, 1 / powerOfSubstitution);
            this->max_ = 0;
            aux        = -aux;
        } else
        {
            return 0;
        }
    } else
    {
        this->min_ = min;
        this->max_ = max;
    }

    if (reverse_)
    {
        reverseX_ = this->min_ + this->max_;
    }

    this->integrand_           = integrand;
    this->powerOfSubstitution_ = powerOfSubstitution;

    if (randomRatio > 1)
    {
        randomRatio = 1;
    }

    randomNumber_ = randomRatio;
    result        = RombergIntegrateOpened();

    if (randomNumber_ < 0)
    {
        randomNumber_ /= -std::abs(result);

        if (randomNumber_ > 1)
        {
            randomNumber_ = 1;
        }

        if (randomNumber_ < 0)
        {
            randomNumber_ = 0;
        }
    }

    savedResult_ = result;
    randomDo_    = true;

    return aux * result;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::IntegrateWithLog(double min, double max, std::function<double(double)> integrand, double randomRatio)
{
    double aux, result;

    reverse_ = false;
    useLog_  = true;

    if (std::abs(max - min) <= std::abs(min) * COMPUTER_PRECISION)
    {
        this->min_ = min;
        this->max_ = max;
        randomDo_  = true;

        return 0;
    } else if (min > max)
    {
        std::swap(min, max);
        aux      = -1;
        reverse_ = !reverse_;
    } else
    {
        aux = 1;
    }

    this->min_ = std::log(min);
    this->max_ = std::log(max);

    if (reverse_)
    {
        reverseX_ = this->min_ + this->max_;
    }

    this->integrand_     = integrand;
    powerOfSubstitution_ = 0;

    if (randomRatio > 1)
    {
        randomRatio = 1;
    }

    randomNumber_ = randomRatio;
    result        = RombergIntegrateOpened();

    if (randomNumber_ < 0)
    {
        randomNumber_ /= -std::abs(result);

        if (randomNumber_ > 1)
        {
            randomNumber_ = 1;
        }

        if (randomNumber_ < 0)
        {
            randomNumber_ = 0;
        }
    }

    savedResult_ = result;
    randomDo_    = true;
    return aux * result;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Integral::Integral()
    : maxSteps_romberg_(12)
    , maxSteps_upper_limit_(20)
    , romberg_(5)
    , precision_(1.e-6)
    , max_(1)
    , min_(0)
    , iX_()
    , iY_()
    , c_()
    , d_()
    , romberg4refine_(2)
    , powerOfSubstitution_(0)
    , randomDo_(false)
    , useLog_(false)
    , randomNumber_(0)
    , randomX_(0)
    , reverse_(false)
    , reverseX_(0)
    , savedResult_(0)
    , q_last_3_results_()
    , q_rlist2_()
    , q_iord_()
{
    int aux;
    if (romberg_ <= 0)
    {
        log_warn("Warning (in Integral/Integral/0): romberg = %i must be > 0, setting to 1", romberg_);
        romberg_ = 1;
    }

    if (maxSteps_upper_limit_ <= 0)
    {
        log_warn("Warning (in Integral/Integral/1): maxSteps = %i must be > 0, setting to 1", maxSteps_upper_limit_);
        maxSteps_upper_limit_ = 1;
    }

    if (precision_ <= 0)
    {
        log_warn("Warning (in Integral/Integral/2): precision = %f must be > 0, setting to 1.e-6", precision_);
        precision_ = 1.e-6;
    }

    iX_.resize(maxSteps_upper_limit_);
    iY_.resize(maxSteps_upper_limit_);

    aux = std::max(romberg_, romberg4refine_);
    c_.resize(aux);
    d_.resize(aux);

    q_last_3_results_.resize(3);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Integral::Integral(const Integral& integral)
    : maxSteps_romberg_(integral.maxSteps_romberg_)
    , maxSteps_upper_limit_(integral.maxSteps_upper_limit_)
    , romberg_(integral.romberg_)
    , precision_(integral.precision_)
    , max_(integral.max_)
    , min_(integral.min_)
    , iX_(integral.iX_)
    , iY_(integral.iY_)
    , c_(integral.c_)
    , d_(integral.d_)
    , romberg4refine_(integral.romberg4refine_)
    , powerOfSubstitution_(integral.powerOfSubstitution_)
    , randomDo_(integral.randomDo_)
    , useLog_(integral.useLog_)
    , randomNumber_(integral.randomNumber_)
    , randomX_(integral.randomX_)
    , reverse_(integral.reverse_)
    , reverseX_(integral.reverseX_)
    , savedResult_(integral.savedResult_)
    , q_last_3_results_(integral.q_last_3_results_)
    , q_rlist2_(integral.q_rlist2_)
    , q_iord_(integral.q_iord_)
{
    integrand_ = std::ref(integral.integrand_);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Integral::Integral(int romberg, int maxSteps, double precision)
    : maxSteps_romberg_(12)
    , max_(1)
    , min_(0)
    , iX_()
    , iY_()
    , c_()
    , d_()
    , romberg4refine_(2)
    , powerOfSubstitution_(0)
    , randomDo_(false)
    , useLog_(false)
    , randomNumber_(0)
    , randomX_(0)
    , reverse_(false)
    , reverseX_(0)
    , savedResult_(0)
    , q_last_3_results_()
    , q_rlist2_()
    , q_iord_()
{
    int aux;
    if (romberg <= 0)
    {
        log_warn("Warning (in Integral/Integral/0): romberg = %i must be > 0, setting to 1", romberg);
        romberg = 1;
    }

    if (maxSteps <= 0)
    {
        log_warn("Warning (in Integral/Integral/1): maxSteps = %i must be > 0, setting to 1", maxSteps);
        maxSteps = 1;
    }

    if (precision <= 0)
    {
        log_warn("Warning (in Integral/Integral/2): precision = %f must be > 0, setting to 1.e-6", precision);
        precision = 1.e-6;
    }

    this->romberg_              = romberg;
    this->maxSteps_upper_limit_ = maxSteps;
    this->precision_            = precision;

    iX_.resize(maxSteps_upper_limit_);
    iY_.resize(maxSteps_upper_limit_);

    aux = std::max(romberg_, romberg4refine_);
    c_.resize(aux);
    d_.resize(aux);

    q_last_3_results_.resize(3);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Integral& Integral::operator=(const Integral& integral)
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

bool Integral::operator==(const Integral& integral) const
{
    // if(integrand_ != integral.integrand_)     return false;
    if (iX_.size() != integral.iX_.size())
        return false;
    if (iY_.size() != integral.iY_.size())
        return false;
    if (c_.size() != integral.c_.size())
        return false;
    if (d_.size() != integral.d_.size())
        return false;

    for (unsigned int i = 0; i < iX_.size(); i++)
    {
        if (iX_.at(i) != integral.iX_.at(i))
            return false;
    }
    for (unsigned int i = 0; i < iY_.size(); i++)
    {
        if (iY_.at(i) != integral.iY_.at(i))
            return false;
    }
    for (unsigned int i = 0; i < c_.size(); i++)
    {
        if (c_.at(i) != integral.c_.at(i))
            return false;
    }
    for (unsigned int i = 0; i < d_.size(); i++)
    {
        if (d_.at(i) != integral.d_.at(i))
            return false;
    }
    if (maxSteps_upper_limit_ != integral.maxSteps_upper_limit_)
        return false;
    if (maxSteps_romberg_ != integral.maxSteps_romberg_)
        return false;
    if (romberg_ != integral.romberg_)
        return false;
    if (precision_ != integral.precision_)
        return false;
    if (max_ != integral.max_)
        return false;
    if (min_ != integral.min_)
        return false;
    if (romberg4refine_ != integral.romberg4refine_)
        return false;
    if (randomDo_ != integral.randomDo_)
        return false;
    if (useLog_ != integral.useLog_)
        return false;
    if (randomNumber_ != integral.randomNumber_)
        return false;
    if (randomX_ != integral.randomX_)
        return false;
    if (reverse_ != integral.reverse_)
        return false;
    if (reverseX_ != integral.reverseX_)
        return false;
    if (savedResult_ != integral.savedResult_)
        return false;

    if (powerOfSubstitution_ != integral.powerOfSubstitution_)
        return false;

    else
        return true;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool Integral::operator!=(const Integral& integral) const
{
    return !(*this == integral);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------/

void Integral::swap(Integral& integral)
{
    using std::swap;

    swap(maxSteps_upper_limit_, integral.maxSteps_upper_limit_);
    swap(maxSteps_romberg_, integral.maxSteps_romberg_);
    swap(romberg_, integral.romberg_);
    swap(precision_, integral.precision_);
    swap(max_, integral.max_);
    swap(min_, integral.min_);

    iX_.swap(integral.iX_);
    iY_.swap(integral.iY_);

    c_.swap(integral.c_);
    d_.swap(integral.d_);

    integrand_ = std::ref(integral.integrand_);

    swap(romberg4refine_, integral.romberg4refine_);
    swap(powerOfSubstitution_, integral.powerOfSubstitution_);
    swap(randomDo_, integral.randomDo_);
    swap(useLog_, integral.useLog_);
    swap(randomNumber_, integral.randomNumber_);
    swap(randomX_, integral.randomX_);
    swap(reverse_, integral.reverse_);
    swap(reverseX_, integral.reverseX_);
    swap(savedResult_, integral.savedResult_);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::Function(double x)
{
    double result, t;
    if (reverse_)
    {
        x = reverseX_ - x;
    }

    if (powerOfSubstitution_ == 0)
    {
        t      = x;
        result = 1;
    } else if (powerOfSubstitution_ > 0)
    {
        t      = std::pow(x, -powerOfSubstitution_);
        result = powerOfSubstitution_ * (t / x);
    } else
    {
        t      = -std::pow(-x, powerOfSubstitution_);
        result = -powerOfSubstitution_ * (t / x);
    }

    if (useLog_)
    {
        t = std::exp(t);
        result *= t;
    }
    result *= integrand_(t);

    if (result != result)
    {
        if (integrand_(t) == 0)
        {
            log_info("substitution not suitable! returning 0!");
            return 0;
        } else
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

    if (n == 1)
    {
        return (Function(max_) + Function(min_)) * (max_ - min_) / 2;
    }

    n /= 2;
    stepSize  = (max_ - min_) / n;
    resultSum = 0;

    for (xStep = min_ + stepSize / 2; xStep < max_; xStep += stepSize)
    {
        resultSum += Function(xStep);
    }

    return (oldSum + resultSum * stepSize) / 2;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::Trapezoid3(int n, double oldSum)
{
    double xStep, stepSize, resultSum;

    if (n == 1)
    {
        return (max_ - min_) * Function((max_ + min_) / 2);
    }

    stepSize  = (max_ - min_) / n;
    resultSum = 0;

    for (xStep = min_ + stepSize / 2; xStep < max_; xStep += stepSize)
    {
        resultSum += Function(xStep);
        xStep += 2 * stepSize;
        resultSum += Function(xStep);
    }

    return oldSum / 3 + resultSum * stepSize;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::Trapezoid3S(int n, double oldSum, int stepNumber)
{
    double xStep, stepSize, resultSum;
    double smallSum = 0, approX = 0, functionValue1, functionValue2, sumDifference;
    double functionDifference, functionSum, aEq, bEq, bEq2, determinant;
    bool flag;

    if (n == 1)
    {
        return (max_ - min_) * Function((max_ + min_) / 2);
    }

    stepSize = (max_ - min_) / n;
    if (stepNumber >= romberg_ - 1)
    {
        if (randomNumber_ >= 0)
        {
            smallSum = randomNumber_ * oldSum / (1.5 * stepSize);
        } else
        {
            smallSum = -randomNumber_ * oldSum / (1.5 * stepSize);
            if (oldSum < 0)
            {
                smallSum *= -1;
            }
        }
    }
    resultSum = 0;
    flag      = false;
    for (xStep = min_ + stepSize / 2; xStep < max_; xStep += stepSize)
    {
        resultSum += functionValue1 = Function(xStep);
        xStep += 2 * stepSize;
        resultSum += functionValue2 = Function(xStep);

        if (!flag)
            if (stepNumber >= romberg_ - 1)
                if (std::abs(resultSum) >= std::abs(smallSum))
                {
                    functionSum = functionValue1 + functionValue2;

                    sumDifference = (smallSum - (resultSum - functionSum));
                    sumDifference *= 1.5 * stepSize;

                    functionDifference = functionValue2 - functionValue1;

                    aEq  = functionDifference / stepSize;
                    bEq  = (functionValue2 - 5 * functionValue1) / 2;
                    bEq2 = bEq * bEq;

                    if (std::abs(aEq * sumDifference) < precision_ * bEq2)
                    {
                        approX = sumDifference * 2 / functionSum;
                    } else
                    {
                        determinant = bEq2 + 4 * aEq * sumDifference;
                        if (determinant >= 0)
                        {
                            determinant = std::sqrt(determinant);
                            approX      = (bEq + determinant) / aEq;

                            if (approX < 0 || approX > 3 * stepSize)
                            {
                                approX = (bEq - determinant) / aEq;
                            } else if (approX < 0 || approX > 3 * stepSize)
                            {
                                approX = sumDifference * 2 / functionSum;
                            }
                        } else
                        {
                            log_warn("Warning (in Integral/trapezoid3S): Determinant is negative, proceeding with "
                                     "linear approximation");
                            approX = sumDifference * 2 / functionSum;
                        }
                    }

                    approX += xStep - 2.5 * stepSize;
                    flag = true;
                }
    }

    if (stepNumber >= romberg_ - 1)
    {
        if (!flag)
        {
            // Set randomX_ to max_, but not exactly to max_ to
            // avoid differegences in Ionization::InelCorrection.
            // Use the next lower representation.
#if __cplusplus > 199711L
            approX = std::nexttoward(max_, min_);
#else
            approX = max_ - std::abs(min_) * COMPUTER_PRECISION;
#endif
        }
        randomX_ = approX;
    }
    return oldSum / 3 + resultSum * stepSize;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Integral::InterpolationResults Integral::Interpolate(int start, double x)
{
    int num, i, k;
    bool dd;
    double error = 0, result = 0;
    double aux, aux2, dx1, dx2;
    Integral::InterpolationResults interpolation_results;

    num = 0;
    aux = std::abs(x - iX_[start + 0]);

    for (i = 0; i < romberg_; i++)
    {
        aux2 = std::abs(x - iX_[start + i]);
        if (aux2 < aux)
        {
            num = i;
            aux = aux2;
        }
        c_[i] = iY_[start + i];
        d_[i] = iY_[start + i];
    }

    if (num == 0)
    {
        dd = true;
    } else if (num == romberg_ - 1)
    {
        dd = false;
    } else
    {
        if (std::abs(x - iX_[start + num - 1]) > std::abs(x - iX_[start + num + 1]))
        {
            dd = true;
        } else
        {
            dd = false;
        }
    }

    result = iY_[start + num];
    for (k = 1; k < romberg_; k++)
    {
        for (i = 0; i < romberg_ - k; i++)
        {
            dx1 = iX_[start + i] - x;
            dx2 = iX_[start + i + k] - x;

            aux  = c_[i + 1] - d_[i];
            aux2 = dx1 - dx2;

            if (aux2 != 0)
            {
                aux   = aux / aux2;
                c_[i] = dx1 * aux;
                d_[i] = dx2 * aux;
            } else
            {
                c_[i] = 0;
                d_[i] = 0;
            }
        }

        if (num == 0)
        {
            dd = true;
        }

        if (num == romberg_ - k)
        {
            dd = false;
        }

        if (dd)
        {
            error = c_[num];
        } else
        {
            num--;
            error = d_[num];
        }

        dd = !dd;
        result += error;
    }

    interpolation_results.Error = error;
    interpolation_results.Value = result;
    return interpolation_results;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::RombergIntegrateClosed()
{
    int k = 1;
    double n = 1;
    double error, result, value;
    Integral::InterpolationResults interpolation_results;

    value  = 0;
    result = 0;

    for (int i = 0; i < maxSteps_romberg_; i++)
    {
        result = Trapezoid(k, result);
        iX_[i] = n;
        iY_[i] = result;

        if (i >= romberg_ - 1)
        {
            interpolation_results = Interpolate(i - (romberg_ - 1), 0);

            error = interpolation_results.Error;
            value = interpolation_results.Value;

            if (value != 0)
            {
                error /= value;
            }

            if (std::abs(error) < precision_)
            {
                return value;
            }
        }

        k = k * 2;
        n = n / 4;
    }

    QuadpackResults q_results = qags();
    log_warn("Precision %e has not been reached after %i steps the value is %e!\nUsing now qags! value = %e, abserr = "
             "%e, neval = %i, ier = %i",
             precision_,
             maxSteps_romberg_,
             value,
             q_results.value,
             q_results.abserr,
             q_results.neval,
             q_results.ier);
    return q_results.value;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::RombergIntegrateOpened()
{
    int i, k;
    double n;
    double error, result, value;
    Integral::InterpolationResults interpolation_results;

    k      = 1;
    n      = 1;
    value  = 0;
    result = 0;

    for (i = 0; i < maxSteps_romberg_; i++)
    {
        if (randomNumber_ == 0 || randomNumber_ == 1)
        {
            result = Trapezoid3(k, result);
            if (result != result)
            {
                log_error("Function Value of is nan! Returning 0!");
                return 0;
            }

            if (randomNumber_ == 0)
            {
                randomX_ = min_;
            } else
            {
                randomX_ = max_;
            }
        } else
        {
            result = Trapezoid3S(k, result, i);
        }

        iX_[i] = n;
        iY_[i] = result;
        if (i >= romberg_ - 1)
        {
            interpolation_results = Interpolate(i - (romberg_ - 1), 0);
            error                 = interpolation_results.Error;
            value                 = interpolation_results.Value;

            if (value != 0)
            {
                error /= value;
            }

            if (std::abs(error) < precision_)
            {
                return value;
            }
        }

        k *= 3;
        n /= 9;
    }

    QuadpackResults q_results = qags();
    log_warn("Precision %e has not been reached after %i steps the value is %e!\nUsing now qags! value = %e, abserr = "
             "%e, neval = %i, ier = %i",
             precision_,
             maxSteps_romberg_,
             value,
             q_results.value,
             q_results.abserr,
             q_results.neval,
             q_results.ier);
    return q_results.value;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::RombergIntegrateOpened(double bigValue)
{
    int i, k;
    double n;
    double error, result, value;
    InterpolationResults interpolation_results;

    k      = 1;
    n      = 1;
    value  = 0;
    result = 0;

    for (i = 0; i < maxSteps_romberg_; i++)
    {
        result = Trapezoid3(k, result);
        iX_[i] = n;
        iY_[i] = result;
        if (i >= romberg_ - 1)
        {
            interpolation_results = Interpolate(i - (romberg_ - 1), 0);
            error                 = interpolation_results.Error;
            value                 = interpolation_results.Value;
            error /= bigValue;

            if (std::abs(error) < precision_)
            {
                return value;
            }
        }
        k *= 3;
        n /= 9;
    }

    QuadpackResults q_results = qags();
    log_warn("Precision %e has not been reached after %i steps the value is %e!\nUsing now qags! value = %e, abserr = "
             "%e, neval = %i, ier = %i",
             precision_,
             maxSteps_romberg_,
             value,
             q_results.value,
             q_results.abserr,
             q_results.neval,
             q_results.ier);
    return q_results.value;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::InitIntegralOpenedAndClosed(double min, double max, std::function<double(double)> integrand)
{
    double aux;

    reverse_ = false;
    useLog_  = false;

    if (min > max)
    {
        std::swap(min, max);
        aux = -1;
    } else
    {
        aux = 1;
    }

    this->min_           = min;
    this->max_           = max;
    integrand_           = integrand;
    powerOfSubstitution_ = 0;
    randomDo_            = false;

    return aux;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::IntegrateClosed(double min, double max, std::function<double(double)> integrand)
{
    double aux;
    aux = InitIntegralOpenedAndClosed(min, max, integrand);

    if (std::abs(max_ - min_) <= std::abs(min_) * COMPUTER_PRECISION)
    {
        return 0;
    }
    return aux * RombergIntegrateClosed();
}

//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::IntegrateOpened(double min, double max, std::function<double(double)> integrand)
{
    double aux;
    aux = InitIntegralOpenedAndClosed(min, max, integrand);

    if (std::abs(max_ - min_) <= std::abs(min_) * COMPUTER_PRECISION)
    {
        return 0;
    }

    return aux * RombergIntegrateOpened();
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::InitIntegralWithSubstitution(double min,
                                              double max,
                                              std::function<double(double)> integrand,
                                              double powerOfSubstitution)
{
    double aux;

    reverse_ = false;
    useLog_  = false;

    if (min > max)
    {
        std::swap(min, max);
        aux = -1;
    } else
    {
        aux = 1;
    }

    if (powerOfSubstitution > 0)
    {
        if (max > 0 && min > 0)
        {
            this->min_ = std::pow(max, -1 / powerOfSubstitution);
            this->max_ = std::pow(min, -1 / powerOfSubstitution);
        } else if (max > 0)
        {
            this->min_ = 0;
            this->max_ = std::pow(max, -1 / powerOfSubstitution);
            aux        = -aux;
        } else
        {
            return 0;
        }
    } else if (powerOfSubstitution < 0)
    {
        if (max < 0 && min < 0)
        {
            this->min_ = -std::pow(-max, 1 / powerOfSubstitution);
            this->max_ = -std::pow(-min, 1 / powerOfSubstitution);
        } else if (min < 0)
        {
            this->min_ = -std::pow(-min, 1 / powerOfSubstitution);
            this->max_ = 0;
            aux        = -aux;
        } else
        {
            return 0;
        }
    } else
    {
        this->min_ = min;
        this->max_ = max;
    }

    this->integrand_           = integrand;
    this->powerOfSubstitution_ = powerOfSubstitution;
    randomNumber_              = 0;
    randomDo_                  = false;

    return aux;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::IntegrateWithSubstitution(double min,
                                           double max,
                                           std::function<double(double)> integrand,
                                           double powerOfSubstitution)
{
    double aux;

    aux = InitIntegralWithSubstitution(min, max, integrand, powerOfSubstitution);

    if (std::abs(max_ - min_) <= std::abs(min_) * COMPUTER_PRECISION)
    {
        return 0;
    }

    return aux * RombergIntegrateOpened();
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Integral::RefineUpperLimit(double result)
{
    int i, rombergStore;
    double maxStore, minStore, functionValue, f, df;
    double deltaX, deltaOld, currentX, aux;
    double xlow, xhi, flow, fhi;

    // This was commented out but is necessary since
    // it ensures to find a value if the searched integral value
    // is out of integral range.
    if (randomNumber_ == 0 || randomNumber_ == 1)
    {
        return;
    }

    xlow = min_;
    xhi  = max_;
    flow = -randomNumber_ * result;
    fhi  = (1 - randomNumber_) * result;

    if (flow * fhi > 0)
    {
        log_error("Root must be bracketed");
        return;
    }

    if (flow > 0)
    {
        std::swap(xlow, xhi);
        // std::swap(flow, fhi);
    }

    deltaX   = max_ - min_;
    deltaOld = max_ - min_;

    if (randomX_ < min_ || randomX_ > max_)
    {
        randomX_ = (min_ + max_) / 2;
    }
    currentX     = randomX_;
    rombergStore = romberg_;
    minStore     = min_;
    maxStore     = max_;
    max_         = randomX_;

    functionValue = RombergIntegrateOpened(result) - randomNumber_ * result;

    romberg_ = romberg4refine_;
    f        = functionValue;
    df       = Function(currentX);

    for (i = 0; i < maxSteps_upper_limit_; i++)
    {
        if (f < 0)
        {
            xlow = currentX;
            // flow = f;
        } else
        {
            xhi = currentX;
            // fhi = f;
        }

        if (((currentX - xhi) * df - f) * ((currentX - xlow) * df - f) > 0 || std::abs(2 * f) > std::abs(deltaOld * df))
        {
            deltaOld = deltaX;
            deltaX   = (xhi - xlow) / 2;
            currentX = xlow + deltaX;

            if (xlow == currentX)
            {
                break;
            }
        } else
        {
            deltaOld = deltaX;
            deltaX   = f / df;
            aux      = currentX;
            currentX -= deltaX;

            if (aux == currentX)
            {
                break;
            }
        }

        if (std::abs(deltaX) < precision_)
        {
            break;
        }

        min_ = randomX_;
        max_ = currentX;

        if (min_ > max_)
        {
            std::swap(min_, max_);
            aux = -1;
        } else
        {
            aux = 1;
        }

        f  = functionValue + aux * RombergIntegrateOpened(result);
        df = Function(currentX);
    }

    if (i == maxSteps_upper_limit_)
    {
        log_warn("Precision %f has not been reached after %i steps!", precision_, maxSteps_upper_limit_);
    }

    randomX_ = currentX;
    romberg_ = rombergStore;
    min_     = minStore;
    max_     = maxStore;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::InitIntegralWithLog(double min, double max, std::function<double(double)> integrand)
{
    double aux;

    reverse_ = false;
    useLog_  = true;

    if (min > max)
    {
        std::swap(min, max);
        aux = -1;

    } else
    {
        aux = 1;
    }

    this->min_           = std::log(min);
    this->max_           = std::log(max);
    this->integrand_     = integrand;
    powerOfSubstitution_ = 0;
    randomNumber_        = 0;
    randomDo_            = false;

    return aux;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::IntegrateWithLog(double min, double max, std::function<double(double)> integrand)
{
    double aux;

    aux = InitIntegralWithLog(min, max, integrand);

    if (std::abs(max_ - min_) <= std::abs(min_) * COMPUTER_PRECISION)
    {
        return 0;
    }

    return aux * RombergIntegrateOpened();
}

//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::InitIntegralWithLogSubstitution(double min,
                                                 double max,
                                                 std::function<double(double)> integrand,
                                                 double powerOfSubstitution)
{
    double aux;

    reverse_ = false;
    useLog_  = true;

    if (min < 0 || max < 0)
    {
        return 0;
    } else if (min > max)
    {
        std::swap(min, max);
        aux = -1;
    } else
    {
        aux = 1;
    }

    if (powerOfSubstitution > 0)
    {
        if (max > 1 && min > 1)
        {
            this->min_ = std::pow(std::log(max), -1 / powerOfSubstitution);
            this->max_ = std::pow(std::log(min), -1 / powerOfSubstitution);
        } else if (max > 1)
        {
            this->min_ = 0;
            this->max_ = std::pow(std::log(max), -1 / powerOfSubstitution);
            aux        = -aux;
        } else
        {
            return 0;
        }
    } else if (powerOfSubstitution < 0)
    {
        if (max < 1 && min < 1)
        {
            this->min_ = -std::pow(-std::log(max), 1 / powerOfSubstitution);
            this->max_ = -std::pow(-std::log(min), 1 / powerOfSubstitution);
        } else if (min < 1)
        {
            this->min_ = -std::pow(-std::log(min), 1 / powerOfSubstitution);
            this->max_ = 0;
            aux        = -aux;
        } else
        {
            return 0;
        }
    } else
    {
        this->min_ = std::log(min);
        this->max_ = std::log(max);
    }

    this->integrand_           = integrand;
    this->powerOfSubstitution_ = powerOfSubstitution;
    randomNumber_              = 0;
    randomDo_                  = false;

    return aux;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Integral::IntegrateWithLogSubstitution(double min,
                                              double max,
                                              std::function<double(double)> integrand,
                                              double powerOfSubstitution)
{

    double aux;

    aux = InitIntegralWithLogSubstitution(min, max, integrand, powerOfSubstitution);

    if (std::abs(max_ - min_) <= std::abs(min_) * COMPUTER_PRECISION)
    {
        return 0;
    }

    return aux * RombergIntegrateOpened();
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Integral::QuadpackResults Integral::qags(double q_limit, double q_epsabs, double q_epsrel)
{

    static const double q_epmach = std::numeric_limits<double>::epsilon(); // machine epsilon
    static const double q_uflow  = std::numeric_limits<double>::min();     // smallest finite value
    static const double q_oflow  = std::numeric_limits<double>::max();     // largest finite value
    // output parameter
    QuadpackResults output;
    // double result = 0.0;
    // double abserr = 0.0;
    // int neval = 0;
    // int ier = 0;

    if (q_limit < 1)
    {
        log_warn("the limit is below 1, returning 0");
        return output;
    }

    if (q_epsabs < 0. && q_epsrel < 0.)
    {
        output.ier = 6;
        return output;
    }
    std::pair<Integral::InterpolationResults, Integral::InterpolationResults> qk21_output;

    // First approximation to the integral.
    qk21_output   = q_gaus_kronrod_21(min_, max_);
    output.value  = qk21_output.first.Value;
    output.abserr = qk21_output.first.Error;
    double defabs = qk21_output.second.Value;
    double resabs = qk21_output.second.Error;

    // Test on accuracy.
    double dres   = std::abs(output.value);
    double errbnd = std::max(q_epsabs, q_epsrel * dres);
    int last      = 1;

    if (output.abserr <= 100 * q_epmach * defabs && output.abserr > errbnd)
        output.ier = 2;

    if (q_limit == 1)
        output.ier = 1;

    if (output.ier != 0 || (output.abserr <= errbnd && output.abserr != resabs) || output.abserr == 0.0)
    {
        output.neval = 42 * last - 21;
        return output;
    }
    std::vector<double> q_alist_;
    std::vector<double> q_blist_;
    std::vector<double> q_elist_;
    std::vector<double> q_rlist_;

    // the maximum number of elements the epsilon table can contain.
    // if this number is reached, the upper diagonal of the epsilon table is deleted.
    const int q_limit_epsilon_table_ = 50;
    q_rlist2_.resize(q_limit_epsilon_table_ + 2);
    q_alist_.resize(q_limit);
    q_blist_.resize(q_limit);
    q_elist_.resize(q_limit);
    q_rlist_.resize(q_limit);
    q_iord_.resize(q_limit);

    int ierro   = 0;
    q_alist_[0] = min_;
    q_blist_[0] = max_;
    q_elist_[0] = output.abserr;
    q_rlist_[0] = output.value;
    q_iord_[0]  = 1;

    // local parameters
    double abseps, area1, area12, area2, a1, a2, b1, b2, defab1, defab2, erlarg, erlast, error1, error2, erro12, reseps;
    int jupbnd;
    Integral::InterpolationResults extrapolation_output;

    // Initialization
    q_rlist2_[0]  = output.value;
    double area   = output.value;
    double errmax = output.abserr;
    double errsum = output.abserr;
    int maxerr    = 1;
    int nrmax     = 1;
    int nres      = 0;
    int numrl2    = 2;
    int ktmin     = 0;
    bool extrap   = false;
    bool noext    = false;
    int iroff1    = 0;
    int iroff2    = 0;
    int iroff3    = 0;

    double small  = std::abs(max_ - min_) * 0.375;
    double ertest = 0.0; // Is set in loop when last == 2
    double correc = 0.0; // Is set in loop if error is not reached

    int ksgn = -1;
    if (dres >= (1.0 - 50 * q_epmach) * defabs)
        ksgn = 1;

    output.abserr = q_oflow;

    // main loop
    for (last = 2; last <= q_limit; last++)
    {
        // Bisect the subinterval with the nrmax-th largest error estimate.
        a1     = q_alist_[maxerr - 1];
        b1     = 0.5 * (q_alist_[maxerr - 1] + q_blist_[maxerr - 1]);
        a2     = b1;
        b2     = q_blist_[maxerr - 1];
        erlast = errmax;

        qk21_output = q_gaus_kronrod_21(a1, b1);
        area1       = qk21_output.first.Value;
        error1      = qk21_output.first.Error;
        resabs      = qk21_output.second.Value;
        defab1      = qk21_output.second.Error;
        qk21_output = q_gaus_kronrod_21(a2, b2);
        area2       = qk21_output.first.Value;
        error2      = qk21_output.first.Error;
        resabs      = qk21_output.second.Value;
        defab2      = qk21_output.second.Error;

        // Improve previous approximations to integral and error and test for accuracy.
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area + area12 - q_rlist_[maxerr - 1];

        if (defab1 != error1 && defab2 != error2)
        {
            if (std::abs(q_rlist_[maxerr - 1] - area12) <= 1.0e-5 * std::abs(area12) && erro12 >= 0.99 * errmax)
            {
                if (extrap)
                {
                    iroff2 = iroff2 + 1;
                } else
                {
                    iroff1 = iroff1 + 1;
                }
            }
            if (last > 10 && erro12 > errmax)
                iroff3 = iroff3 + 1;
        }

        q_rlist_[maxerr - 1] = area1;
        q_rlist_[last - 1]   = area2;
        errbnd               = std::max(q_epsabs, q_epsrel * std::abs(area));

        // Test for roundoff error and eventually set error flag.

        if (iroff1 + iroff2 >= 10 || iroff3 >= 20)
            output.ier = 2;

        if (iroff2 >= 5)
            ierro = 3;

        // Set error flag in the case that the number of subintervals equals limit.
        if (last == q_limit)
            output.ier = 1;

        // Set error flag in the case of bad integrand behavior at a point of the integration range.
        // if(std::max(std::abs(a1), std::abs(b2)) <= (1 + 1000*q_epmach) * (std::abs(a2) + 1000*q_uflow)) // diff
        // between f90 & f77
        if (std::max(std::abs(a1), std::abs(b2)) <= (1 + 100 * q_epmach) * (std::abs(a2) + 1000 * q_uflow))
            output.ier = 4;

        // Append the newly-created intervals to the list.
        if (error2 <= error1)
        {
            q_alist_[last - 1]   = a2;
            q_blist_[maxerr - 1] = b1;
            q_blist_[last - 1]   = b2;
            q_elist_[maxerr - 1] = error1;
            q_elist_[last - 1]   = error2;
        } else
        {
            q_alist_[maxerr - 1] = a2;
            q_alist_[last - 1]   = a1;
            q_blist_[last - 1]   = b1;
            q_rlist_[maxerr - 1] = area2;
            q_rlist_[last - 1]   = area1;
            q_elist_[maxerr - 1] = error2;
            q_elist_[last - 1]   = error1;
        }

        // Call QSORT to maintain the descending ordering in the list of error estimates
        // and select the subinterval with nrmax-th largest error estimate (to be bisected next).
        nrmax  = q_sort(q_limit, last, maxerr, nrmax, q_elist_);
        maxerr = q_iord_[nrmax - 1];
        errmax = q_elist_[maxerr - 1];
        if (errsum <= errbnd)
        {
            // Compute global integral sum.
            output.value  = std::accumulate(q_rlist_.begin(), q_rlist_.begin() + (last), 0.0);
            output.abserr = errsum;

            if (2 < output.ier)
                output.ier -= 1;

            output.neval = 42 * last - 21;
            return output;
        }

        if (output.ier != 0)
            break;

        if (last == 2)
        {
            erlarg       = errsum;
            ertest       = errbnd;
            q_rlist2_[1] = area;
            continue;
        }

        if (noext)
            continue;

        erlarg = erlarg - erlast;

        if (std::abs(b1 - a1) > small)
            erlarg = erlarg + erro12;

        // Test whether the interval to be bisected next is the smallest interval.
        if (!extrap)
        {
            if (std::abs(q_blist_[maxerr - 1] - q_alist_[maxerr - 1]) > small)
                continue;

            extrap = true;
            nrmax  = 2;
        }

        // The smallest interval has the largest error.
        // Before bisecting decrease the sum of the errors over
        // the larger intervals (erlarg) and perform extrapolation.
        if (ierro != 3 && erlarg > ertest)
        {
            if (last > (2 + q_limit / 2))
            {
                jupbnd = q_limit + 3 - last;
            } else
            {
                jupbnd = last;
            }

            bool go_to_main_loop_start = false;
            int id                     = nrmax;
            for (int idx = id; idx <= jupbnd; idx++)
            {
                maxerr = q_iord_[nrmax - 1];
                errmax = q_elist_[maxerr - 1];
                if (std::abs(q_blist_[maxerr - 1] - q_alist_[maxerr - 1]) > small)
                {
                    go_to_main_loop_start = true;
                    break;
                }
                nrmax = nrmax + 1;
            }
            // while(nrmax <= jupbnd)
            // {
            //     maxerr = q_iord_[nrmax-1];
            //     errmax = q_elist_[maxerr-1];
            //     if(std::abs(q_blist_[maxerr-1] - q_alist_[maxerr-1]) > small)
            //     {
            //         go_to_main_loop_start = true;
            //         break;
            //     }
            //     nrmax = nrmax+1;
            // }
            if (go_to_main_loop_start)
                continue;
        }

        // Perform extrapolation.
        numrl2                = numrl2 + 1;
        q_rlist2_[numrl2 - 1] = area;
        nres += 1;
        extrapolation_output = q_epsilon_extrapolation(q_limit_epsilon_table_, numrl2, nres); // f90 q_rlist2_ -> epstab
        reseps               = extrapolation_output.Value;
        abseps               = extrapolation_output.Error;
        ktmin                = ktmin + 1;

        if (ktmin > 5 && output.abserr < 1.0e-3 * errsum)
            output.ier = 5;

        if (abseps < output.abserr)
        {
            ktmin         = 0;
            output.abserr = abseps;
            output.value  = reseps;
            correc        = erlarg;
            ertest        = std::max(q_epsabs, q_epsrel * std::abs(reseps));

            if (output.abserr <= ertest)
                break;
        }

        // Prepare bisection of the smallest interval.
        if (numrl2 == 1)
            noext = true;

        if (output.ier == 5)
            break;

        maxerr = q_iord_[0];
        errmax = q_elist_[maxerr - 1];
        nrmax  = 1;
        extrap = false;
        small  = small * 0.5;
        erlarg = errsum;
    } // end of main loop

    // Set final result and error estimate.
    if (output.abserr == q_oflow)
    {
        // Compute global integral sum.
        output.value  = std::accumulate(q_rlist_.begin(), q_rlist_.begin() + (last), 0.0);
        output.abserr = errsum;
    } else
    {
        // Test on divergence.
        if (output.ier + ierro == 0)
        {
            if (ksgn != -1 || std::max(std::abs(output.value), std::abs(area)) > defabs * 0.01)
            {
                if (0.01 > (output.value / area) || (output.value / area) > 100 || errsum > std::abs(area))
                    output.ier = 6;
            }
        }

        if (ierro == 3)
            output.abserr += correc;

        if (output.ier == 0)
            output.ier = 3;

        if (output.value != 0.0 && area != 0.0)
        {
            if (output.abserr / std::abs(output.value) > errsum / std::abs(area))
            {
                // // Compute global integral sum.
                output.value  = std::accumulate(q_rlist_.begin(), q_rlist_.begin() + (last), 0.0);
                output.abserr = errsum;
            }
        } else
        {
            if (output.abserr > errsum)
            {
                // Compute global integral sum.
                output.value  = std::accumulate(q_rlist_.begin(), q_rlist_.begin() + (last), 0.0);
                output.abserr = errsum;
            } else
            {
                if (area != 0.0)
                {
                    if (ksgn != -1 || std::max(std::abs(output.value), std::abs(area)) > defabs * 0.01)
                    {
                        if (0.01 > (output.value / area) || (output.value / area) > 100 || errsum > std::abs(area))
                            output.ier = 6;
                    }
                }
            }
        }
    }

    if (2 < output.ier)
        output.ier -= 1;

    output.neval = 42 * last - 21;
    return output;
}

std::pair<Integral::InterpolationResults, Integral::InterpolationResults> Integral::q_gaus_kronrod_21(double min_lim,
                                                                                                      double max_lim)
{
    double q_fv1[10];
    double q_fv2[10];
    const double q_epmach                  = std::numeric_limits<double>::epsilon(); // machine epsilon
    const double q_uflow                   = std::numeric_limits<double>::min();     // smallest finite value
    const double q_weights_gaus_10p[5]     = { 0.066671344308688137593568809893332,
                                           0.149451349150580593145776339657697,
                                           0.219086362515982043995534934228163,
                                           0.269266719309996355091226921569469,
                                           0.295524224714752870173892994651338 };
    const double q_weights_kronrod_21p[11] = { 0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
                                               0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
                                               0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
                                               0.123491976262065851077958109831074, 0.134709217311473325928054001771707,
                                               0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
                                               0.149445554002916905664936468389821 };
    const double q_abscissae_kronrod_21p[11] = {
        0.995657163025808080735527280689003, 0.973906528517171720077964012084452, 0.930157491355708226001207180059508,
        0.865063366688984510732096688423493, 0.780817726586416897063717578345042, 0.679409568299024406234327365114874,
        0.562757134668604683339000099272694, 0.433395394129247190799265943165784, 0.294392862701460198131126603103866,
        0.148874338981631210884826001129720, 0.000000000000000000000000000000000
    };
    // output parameters
    double result, abserr, resabs, resasc;
    Integral::InterpolationResults qk21_output;
    Integral::InterpolationResults qk21_abs_output;
    // internal parameters
    double absc, fval1, fval2, reskh;

    double centr  = 0.5 * (min_lim + max_lim);
    double hlgth  = 0.5 * (max_lim - min_lim);
    double dhlgth = std::abs(hlgth);

    // Compute the 21-point Kronrod approximation to the
    // integral, and estimate the absolute error.

    double resg = 0.0;
    double fc   = Function(centr);
    double resk = q_weights_kronrod_21p[10] * fc;
    resabs      = std::abs(resk);

    for (int idx = 1; idx <= 9; idx = idx + 2)
    {
        absc       = hlgth * q_abscissae_kronrod_21p[idx];
        fval1      = Function(centr - absc);
        fval2      = Function(centr + absc);
        q_fv1[idx] = fval1;
        q_fv2[idx] = fval2;
        resg       = resg + q_weights_gaus_10p[(idx - 1) / 2] * (fval1 + fval2);
        resk       = resk + q_weights_kronrod_21p[idx] * (fval1 + fval2);
        resabs     = resabs + q_weights_kronrod_21p[idx] * (std::abs(fval1) + std::abs(fval2));
    }

    for (int idx = 0; idx <= 8; idx = idx + 2)
    {
        absc       = hlgth * q_abscissae_kronrod_21p[idx];
        fval1      = Function(centr - absc);
        fval2      = Function(centr + absc);
        q_fv1[idx] = fval1;
        q_fv2[idx] = fval2;
        resk       = resk + q_weights_kronrod_21p[idx] * (fval1 + fval2);
        resabs     = resabs + q_weights_kronrod_21p[idx] * (std::abs(fval1) + std::abs(fval2));
    }

    reskh  = resk * 0.5;
    resasc = q_weights_kronrod_21p[10] * std::abs(fc - reskh);

    for (int idx = 0; idx < 10; idx++)
    {
        resasc = resasc + q_weights_kronrod_21p[idx] * (std::abs(q_fv1[idx] - reskh) + std::abs(q_fv2[idx] - reskh));
    }
    result = resk * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);

    if (resasc != 0.0 && abserr != 0.0)
        abserr = resasc * std::min(1., std::pow(200.0 * abserr / resasc, 1.5));

    if (resabs > q_uflow / (50. * q_epmach))
        abserr = std::max((q_epmach * 50.0) * resabs, abserr);

    qk21_output.Value     = result;
    qk21_output.Error     = abserr;
    qk21_abs_output.Value = resabs;
    qk21_abs_output.Error = resasc;
    return std::make_pair(qk21_output, qk21_abs_output);
}

int Integral::q_sort(int q_limit, int last, int maxerr, int nrmax, const std::vector<double>& q_elist)
{
    // Check whether the list contains more than two error estimates.
    if (last <= 2)
    {
        q_iord_[0] = 1;
        q_iord_[1] = 2;
        return nrmax;
    }

    double errmax, errmin;
    int i, ibeg, isucc, jbnd, jupbn, k, nrmax_intern;

    nrmax_intern = nrmax;

    // This part of the routine is only executed if, due to a difficult integrand,
    // subdivision increased the error estimate. in the normal case the insert procedure should
    // start after the nrmax-th largest error estimate.

    errmax = q_elist[maxerr - 1];

    if (nrmax_intern != 1)
    {
        int max_iter = nrmax_intern - 1;
        for (int idx = 1; idx <= max_iter; idx++)
        {
            isucc = q_iord_[nrmax_intern - 2];

            if (q_elist[isucc - 1] <= errmax)
                break;

            q_iord_[nrmax_intern - 1] = isucc;
            nrmax_intern              = nrmax_intern - 1;
        }
        // while(nrmax_intern > 1)
        // {
        //     isucc = q_iord_[nrmax_intern-2];

        //     if(q_elist[isucc-1] <= errmax)
        //         break;

        //     q_iord_[nrmax_intern-1] = isucc;
        //     nrmax_intern = nrmax_intern - 1;
        // }
    }

    // Compute the number of elements in the list to be maintained in descending order.
    // This number depends on the number of subdivisions still allowed.

    jupbn = last;
    if ((q_limit / 2 + 2) < last)
        jupbn = q_limit + 3 - last;

    errmin = q_elist[last - 1];

    // Insert errmax by traversing the list top-down, starting
    // comparison from the element elist(iord(nrmax+1)).

    jbnd = jupbn - 1;
    ibeg = nrmax_intern + 1;

    if (ibeg > jbnd)
    {
        q_iord_[jbnd - 1]  = maxerr;
        q_iord_[jupbn - 1] = last;
        return nrmax_intern;
    }

    bool jump_to_60 = false;
    for (i = ibeg; i <= jbnd; i++)
    {
        isucc = q_iord_[i - 1];
        if (q_elist[isucc - 1] <= errmax)
        {
            jump_to_60 = true;
            break;
        }

        q_iord_[i - 2] = isucc;
    }

    if (!jump_to_60)
    {
        q_iord_[jbnd - 1]  = maxerr;
        q_iord_[jupbn - 1] = last;
        return nrmax_intern;
    }

    // Insert errmin by traversing the list bottom-up.
    q_iord_[i - 2] = maxerr;
    k              = jbnd;

    for (int idx = i; idx <= jbnd; idx++)
    {
        isucc = q_iord_[k - 1];
        if (errmin < q_elist[isucc - 1])
        {
            q_iord_[k] = last;
            return nrmax_intern;
        }
        q_iord_[k] = isucc;
        k          = k - 1;
    }
    // while(k > 0):
    // {
    //     isucc = q_iord_[k-1];
    //     if(errmin < q_elist[isucc-1])
    //     {
    //         q_iord_[k] = last;
    //         return nrmax_intern;
    //     }
    //     q_iord_[k] = isucc;
    //     k = k-1;
    // }
    q_iord_[i - 1] = last;

    return nrmax_intern;
}

Integral::InterpolationResults Integral::q_epsilon_extrapolation(int q_limit_epsilon_table, int numrl2, int nres)
{
    const double q_epmach = std::numeric_limits<double>::epsilon(); // machine epsilon
    const double q_oflow  = std::numeric_limits<double>::max();     // largest finite value
    // output
    Integral::InterpolationResults extrapolation_output;
    double result, abserr;
    abserr = q_oflow;
    result = q_rlist2_[numrl2 - 1];
    if (numrl2 < 3)
    {
        abserr = std::max(abserr, 5.0 * q_epmach * std::abs(result));
        // abserr = std::max(abserr, 0.5 * q_epmach * std::abs(result)); diff between f90 & f77
        extrapolation_output.Value = result;
        extrapolation_output.Error = abserr;
        return extrapolation_output;
    }

    double delta1, delta2, delta3, epsinf, error, err1, err2, err3, e0, e1, e1abs, e2, e3, res, ss, tol1, tol2, tol3;
    int ib, k1, k2, k3, newelm, num;

    q_rlist2_[numrl2 + 1] = q_rlist2_[numrl2 - 1];
    newelm                = (numrl2 - 1) / 2;
    q_rlist2_[numrl2 - 1] = q_oflow;
    num                   = numrl2;
    k1                    = numrl2;

    for (int idx = 1; idx <= newelm; idx++)
    {
        k2     = k1 - 1;
        k3     = k1 - 2;
        res    = q_rlist2_[k1 + 1];
        e0     = q_rlist2_[k3 - 1];
        e1     = q_rlist2_[k2 - 1];
        e2     = res;
        e1abs  = std::abs(e1);
        delta2 = e2 - e1;
        err2   = std::abs(delta2);
        tol2   = std::max(std::abs(e2), e1abs) * q_epmach;
        delta3 = e1 - e0;
        err3   = std::abs(delta3);
        tol3   = std::max(e1abs, std::abs(e0)) * q_epmach;

        // If e0, e1 and e2 are equal to within machine accuracy, convergence is assumed.
        if (err2 <= tol2 && err3 <= tol3)
        {
            result = res;
            abserr = err2 + err3;

            abserr = std::max(abserr, 5.0 * q_epmach * std::abs(result));
            // abserr = std::max(abserr, 0.5 * q_epmach * std::abs(result)); diff between f90 & f77
            extrapolation_output.Value = result;
            extrapolation_output.Error = abserr;
            return extrapolation_output;
        }

        e3                = q_rlist2_[k1 - 1];
        q_rlist2_[k1 - 1] = e1;
        delta1            = e1 - e3;
        err1              = std::abs(delta1);
        tol1              = std::max(e1abs, std::abs(e3)) * q_epmach;

        // If two elements are very close to each other, omit a part
        // of the table by adjusting the value of N.
        if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3)
        {
            numrl2 = idx + idx - 1;
            break;
        }

        ss     = 1.0 / delta1 + 1.0 / delta2 - 1.0 / delta3;
        epsinf = std::abs(ss * e1);

        if (epsinf <= 1.0e-4)
        {
            numrl2 = idx + idx - 1;
            break;
        }

        res               = e1 + 1.0 / ss;
        q_rlist2_[k1 - 1] = res;
        k1                = k1 - 2;
        error             = err2 + std::abs(res - e2) + err3;

        if (error <= abserr)
        {
            abserr = error;
            result = res;
        }
    }

    // Shift the table.
    if (numrl2 == q_limit_epsilon_table)
        numrl2 = 2 * (q_limit_epsilon_table / 2) - 1;

    ib = 1;
    if ((num / 2) * 2 == num)
        ib = 2;

    for (int idx = 1; idx <= newelm + 1; idx++)
    {
        q_rlist2_[ib - 1] = q_rlist2_[ib + 1];
        ib                = ib + 2;
    }

    if (num != numrl2)
    {
        for (int idx = 0; idx < numrl2; idx++)
        {
            q_rlist2_[idx] = q_rlist2_[num - numrl2 + idx];
        }
    }

    if (nres < 4)
    {
        q_last_3_results_[nres - 1] = result;
        abserr                      = q_oflow;
    } else
    {
        abserr = std::abs(result - q_last_3_results_[2]) + std::abs(result - q_last_3_results_[1]) +
                 std::abs(result - q_last_3_results_[0]);
        q_last_3_results_[0] = q_last_3_results_[1];
        q_last_3_results_[1] = q_last_3_results_[2];
        q_last_3_results_[2] = result;
    }

    abserr = std::max(abserr, 5.0 * q_epmach * std::abs(result));
    // abserr = std::max(abserr, 0.5 * q_epmach * std::abs(result)); diff between f90 & f77
    extrapolation_output.Value = result;
    extrapolation_output.Error = abserr;
    return extrapolation_output;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Integral::SetIntegrand(std::function<double(double)> integrand)
{
    integrand_ = integrand;
}

void Integral::SetMax(double max)
{
    max_ = max;
}

void Integral::SetMaxStepsUpperLimit(int maxSteps)
{
    maxSteps_upper_limit_ = maxSteps;
}

void Integral::SetMaxStepsRomberg(int maxSteps)
{
    maxSteps_romberg_ = maxSteps;
}

void Integral::SetMin(double min)
{
    min_ = min;
}

void Integral::SetPowerOfSubstitution(double powerOfSubstitution)
{
    powerOfSubstitution_ = powerOfSubstitution;
}

void Integral::SetPrecision(double precision)
{
    precision_ = precision;
}

void Integral::SetRandomDo(bool randomDo)
{
    randomDo_ = randomDo;
}

void Integral::SetRandomNumber(double randomNumber)
{
    randomNumber_ = randomNumber;
}

void Integral::SetRandomX(double randomX)
{
    randomX_ = randomX;
}

void Integral::SetReverse(bool reverse)
{
    reverse_ = reverse;
}

void Integral::SetReverseX(double reverseX)
{
    reverseX_ = reverseX;
}

void Integral::SetRomberg(int romberg)
{
    romberg_ = romberg;
}

void Integral::SetRomberg4refine(int romberg4refine)
{
    romberg4refine_ = romberg4refine;
}

void Integral::SetSavedResult(double savedResult)
{
    savedResult_ = savedResult;
}

void Integral::SetUseLog(bool useLog)
{
    useLog_ = useLog;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Integral::~Integral() {}
