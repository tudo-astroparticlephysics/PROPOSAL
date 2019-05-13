/*! \file   Interpolant.cxx
 *   \brief  Source file for the interpolation routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   02.08.2010
 *   \author Jan-Hendrik Koehne
 */

/**
 * This class provides routines for function interpolation. Include the function to be interpolated
 * in a class that implements the interface FunctionInt or FunctionInt2 (defined below).
 * Methods contained here are based on the Numerical Recipes (W. H. Press et al.).
 * <pre>
 * interface FunctionInt{
 *     double functionInt(double x);
 * }
 *
 * interface FunctionInt2{
 *     double functionInt(double x1, double x2);
 * }
 * </pre>
 * @author Dmitry Chirkin
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

const double Interpolant::bigNumber_  = -300;
const double Interpolant::aBigNumber_ = -299;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::Interpolate(double x)
{
    int start;
    double result, aux;

    if (isLog_)
    {
        x = Slog(x);
    }

    reverse_ = true;
    aux      = (x - xmin_) / step_;
    starti_  = (int)aux;

    if (starti_ < 0)
    {
        starti_ = 0;
    } else if (starti_ >= max_)
    {
        starti_ = max_ - 1;
    }

    start = (int)(aux - 0.5 * (romberg_ - 1));

    if (start < 0)
    {
        start = 0;
    } else if (start + romberg_ > max_ || start > max_)
    {
        start = max_ - romberg_;
    }
    result = Interpolate(x, start);

    if (logSubst_)
    {
        if (self_)
        {
            result = Exp(result);
        }
    }

    return result;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::Interpolate(double x1, double x2)
{
    int i, start;
    double aux, aux2 = 0, result;

    if (isLog_)
    {
        x2 = std::log(x2);
    }

    reverse_ = true;

    aux     = (x2 - xmin_) / step_;
    starti_ = (int)aux;

    if (starti_ < 0)
    {
        starti_ = 0;
    } else if (starti_ >= max_)
    {
        starti_ = max_ - 1;
    }

    start = (int)(aux - 0.5 * (romberg_ - 1));

    if (start < 0)
    {
        start = 0;
    } else if (start + romberg_ > max_ || start > max_)
    {
        start = max_ - romberg_;
    }

    for (i = start; i < start + romberg_; i++)
    {
        iY_.at(i) = Interpolant_.at(i)->Interpolate(x1);
    }

    if (!fast_)
    {
        aux = 0;

        for (i = start; i < start + romberg_; i++)
        {
            if (Interpolant_.at(i)->precision_ > aux)
            {
                aux  = Interpolant_.at(i)->precision_;
                aux2 = Interpolant_.at(i)->worstX_;
            }
        }

        if (aux > precision2_)
        {
            precision2_ = aux;
            worstX2_    = aux2;
        }
    }

    result = Interpolate(x2, start);

    if (logSubst_)
    {
        result = Exp(result);
    }

    return result;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::InterpolateArray(double x)
{
    int i, j, m, start, auxdir;
    bool dir;

    reverse_ = false;
    i        = 0;
    j        = max_ - 1;
    dir      = iX_.at(max_ - 1) > iX_.at(0);

    while (j - i > 1)
    {
        m = (i + j) / 2;

        if ((x > iX_.at(m)) == dir)
        {
            i = m;
        } else
        {
            j = m;
        }
    }

    if (i + 1 < max_)
    {
        if (((x - iX_.at(i)) < (iX_.at(i + 1) - x)) == dir)
        {
            auxdir = 0;
        } else
        {
            auxdir = 1;
        }
    } else
    {
        auxdir = 0;
    }

    starti_ = i + auxdir;
    start   = i - (int)(0.5 * (romberg_ - 1 - auxdir));

    if (start < 0)
    {
        start = 0;
    }

    if (start + romberg_ > max_ || start > max_)
    {
        start = max_ - romberg_;
    }

    return Interpolate(x, start);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::InterpolateArray(double x1, double x2)
{
    int i, j, m, start, auxdir, aux, aux2;
    bool dir;

    reverse_ = false;
    i        = 0;
    j        = max_ - 1;
    dir      = iX_.at(max_ - 1) > iX_.at(0);

    while (j - i > 1)
    {
        m = (i + j) / 2;

        if ((x1 > iX_.at(m)) == dir)
        {
            i = m;
        } else
        {
            j = m;
        }
    }

    if (i + 1 < max_)
    {
        if (((x1 - iX_.at(i)) < (iX_.at(i + 1) - x1)) == dir)
        {
            auxdir = 0;
        } else
        {
            auxdir = 1;
        }
    } else
    {
        auxdir = 0;
    }

    starti_ = i + auxdir;
    start   = i - (int)(0.5 * (romberg_ - 1 - auxdir));

    if (start < 0)
    {
        start = 0;
    }

    if (start + romberg_ > max_ || start > max_)
    {
        start = max_ - romberg_;
    }

    for (i = start; i < start + romberg_; i++)
    {
        iY_.at(i) = Interpolant_.at(i)->InterpolateArray(x2);
    }

    if (!fast_)
    {
        aux = 0;
        aux2 = 0;

        for (i = start; i < start + romberg_; i++)
        {
            if (Interpolant_.at(i)->precision_ > aux)
            {
                aux  = Interpolant_.at(i)->precision_;
                aux2 = Interpolant_.at(i)->worstX_;
            }
        }

        if (aux > precision2_)
        {
            precision2_ = aux;
            worstX2_    = aux2;
        }
    }

    double result = Interpolate(x1, start);

    return result;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::FindLimit(double y)
{
    int i, j, m, start, auxdir;
    bool dir, rat, rel;
    double result;

    rat      = false;
    reverse_ = false;

    if (logSubst_)
    {
        y = Log(y);
    }

    i   = 0;
    j   = max_ - 1;
    dir = iY_.at(max_ - 1) > iY_.at(0);

    while (j - i > 1)
    {
        m = (i + j) / 2;

        if ((y > iY_.at(m)) == dir)
        {
            i = m;
        } else
        {
            j = m;
        }
    }

    iX_.swap(iY_);
    if (!fast_)
    {
        std::swap(precision_, precisionY_);
        std::swap(worstX_, worstY_);

        rat       = rational_;
        rational_ = rationalY_;
    }

    std::swap(romberg_, rombergY_);

    rel       = relative_;
    relative_ = relativeY_;

    if (i + 1 < max_)
    {
        if (((y - iX_.at(i)) < (iX_.at(i + 1) - y)) == dir)
        {
            auxdir = 0;
        } else
        {
            auxdir = 1;
        }
    } else
    {
        auxdir = 0;
    }

    starti_ = i + auxdir;
    start   = i - (int)(0.5 * (romberg_ - 1 - auxdir));

    if (start < 0)
    {
        start = 0;
    }

    if (start + romberg_ > max_ || start > max_)
    {
        start = max_ - romberg_;
    }

    result = Interpolate(y, start);

    iX_.swap(iY_);

    if (!fast_)
    {
        std::swap(precision_, precisionY_);
        std::swap(worstX_, worstY_);

        rational_ = rat;
    }

    std::swap(romberg_, rombergY_);

    relative_ = rel;

    if (result < xmin_)
    {
        result = xmin_;
    } else if (result > xmax_)
    {
        result = xmax_;
    }

    if (isLog_)
    {
        result = Exp(result);
    }

    return result;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::FindLimit(double x1, double y)
{
    int i, j, m, start, auxdir;
    bool dir;
    double result, aux, aux2 = 0;

    reverse_ = false;

    if (logSubst_)
    {
        y = Log(y);
    }

    if (!flag_)
    {
        for (i = 0; i < max_; i++)
        {
            iY_.at(i) = Interpolant_.at(i)->Interpolate(x1);
        }
    }

    i = 0;
    j = max_ - 1;

    if (flag_)
    {
        dir = Interpolant_.at(max_ - 1)->Interpolate(x1) > Interpolant_.at(0)->Interpolate(x1);
    } else
    {
        dir = iY_.at(max_ - 1) > iY_.at(0);
    }

    while (j - i > 1)
    {
        m = (i + j) / 2;

        if (flag_)
        {
            aux = Interpolant_.at(m)->Interpolate(x1);
        } else
        {
            aux = iY_.at(m);
        }

        if ((y > aux) == dir)
        {
            i = m;
        } else
        {
            j = m;
        }
    }

    iX_.swap(iY_);

    if (!fast_)
    {
        std::swap(precision_, precisionY_);
        std::swap(worstX_, worstY_);
        std::swap(rational_, rationalY_);
    }

    std::swap(romberg_, rombergY_);
    std::swap(relative_, relativeY_);

    if (i + 1 < max_)
    {
        // iX_ is not pre-filled if flag is true. Fill the bits we're about to use.
        if (flag_)
        {
            iX_.at(i) = Interpolant_.at(i)->Interpolate(x1);
            iX_.at(i+1) = Interpolant_.at(i+1)->Interpolate(x1);
        }
        if (((y - iX_.at(i)) < (iX_.at(i + 1) - y)) == dir)
        {
            auxdir = 0;
        } else
        {
            auxdir = 1;
        }
    } else
    {
        auxdir = 0;
    }

    starti_ = i + auxdir;
    start   = i - (int)(0.5 * (romberg_ - 1 - auxdir));

    if (start < 0)
    {
        start = 0;
    }

    if (start + romberg_ > max_ || start > max_)
    {
        start = max_ - romberg_;
    }

    if (flag_)
    {
        for (i = start; i < start + romberg_; i++)
        {
            iX_.at(i) = Interpolant_.at(i)->Interpolate(x1);
        }
    }

    result = Interpolate(y, start);

    iX_.swap(iY_);

    if (!fast_)
    {
        std::swap(precision_, precisionY_);
        std::swap(worstX_, worstY_);
        std::swap(rational_, rationalY_);
    }

    std::swap(romberg_, rombergY_);
    std::swap(relative_, relativeY_);

    if (result < xmin_)
    {
        result = xmin_;
    } else if (result > xmax_)
    {
        result = xmax_;
    }

    if (!fast_)
    {
        aux = 0;

        for (i = start; i < start + romberg_; i++)
        {
            if (Interpolant_.at(i)->precision_ > aux)
            {
                aux  = Interpolant_.at(i)->precision_;
                aux2 = Interpolant_.at(i)->worstX_;
            }
        }

        if (aux > precision2_)
        {
            precision2_ = aux;
            worstX2_    = aux2;
        }
    }

    if (isLog_)
    {
        result = std::exp(result);
    }

    return result;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------Save and Load-------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool Interpolant::Save(std::string Path, bool binary_tables)
{
    std::ofstream out;

    if (!binary_tables)
    {
        std::stringstream ss;
        ss << Path << ".txt";

        out.open(ss.str().c_str());

        out.precision(16);

        if (!out.good())
        {
            log_error("Can not open file %s for writing", Path.c_str());
            return 0;
        }
    } else
    {

        out.open(Path.c_str(), std::ios::binary);

        if (!out.good())
        {
            log_error("Can not open file %s for writing", Path.c_str());
            return 0;
        }
    }

    Save(out, binary_tables);

    out.close();
    return 1;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool Interpolant::Save(std::ofstream& out, bool binary_tables)
{
    if (!out.good())
    {
        log_error("Can not open file for writing");
        return 0;
    }
    bool D2 = false;
    if (function2d_ != NULL)
        D2 = true;

    if (binary_tables)
    {
        out.write(reinterpret_cast<char*>(&D2), sizeof D2);

        if (D2)
        {
            if (isLog_)
            {
                double xmax = std::exp(xmax_);
                double xmin = std::exp(xmin_);

                out.write(reinterpret_cast<char*>(&max_), sizeof max_);
                out.write(reinterpret_cast<char*>(&xmin), sizeof xmin);
                out.write(reinterpret_cast<char*>(&xmax), sizeof xmax);
            } else
            {
                out.write(reinterpret_cast<char*>(&max_), sizeof max_);
                out.write(reinterpret_cast<char*>(&xmin_), sizeof xmin_);
                out.write(reinterpret_cast<char*>(&xmax_), sizeof xmax_);
            }
            out.write(reinterpret_cast<char*>(&romberg_), sizeof romberg_);
            out.write(reinterpret_cast<char*>(&rational_), sizeof rational_);
            out.write(reinterpret_cast<char*>(&relative_), sizeof relative_);
            out.write(reinterpret_cast<char*>(&isLog_), sizeof isLog_);
            out.write(reinterpret_cast<char*>(&rombergY_), sizeof rombergY_);
            out.write(reinterpret_cast<char*>(&rationalY_), sizeof rationalY_);
            out.write(reinterpret_cast<char*>(&relativeY_), sizeof relativeY_);
            out.write(reinterpret_cast<char*>(&logSubst_), sizeof logSubst_);

            for (int i = 0; i < max_; i++)
            {
                out.write(reinterpret_cast<char*>(&iX_.at(i)), sizeof iX_.at(i));
                Interpolant_.at(i)->Save(out, binary_tables);
            }
        } else
        {
            if (isLog_)
            {
                double xmax = std::exp(xmax_);
                double xmin = std::exp(xmin_);

                out.write(reinterpret_cast<char*>(&max_), sizeof max_);
                out.write(reinterpret_cast<char*>(&xmin), sizeof xmin);
                out.write(reinterpret_cast<char*>(&xmax), sizeof xmax);
            } else
            {
                out.write(reinterpret_cast<char*>(&max_), sizeof max_);
                out.write(reinterpret_cast<char*>(&xmin_), sizeof xmin_);
                out.write(reinterpret_cast<char*>(&xmax_), sizeof xmax_);
            }

            out.write(reinterpret_cast<char*>(&romberg_), sizeof romberg_);
            out.write(reinterpret_cast<char*>(&rational_), sizeof rational_);
            out.write(reinterpret_cast<char*>(&relative_), sizeof relative_);
            out.write(reinterpret_cast<char*>(&isLog_), sizeof isLog_);
            out.write(reinterpret_cast<char*>(&rombergY_), sizeof rombergY_);
            out.write(reinterpret_cast<char*>(&rationalY_), sizeof rationalY_);
            out.write(reinterpret_cast<char*>(&relativeY_), sizeof relativeY_);
            out.write(reinterpret_cast<char*>(&logSubst_), sizeof logSubst_);

            for (int i = 0; i < max_; i++)
            {
                out.write(reinterpret_cast<char*>(&iX_.at(i)), sizeof iX_.at(i));
                out.write(reinterpret_cast<char*>(&iY_.at(i)), sizeof iY_.at(i));
            }
        }
    } else
    {
        out << D2 << std::endl;

        if (D2)
        {
            if (isLog_)
            {
                out << max_ << "\t" << std::exp(xmin_) << "\t" << std::exp(xmax_) << std::endl;
            } else
            {
                out << max_ << "\t" << xmin_ << "\t" << xmax_ << std::endl;
            }
            out << romberg_ << "\t" << rational_ << "\t" << relative_ << "\t" << isLog_ << std::endl;
            out << rombergY_ << "\t" << rationalY_ << "\t" << relativeY_ << "\t" << logSubst_ << std::endl;

            for (int i = 0; i < max_; i++)
            {
                out << iX_.at(i) << std::endl;
                Interpolant_.at(i)->Save(out, binary_tables);
            }
        } else
        {
            if (isLog_)
            {
                out << max_ << "\t" << std::exp(xmin_) << "\t" << std::exp(xmax_) << std::endl;
            } else
            {
                out << max_ << "\t" << xmin_ << "\t" << xmax_ << std::endl;
            }
            out << romberg_ << "\t" << rational_ << "\t" << relative_ << "\t" << isLog_ << std::endl;
            out << rombergY_ << "\t" << rationalY_ << "\t" << relativeY_ << "\t" << logSubst_ << std::endl;

            for (int i = 0; i < max_; i++)
            {
                out << iX_.at(i) << "\t" << iY_.at(i) << std::endl;
            }
        }
    }
    out.flush();
    return 1;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool Interpolant::Load(std::string Path, bool binary_tables)
{
    bool success;
    std::ifstream in;

    if (!binary_tables)
    {
        std::stringstream ss;
        ss << Path << ".txt";

        in.open(ss.str().c_str());
    } else
    {
        in.open(Path.c_str(), std::ios::binary);
    }

    success = Load(in, binary_tables);

    in.close();
    return success;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool Interpolant::Load(std::ifstream& in, bool binary_tables)
{
    bool D2;

    int max;
    double xmin, xmax;
    int romberg, rombergY;
    bool rational, rationalY, relative, relativeY, isLog, logSubst;

    if (binary_tables)
    {
        in.read(reinterpret_cast<char*>(&D2), sizeof D2);

        if (D2)
        {
            if (!in.good())
                return 0;

            in.read(reinterpret_cast<char*>(&max), sizeof max);
            in.read(reinterpret_cast<char*>(&xmin), sizeof xmin);
            in.read(reinterpret_cast<char*>(&xmax), sizeof xmax);
            in.read(reinterpret_cast<char*>(&romberg), sizeof romberg);
            in.read(reinterpret_cast<char*>(&rational), sizeof rational);
            in.read(reinterpret_cast<char*>(&relative), sizeof relative);
            in.read(reinterpret_cast<char*>(&isLog), sizeof isLog);
            in.read(reinterpret_cast<char*>(&rombergY), sizeof rombergY);
            in.read(reinterpret_cast<char*>(&rationalY), sizeof rationalY);
            in.read(reinterpret_cast<char*>(&relativeY), sizeof relativeY);
            in.read(reinterpret_cast<char*>(&logSubst), sizeof logSubst);

            InitInterpolant(
                max, xmin, xmax, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);

            Interpolant_.resize(max_);

            for (int i = 0; i < max_; i++)
            {
                in.read(reinterpret_cast<char*>(&iX_.at(i)), sizeof iX_.at(i));

                if (!in.good())
                    return 0;
                Interpolant_.at(i) = new Interpolant();
                Interpolant_.at(i)->Load(in, binary_tables);
                Interpolant_.at(i)->self_ = false;
            }

        } else
        {
            if (!in.good())
                return 0;
            in.read(reinterpret_cast<char*>(&max), sizeof max);
            in.read(reinterpret_cast<char*>(&xmin), sizeof xmin);
            in.read(reinterpret_cast<char*>(&xmax), sizeof xmax);
            in.read(reinterpret_cast<char*>(&romberg), sizeof romberg);
            in.read(reinterpret_cast<char*>(&rational), sizeof rational);
            in.read(reinterpret_cast<char*>(&relative), sizeof relative);
            in.read(reinterpret_cast<char*>(&isLog), sizeof isLog);
            in.read(reinterpret_cast<char*>(&rombergY), sizeof rombergY);
            in.read(reinterpret_cast<char*>(&rationalY), sizeof rationalY);
            in.read(reinterpret_cast<char*>(&relativeY), sizeof relativeY);
            in.read(reinterpret_cast<char*>(&logSubst), sizeof logSubst);

            InitInterpolant(
                max, xmin, xmax, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);

            for (int i = 0; i < max_; i++)
            {
                in.read(reinterpret_cast<char*>(&iX_.at(i)), sizeof iX_.at(i));
                in.read(reinterpret_cast<char*>(&iY_.at(i)), sizeof iY_.at(i));
                if (!in.good())
                    return 0;
            }
        }
    } else
    {
        in >> D2;

        if (D2)
        {
            if (!in.good())
                return 0;
            in >> max >> xmin >> xmax;
            in >> romberg >> rational >> relative >> isLog;
            in >> rombergY >> rationalY >> relativeY >> logSubst;

            InitInterpolant(
                max, xmin, xmax, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);

            Interpolant_.resize(max_);

            for (int i = 0; i < max_; i++)
            {
                in >> iX_.at(i);
                if (!in.good())
                    return 0;
                Interpolant_.at(i) = new Interpolant();
                Interpolant_.at(i)->Load(in, binary_tables);
                Interpolant_.at(i)->self_ = false;
            }

        } else
        {
            if (!in.good())
                return 0;
            in >> max >> xmin >> xmax;
            in >> romberg >> rational >> relative >> isLog;
            in >> rombergY >> rationalY >> relativeY >> logSubst;

            InitInterpolant(
                max, xmin, xmax, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);

            for (int i = 0; i < max_; i++)
            {
                in >> iX_.at(i) >> iY_.at(i);
                if (!in.good())
                    return 0;
            }
        }
    }
    return 1;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Interpolant::Interpolant()
    : romberg_(1.)
    , rombergY_(1.)
    , iX_()
    , iY_()
    , c_()
    , d_()
    , max_(1.)
    , xmin_(1.)
    , xmax_(1.)
    , step_(0)
    , rational_(false)
    , relative_(false)
    , function1d_(NULL)
    , function2d_(NULL)
    , Interpolant_()
    , row_(0)
    , starti_(0)
    , rationalY_(false)
    , relativeY_(false)
    , reverse_(false)
    , self_(true)
    , flag_(false)
    , isLog_(false)
    , logSubst_(false)
    , precision_(0)
    , worstX_(0)
    , precision2_(0)
    , worstX2_(0)
    , precisionY_(0)
    , worstY_(0)
    , fast_(true)
    , x_save_(0)
    , y_save_(0)
{
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Interpolant::Interpolant(const Interpolant& interpolant)
    : romberg_(interpolant.romberg_)
    , rombergY_(interpolant.rombergY_)
    , iX_(interpolant.iX_)
    , iY_(interpolant.iY_)
    , c_(interpolant.c_)
    , d_(interpolant.d_)
    , max_(interpolant.max_)
    , xmin_(interpolant.xmin_)
    , xmax_(interpolant.xmax_)
    , step_(interpolant.step_)
    , rational_(interpolant.rational_)
    , relative_(interpolant.relative_)
    , row_(interpolant.row_)
    , starti_(interpolant.starti_)
    , rationalY_(interpolant.rationalY_)
    , relativeY_(interpolant.relativeY_)
    , reverse_(interpolant.reverse_)
    , self_(interpolant.self_)
    , flag_(interpolant.flag_)
    , isLog_(interpolant.isLog_)
    , logSubst_(interpolant.logSubst_)
    , precision_(interpolant.precision_)
    , worstX_(interpolant.worstX_)
    , precision2_(interpolant.precision2_)
    , worstX2_(interpolant.worstX2_)
    , precisionY_(interpolant.precisionY_)
    , worstY_(interpolant.worstY_)
    , fast_(interpolant.fast_)
    , x_save_(interpolant.x_save_)
    , y_save_(interpolant.y_save_)

{
    Interpolant_.resize(interpolant.Interpolant_.size());

    for (unsigned int i = 0; i < interpolant.Interpolant_.size(); i++)
    {
        Interpolant_.at(i) = new Interpolant(*interpolant.Interpolant_.at(i));
    }

    function1d_ = std::ref(interpolant.function1d_);
    function2d_ = std::ref(interpolant.function2d_);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Interpolant::Interpolant(int max,
                         double xmin,
                         double xmax,
                         std::function<double(double)> function1d,
                         int romberg,
                         bool rational,
                         bool relative,
                         bool isLog,
                         int rombergY,
                         bool rationalY,
                         bool relativeY,
                         bool logSubst)
    : romberg_(1.)
    , rombergY_(1.)
    , iX_()
    , iY_()
    , c_()
    , d_()
    , max_(1.)
    , xmin_(1.)
    , xmax_(1.)
    , step_(0)
    , rational_(false)
    , relative_(false)
    , function1d_(NULL)
    , function2d_(NULL)
    , Interpolant_()
    , row_(0)
    , starti_(0)
    , rationalY_(false)
    , relativeY_(false)
    , reverse_(false)
    , self_(true)
    , flag_(false)
    , isLog_(false)
    , logSubst_(false)
    , precision_(0)
    , worstX_(0)
    , precision2_(0)
    , worstX2_(0)
    , precisionY_(0)
    , worstY_(0)
    , fast_(true)
    , x_save_(1)
    , y_save_(0)
{
    InitInterpolant(max, xmin, xmax, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);

    int i;
    double aux, xaux;

    function1d_ = function1d;

    aux = this->xmin_ + step_ / 2;

    for (i = 0; i < max_; i++)
    {
        iX_.at(i) = aux;

        if (isLog_)
        {
            xaux = std::exp(aux);
        } else
        {
            xaux = aux;
        }

        iY_.at(i) = function1d_(xaux);

        if (logSubst_)
        {
            iY_.at(i) = Log(iY_.at(i));
        }

        aux += step_;
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Interpolant::Interpolant(int max1,
                         double x1min,
                         double x1max,
                         int max2,
                         double x2min,
                         double x2max,
                         std::function<double(double, double)> function2d,
                         int romberg1,
                         bool rational1,
                         bool relative1,
                         bool isLog1,
                         int romberg2,
                         bool rational2,
                         bool relative2,
                         bool isLog2,
                         int rombergY,
                         bool rationalY,
                         bool relativeY,
                         bool logSubst)
    : romberg_(1.)
    , rombergY_(1.)
    , iX_()
    , iY_()
    , c_()
    , d_()
    , max_(1.)
    , xmin_(1.)
    , xmax_(1.)
    , step_(0)
    , rational_(false)
    , relative_(false)
    , function1d_(NULL)
    , function2d_(NULL)
    , Interpolant_()
    , row_(0)
    , starti_(0)
    , rationalY_(false)
    , relativeY_(false)
    , reverse_(false)
    , self_(true)
    , flag_(false)
    , isLog_(false)
    , logSubst_(false)
    , precision_(0)
    , worstX_(0)
    , precision2_(0)
    , worstX2_(0)
    , precisionY_(0)
    , worstY_(0)
    , fast_(true)
    , x_save_(1)
    , y_save_(0)
{
    InitInterpolant(
        max2, x2min, x2max, romberg2, rational2, relative2, isLog2, rombergY, rationalY, relativeY, logSubst);

    int i;
    double aux;

    function2d_ = function2d;
    function1d_ = std::bind(&Interpolant::Get2dFunctionFixedY, this, std::placeholders::_1);

    Interpolant_.resize(max_);

    for (i = 0, aux = xmin_ + step_ / 2; i < max_; i++, aux += step_)
    {
        iX_.at(i) = aux;
        row_      = i;

        Interpolant_.at(i) = new Interpolant(max1,
                                             x1min,
                                             x1max,
                                             function1d_,
                                             romberg1,
                                             rational1,
                                             relative1,
                                             isLog1,
                                             rombergY,
                                             rationalY,
                                             relativeY,
                                             logSubst_);

        Interpolant_.at(i)->self_ = false;
    }

    precision2_ = 0;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Interpolant::Interpolant(std::vector<double> x, std::vector<double> y, int romberg, bool rational, bool relative)
    : romberg_(1.)
    , rombergY_(1.)
    , iX_()
    , iY_()
    , c_()
    , d_()
    , max_(1.)
    , xmin_(1.)
    , xmax_(1.)
    , step_(0)
    , rational_(false)
    , relative_(false)
    , function1d_(NULL)
    , function2d_(NULL)
    , Interpolant_()
    , row_(0)
    , starti_(0)
    , rationalY_(false)
    , relativeY_(false)
    , reverse_(false)
    , self_(true)
    , flag_(false)
    , isLog_(false)
    , logSubst_(false)
    , precision_(0)
    , worstX_(0)
    , precision2_(0)
    , worstX2_(0)
    , precisionY_(0)
    , worstY_(0)
    , fast_(true)
    , x_save_(1)
    , y_save_(0)
{
    InitInterpolant(std::min(x.size(), y.size()),
                    x.at(0),
                    x.at(x.size() - 1),
                    romberg,
                    rational,
                    relative,
                    false,
                    romberg,
                    rational,
                    relative,
                    false);

    if (x.size() != y.size())
    {
        log_fatal("size of x(%i) and y(%i) do not match!", (int)(x.size()), (int)(y.size()));
    }

    for (int i = 0; i < (int)x.size(); i++)
    {
        iX_.at(i) = x.at(i);
    }

    for (int i = 0; i < (int)y.size(); i++)
    {
        iY_.at(i) = y.at(i);
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Interpolant::Interpolant(std::vector<double> x1, std::vector<double> x2, std::vector< std::vector<double> > y, int romberg1, bool rational1, bool relative1 , int romberg2, bool rational2, bool relative2)
        : romberg_(1.)
        , rombergY_(1.)
        , iX_()
        , iY_()
        , iY2_()
        , c_()
        , d_()
        , max_(1.)
        , xmin_(1.)
        , xmax_(1.)
        , step_(0)
        , rational_(false)
        , relative_(false)
        , function1d_(NULL)
        , function2d_(NULL)
        , Interpolant_()
        , row_(0)
        , starti_(0)
        , rationalY_(false)
        , relativeY_(false)
        , reverse_(false)
        , self_(true)
        , flag_(false)
        , isLog_(false)
        , logSubst_(false)
        , precision_(0)
        , worstX_(0)
        , precision2_(0)
        , worstX2_(0)
        , precisionY_(0)
        , worstY_(0)
        , fast_(true)
        , x_save_(1)
        , y_save_(0)
{

    //TODO: Not sure what is happening in the romberg=0 case
    if (romberg1 <= 0)
    {
        log_warn("romberg1 = %i must be > 0! setting to 1!", romberg1);
        romberg1 = 1;
    }

    if (romberg2 <= 0)
    {
        log_warn("romberg2 = %i must be > 0! setting to 1!", romberg2);
        romberg2 = 1;
    }

    InitInterpolant(std::min(x1.size(), y.size()),
                    x1.at(0),
                    x1.at(x1.size() - 1),
                    romberg1,
                    rational1,
                    relative1,
                    false,
                    romberg1,
                    rational1,
                    relative1,
                    false);

    if (x2.size() != y[0].size())
    {
        log_fatal("size of x2(%i) and y(%i) do not match!", (int)(x2.size()), (int)(y[0].size()));
    }

    if (x1.size() != y.size())
    {
        log_fatal("size of x1(%i) and y(%i) do not match!", (int)(x1.size()), (int)(y.size()));
    }


    Interpolant_.resize((int)x1.size());

    iY2_.resize(y.size());
    for (int i = 0; i < (int)iY2_.size(); i++)
    {
        iY2_[i].resize(y[i].size());
        for(int j=0; j < (int)iY2_[i].size(); j++){
            iY2_.at(i).at(j) = y.at(i).at(j);
        }
    }

    for (int i = 0; i < (int)x1.size(); i++)
    {
        iX_.at(i) = x1.at(i);
        row_      = i;

        Interpolant_.at(i) = new Interpolant(x2,
                                             iY2_[i],
                                             romberg2,
                                             rational2,
                                             relative2);

        Interpolant_.at(i)->self_ = false;
    }

    precision2_ = 0;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Interpolant::Interpolant(std::vector<double> x1, std::vector< std::vector<double> > x2, std::vector< std::vector<double> > y, int romberg1, bool rational1, bool relative1 , int romberg2, bool rational2, bool relative2)
        : romberg_(1.)
        , rombergY_(1.)
        , iX_()
        , iY_()
        , c_()
        , d_()
        , max_(1.)
        , xmin_(1.)
        , xmax_(1.)
        , step_(0)
        , rational_(false)
        , relative_(false)
        , function1d_(NULL)
        , function2d_(NULL)
        , Interpolant_()
        , row_(0)
        , starti_(0)
        , rationalY_(false)
        , relativeY_(false)
        , reverse_(false)
        , self_(true)
        , flag_(false)
        , isLog_(false)
        , logSubst_(false)
        , precision_(0)
        , worstX_(0)
        , precision2_(0)
        , worstX2_(0)
        , precisionY_(0)
        , worstY_(0)
        , fast_(true)
        , x_save_(1)
        , y_save_(0)
{

    //TODO: Not sure what is happening in the romberg=0 case
    if (romberg1 <= 0)
    {
        log_warn("romberg1 = %i must be > 0! setting to 1!", romberg1);
        romberg1 = 1;
    }

    if (romberg2 <= 0)
    {
        log_warn("romberg2 = %i must be > 0! setting to 1!", romberg2);
        romberg2 = 1;
    }

    InitInterpolant(std::min(x1.size(), y.size()),
                    x1.at(0),
                    x1.at(x1.size() - 1),
                    romberg1,
                    rational1,
                    relative1,
                    false,
                    romberg1,
                    rational1,
                    relative1,
                    false);

    for(int i = 0; i < (int)x1.size(); i++){
        if (x2[i].size() != y[i].size())
        {
            log_fatal("size of x2(%i) and y(%i) do not match!", (int)(x2[i].size()), (int)(y[i].size()));
        }
    }

    if (x1.size() != y.size())
    {
        log_fatal("size of x1(%i) and y(%i) do not match!", (int)(x1.size()), (int)(y.size()));
    }

    if (x2.size() != x1.size())
    {
        log_fatal("size of x2(%i) and x1(%i) do not match!", (int)(x2.size()), (int)(x1.size()));
    }


    Interpolant_.resize((int)x1.size());

    iY2_.resize(y.size());
    for (int i = 0; i < (int)iY2_.size(); i++)
    {
        iY2_[i].resize(y[i].size());
        for(int j=0; j < (int)iY2_[i].size(); j++){
            iY2_.at(i).at(j) = y.at(i).at(j);
        }
    }

    for (int i = 0; i < (int)x1.size(); i++)
    {
        iX_.at(i) = x1.at(i);
        row_      = i;

        Interpolant_.at(i) = new Interpolant(x2[i],
                                             iY2_[i],
                                             romberg2,
                                             rational2,
                                             relative2);

        Interpolant_.at(i)->self_ = false;
    }

    precision2_ = 0;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Interpolant& Interpolant::operator=(const Interpolant& interpolant)
{
    if (this != &interpolant)
    {
        Interpolant tmp(interpolant);
        swap(tmp);
    }
    return *this;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool Interpolant::operator==(const Interpolant& interpolant) const
{
    if (romberg_ != interpolant.romberg_)
        return false;
    if (rombergY_ != interpolant.rombergY_)
        return false;
    if (max_ != interpolant.max_)
        return false;
    if (xmin_ != interpolant.xmin_)
        return false;
    if (xmax_ != interpolant.xmax_)
        return false;
    if (step_ != interpolant.step_)
        return false;
    if (rational_ != interpolant.rational_)
        return false;
    if (relative_ != interpolant.relative_)
        return false;
    if (row_ != interpolant.row_)
        return false;
    if (starti_ != interpolant.starti_)
        return false;
    if (rationalY_ != interpolant.rationalY_)
        return false;
    if (relativeY_ != interpolant.relativeY_)
        return false;
    if (reverse_ != interpolant.reverse_)
        return false;
    if (self_ != interpolant.self_)
        return false;
    if (flag_ != interpolant.flag_)
        return false;
    if (isLog_ != interpolant.isLog_)
        return false;
    if (logSubst_ != interpolant.logSubst_)
        return false;
    if (precision_ != interpolant.precision_)
        return false;
    if (worstX_ != interpolant.worstX_)
        return false;
    if (precision2_ != interpolant.precision2_)
        return false;
    if (worstX2_ != interpolant.worstX2_)
        return false;
    if (precisionY_ != interpolant.precisionY_)
        return false;
    if (worstY_ != interpolant.worstY_)
        return false;
    if (fast_ != interpolant.fast_)
        return false;
    if (x_save_ != interpolant.x_save_)
        return false;
    if (y_save_ != interpolant.y_save_)
        return false;

    if (iX_.size() != interpolant.iX_.size())
        return false;
    if (iY_.size() != interpolant.iY_.size())
        return false;
    if (c_.size() != interpolant.c_.size())
        return false;
    if (d_.size() != interpolant.d_.size())
        return false;

    if (Interpolant_.size() != interpolant.Interpolant_.size())
        return false;

    for (unsigned int i = 0; i < iX_.size(); i++)
    {
        if (iX_.at(i) != interpolant.iX_.at(i))
            return false;
    }
    for (unsigned int i = 0; i < iY_.size(); i++)
    {
        if (iY_.at(i) != interpolant.iY_.at(i))
            return false;
    }
    for (unsigned int i = 0; i < c_.size(); i++)
    {
        if (c_.at(i) != interpolant.c_.at(i))
            return false;
    }
    for (unsigned int i = 0; i < d_.size(); i++)
    {
        if (d_.at(i) != interpolant.d_.at(i))
            return false;
    }
    for (unsigned int i = 0; i < interpolant.Interpolant_.size(); i++)
    {
        if (*Interpolant_.at(i) != *interpolant.Interpolant_.at(i))
            return false;
    }
    // else
    return true;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool Interpolant::operator!=(const Interpolant& interpolant) const
{
    return !(*this == interpolant);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Interpolant::swap(Interpolant& interpolant)
{

    using std::swap;

    swap(romberg_, interpolant.romberg_);
    swap(rombergY_, interpolant.rombergY_);
    swap(max_, interpolant.max_);
    swap(xmin_, interpolant.xmin_);
    swap(xmax_, interpolant.xmax_);
    swap(step_, interpolant.step_);
    swap(rational_, interpolant.rational_);
    swap(relative_, interpolant.relative_);
    swap(row_, interpolant.row_);
    swap(starti_, interpolant.starti_);
    swap(rationalY_, interpolant.rationalY_);
    swap(relativeY_, interpolant.relativeY_);
    swap(reverse_, interpolant.reverse_);
    swap(self_, interpolant.self_);
    swap(flag_, interpolant.flag_);
    swap(isLog_, interpolant.isLog_);
    swap(logSubst_, interpolant.logSubst_);
    swap(precision_, interpolant.precision_);
    swap(worstX_, interpolant.worstX_);
    swap(precision2_, interpolant.precision2_);
    swap(worstX2_, interpolant.worstX2_);
    swap(precisionY_, interpolant.precisionY_);
    swap(worstY_, interpolant.worstY_);
    swap(fast_, interpolant.fast_);
    swap(x_save_, interpolant.x_save_);
    swap(y_save_, interpolant.y_save_);

    iX_.swap(interpolant.iX_);
    iY_.swap(interpolant.iY_);

    c_.swap(interpolant.c_);
    d_.swap(interpolant.d_);

    Interpolant_.swap(interpolant.Interpolant_);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Interpolant::InitInterpolant(int max,
                                  double xmin,
                                  double xmax,
                                  int romberg,
                                  bool rational,
                                  bool relative,
                                  bool isLog,
                                  int rombergY,
                                  bool rationalY,
                                  bool relativeY,
                                  bool logSubst)
{

    self_   = true;
    fast_   = true;
    x_save_ = 1;
    y_save_ = 0;

    if (max <= 0)
    {
        log_warn("max = %i must be > 0! setting to 1!", max);
        max = 1;
    }

    if (isLog)
    {
        if (xmin <= 0)
        {
            log_warn("xmin = %f must be > 0! setting to 1!", xmin);
            xmin = 1;
        }

        if (xmax <= 0)
        {
            log_warn("xmax = %f must be > 0, setting to 1", xmax);
            xmax = 1;
        }
    }

    if (xmin == xmax)
    {
        max = 1;
    } else if (xmin > xmax)
    {
        std::swap(xmin, xmax);
    }

    if (romberg <= 0)
    {
        log_warn("romberg = %i must be > 0! setting to 1!", romberg);
        romberg = 1;
    }

    if (romberg > max)
    {
        log_warn("romberg = %i must be <= max = %i! setting to %i!", romberg, max, max);
        romberg = max;
    }

    if (rombergY <= 0)
    {
        log_warn("rombergY = %i must be > 0! setting to 1!", rombergY);
        rombergY = 1;
    }

    if (rombergY > max)
    {
        log_warn("rombergY = %i must be < %i! setting to %i!", rombergY, max, max);
        rombergY = max;
    }

    this->max_ = max;

    if ((std::log(max) / std::log(2) + romberg) < max)
    {
        flag_ = true;
    } else
    {
        flag_ = false;
    }

    if (isLog)
    {
        this->xmin_ = std::log(xmin);
        this->xmax_ = std::log(xmax);
    } else
    {
        this->xmin_ = xmin;
        this->xmax_ = xmax;
    }

    this->romberg_  = romberg;
    this->rombergY_ = rombergY;

    iX_.resize(max);
    iY_.resize(max);

    step_ = (this->xmax_ - this->xmin_) / max;

    int VectorMax_ = std::max(romberg, rombergY);

    c_.resize(VectorMax_);
    d_.resize(VectorMax_);

    precision_       = 0;
    this->isLog_     = isLog;
    this->logSubst_  = logSubst;
    this->rational_  = rational;
    this->relative_  = relative;
    this->rationalY_ = rationalY;
    this->relativeY_ = relativeY;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::Get2dFunctionFixedY(double x)
{
    if (isLog_)
    {
        return function2d_(x, std::exp(iX_.at(row_)));
    } else
    {
        return function2d_(x, iX_.at(row_));
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::Interpolate(double x, int start)
{
    int num, i, k;
    bool dd, doLog;
    double error = 0, result = 0;
    double aux, aux2, dx1, dx2;

    doLog = false;

    if (logSubst_)
    {
        if (reverse_)
        {
            for (i = 0; i < romberg_; i++)
            {
                if (iY_.at(start + i) == bigNumber_)
                {
                    doLog = true;
                    break;
                }
            }
        }
    }

    if (fast_)
    {
        num = starti_ - start;

        if (x == iX_.at(starti_))
        {
            return iY_.at(starti_);
        }

        if (doLog)
        {
            for (i = 0; i < romberg_; i++)
            {
                c_.at(i) = Exp(iY_.at(start + i));
                d_.at(i) = c_.at(i);
            }
        } else
        {
            for (i = 0; i < romberg_; i++)
            {
                c_.at(i) = iY_.at(start + i);
                d_.at(i) = c_.at(i);
            }
        }
    } else
    {
        num = 0;
        aux = std::abs(x - iX_.at(start + 0));

        for (i = 0; i < romberg_; i++)
        {
            aux2 = std::abs(x - iX_.at(start + i));

            if (aux2 == 0)
            {
                return iY_.at(start + i);
            }

            if (aux2 < aux)
            {
                num = i;
                aux = aux2;
            }

            if (doLog)
            {
                c_.at(i) = Exp(iY_.at(start + i));
                d_.at(i) = c_.at(i);
            } else
            {
                c_.at(i) = iY_.at(start + i);
                d_.at(i) = c_.at(i);
            }
        }
    }

    if (num == 0)
    {
        dd = true;
    } else if (num == romberg_ - 1)
    {
        dd = false;
    } else
    {
        k    = start + num;
        aux  = iX_.at(k - 1);
        aux2 = iX_.at(k + 1);

        if (fast_)
        {
            if (((x - aux) > (aux2 - x)) == (aux2 > aux))
            {
                dd = true;
            } else
            {
                dd = false;
            }
        } else
        {
            if (std::abs(x - aux) > std::abs(x - aux2))
            {
                dd = true;
            } else
            {
                dd = false;
            }
        }
    }

    result = iY_.at(start + num);

    if (doLog)
    {
        result = Exp(result);
    }

    for (k = 1; k < romberg_; k++)
    {
        for (i = 0; i < romberg_ - k; i++)
        {
            if (rational_)
            {
                aux  = c_.at(i + 1) - d_.at(i);
                dx2  = iX_.at(start + i + k) - x;
                dx1  = d_.at(i) * (iX_.at(start + i) - x) / dx2;
                aux2 = dx1 - c_.at(i + 1);

                if (aux2 != 0)
                {
                    aux      = aux / aux2;
                    d_.at(i) = c_.at(i + 1) * aux;
                    c_.at(i) = dx1 * aux;
                } else
                {
                    c_.at(i) = 0;
                    d_.at(i) = 0;
                }
            } else
            {
                dx1  = iX_.at(start + i) - x;
                dx2  = iX_.at(start + i + k) - x;
                aux  = c_.at(i + 1) - d_.at(i);
                aux2 = dx1 - dx2;

                if (aux2 != 0)
                {
                    aux      = aux / aux2;
                    c_.at(i) = dx1 * aux;
                    d_.at(i) = dx2 * aux;
                } else
                {
                    c_.at(i) = 0;
                    d_.at(i) = 0;
                }
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
            error = c_.at(num);
        } else
        {
            num--;
            error = d_.at(num);
        }

        dd = !dd;
        result += error;
    }

    if (!fast_)
    {
        if (relative_)
        {
            if (result != 0)
            {
                aux = std::abs(error / result);
            } else
            {
                aux = 0;
            }
        } else
        {
            aux = std::abs(error);
        }

        if (aux > precision_)
        {
            precision_ = aux;
            worstX_    = x;
        }
    }

    if (doLog)
    {
        result = Log(result);
    }

    return result;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::Exp(double x)
{
    if (x <= aBigNumber_)
    {
        return 0;
    } else
    {
        return std::exp(x);
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::Log(double x)
{
    if (x <= 0)
    {
        return bigNumber_;
    } else
    {
        return std::log(x);
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Interpolant::Slog(double x)
{
    if (x == x_save_)
    {
        return y_save_;
    } else
    {
        x_save_ = x;
        y_save_ = Log(x_save_);

        return y_save_;
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Interpolant::SetRombergY(int rombergY)
{
    rombergY_ = rombergY;
}

void Interpolant::SetRomberg(int romberg)
{
    romberg_ = romberg;
}

void Interpolant::SetIX(const std::vector<double>& iX)
{
    iX_ = iX;
}

void Interpolant::SetIY(const std::vector<double>& iY)
{
    iY_ = iY;
}

void Interpolant::SetC(const std::vector<double>& c)
{
    c_ = c;
}

void Interpolant::SetD(const std::vector<double>& d)
{
    d_ = d;
}

void Interpolant::SetMax(int max)
{
    max_ = max;
}

void Interpolant::SetXmin(double xmin)
{
    xmin_ = xmin;
}

void Interpolant::SetXmax(double xmax)
{
    xmax_ = xmax;
}

void Interpolant::SetStep(double step)
{
    step_ = step;
}

void Interpolant::SetRelative(bool relative)
{
    relative_ = relative;
}

void Interpolant::SetRational(bool rational)
{
    rational_ = rational;
}

void Interpolant::SetRow(int row)
{
    row_ = row;
}

void Interpolant::SetStarti(int starti)
{
    starti_ = starti;
}

void Interpolant::SetRationalY(bool rationalY)
{
    rationalY_ = rationalY;
}

void Interpolant::SetRelativeY(bool relativeY)
{
    relativeY_ = relativeY;
}

void Interpolant::SetSelf(bool self)
{
    self_ = self;
}

void Interpolant::SetFlag(bool flag)
{
    flag_ = flag;
}

void Interpolant::SetReverse(bool reverse)
{
    reverse_ = reverse;
}

void Interpolant::SetIsLog(bool isLog)
{
    isLog_ = isLog;
}

void Interpolant::SetLogSubst(bool logSubst)
{
    logSubst_ = logSubst;
}

void Interpolant::SetPrecision(double precision)
{
    precision_ = precision;
}

void Interpolant::SetWorstX(double worstX)
{
    worstX_ = worstX;
}

void Interpolant::SetWorstX2(double worstX2)
{
    worstX2_ = worstX2;
}

void Interpolant::SetWorstY(double worstY)
{
    worstY_ = worstY;
}

void Interpolant::SetPrecisionY(double precisionY)
{
    precisionY_ = precisionY;
}

void Interpolant::SetFast(bool fast)
{
    fast_ = fast;
}

void Interpolant::SetX_save(double x_save)
{
    x_save_ = x_save;
}

void Interpolant::SetY_save(double y_save)
{
    y_save_ = y_save;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Interpolant::~Interpolant()
{
    iX_.clear();
    iY_.clear();
    c_.clear();
    d_.clear();

    for (unsigned int i = 0; i < Interpolant_.size(); i++)
    {
        delete Interpolant_.at(i);
    }

    Interpolant_.clear();
}
