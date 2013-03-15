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

#include "PROPOSAL/Interpolant.h"
#include <algorithm>


using namespace std;

const double Interpolant::bigNumber   =   -300;
const double Interpolant::aBigNumber  =   -299;

//----------------------------------------------------------------------------//

Interpolant::Interpolant() { }

//----------------------------------------------------------------------------//

Interpolant::Interpolant(int max, double xmin, double xmax, boost::function<double (double)> function1d,
                         int romberg, bool rational, bool relative, bool isLog,
                         int rombergY, bool rationalY, bool relativeY, bool logSubst)
:self   (true)
,fast   (true)
,x_save (1)
,y_save (0)
{
    InitInterpolant(max, xmin, xmax, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);

    int     i;
    double  aux, xaux;


    function1d_     = function1d;

    aux     =   this->xmin+step/2;

    for(i=0; i<max; i++)
    {
        iX[i]   =   aux;

        if(Output::inf)
        {
            try
            {
                iY[i]   =   Output::read();
            }
            catch(int a)
            {
                throw 0;
            }
        }
        else
        {

            if(isLog)
            {
                xaux    =   std::exp(aux);
            }
            else
            {
                xaux    =   aux;
            }

            iY[i]=function1d_(xaux);

            if(logSubst)
            {
                iY[i]=log(iY[i]);
            }

        }

        if(output->outf)
        {
            output->write(iY[i]);
        }

        aux +=  step;
    }

}

//----------------------------------------------------------------------------//

Interpolant::Interpolant(int max1, double x1min, double x1max, int max2, double x2min, double x2max, boost::function<double (double,double)> function2d,
                         int romberg1, bool rational1, bool relative1, bool isLog1,
                         int romberg2, bool rational2, bool relative2, bool isLog2,
                         int rombergY, bool rationalY, bool relativeY, bool logSubst)
:self   (true)
,fast   (true)
,x_save (1)
,y_save (0)
{
    InitInterpolant(max2, x2min, x2max, romberg2, rational2, relative2, isLog2, rombergY, rationalY, relativeY, logSubst);

    int     i;
    double  aux;


    function2d_     = function2d;
    function1d_     = boost::bind(&Interpolant::Get2dFunctionFixedY, this, _1);

    Interpolant_ = (Interpolant *)calloc(max,sizeof(Interpolant));

    for(i=0,aux=xmin+step/2 ; i<max ; i++,aux+=step)
    {
        iX[i]   =   aux;
        row     =   i;

        Interpolant_[i]     =   Interpolant(max1, x1min, x1max, function1d_,
                                       romberg1, rational1, relative1, isLog1,
                                       rombergY, rationalY, relativeY, logSubst);

        Interpolant_[i].self=false;
    }

    precision2  =   0;
}

//----------------------------------------------------------------------------//

Interpolant::Interpolant(vector<double> x, vector<double> y, int romberg, bool rational, bool relative)
:self   (true)
,fast   (true)
,x_save (1)
,y_save (0)
{
    InitInterpolant(min(x.size(), y.size()), x.at(0), x.at(x.size()-1),
                    romberg, rational, relative, false,
                    romberg, rational, relative, false);

    if(x.size() != y.size())
    {
        cerr<<"Warning (in ): Interpolant/Interpolant: x and y do not match"<<endl;
    }

    for(int i=0 ; i<(int)x.size() ; i++)
    {
        iX[i]=x.at(i);
    }

    for(int i =0; i<(int)y.size(); i++)
    {
        iY[i]=y.at(i);
    }
}

//----------------------------------------------------------------------------//

void Interpolant::InitInterpolant(int max, double xmin, double xmax,
                                  int romberg, bool rational, bool relative, bool isLog,
                                  int rombergY, bool rationalY, bool relativeY, bool logSubst)
{
    double aux;

    self   =    true;
    fast   =    true;
    x_save =    1;
    y_save =    0;

    if(max<=0)
    {
        cerr<<"Warning (in Interpolant/Interpolant/0): max = "<<max<<" must be > 0, setting to 1"<<endl;
        max =   1;
    }

    if(isLog)
    {
        if(xmin<=0)
        {
            cerr<<"Warning (in Interpolant/Interpolant/1): xmin = "<<xmin<<" must be > 0, setting to 1"<<endl;
            xmin    =   1;
        }

        if(xmax<=0)
        {
            cerr<<"Warning (in Interpolant/Interpolant/2): xmax = "<<xmax<<" must be > 0, setting to 1"<<endl;
            xmax    =   1;
        }
    }

    if(xmin==xmax)
    {
        max     =   1;
    }
    else if(xmin>xmax)
    {
        aux     =   xmin;
        xmin    =   xmax;
        xmax    =   aux;
    }

    if(romberg<=0)
    {
        cerr<<"Warning (in Interpolant/Interpolant/3): romberg = "<<romberg<<" must be > 0, setting to 1"<<endl;
        cout<<"Warning (in Interpolant/Interpolant/3): romberg = "<<romberg<<" must be > 0, setting to 1"<<endl;
        romberg =   1;
    }

    if(romberg>max)
    {
        cerr<<"Warning (in Interpolant/Interpolant/4): romberg = "<<romberg<<" must be <= max = "<<max<<", setting to "<<max<<endl;
        cout<<"Warning (in Interpolant/Interpolant/4): romberg = "<<romberg<<" must be <= max = "<<max<<", setting to "<<max<<endl;
        romberg =   max;
    }

    if(rombergY<=0)
    {
        cerr<<"Warning (in Interpolant/Interpolant/5): rombergY = "<<rombergY<<" must be > 0, setting to 1"<<endl;
        cout<<"Warning (in Interpolant/Interpolant/5): rombergY = "<<rombergY<<" must be > 0, setting to 1"<<endl;
        rombergY    =   1;
    }

    if(rombergY>max)
    {
        cerr<<"Warning (in Interpolant/Interpolant/6): rombergY = "<<rombergY<<" must be <= max = "<<max<<", setting to "<<max<<endl;
        cout<<"Warning (in Interpolant/Interpolant/6): rombergY = "<<rombergY<<" must be <= max = "<<max<<", setting to "<<max<<endl;
        rombergY    =   max;
    }

    this->max   =   max;

    if((std::log(max)/std::log(2)+romberg)<max)
    {
        flag    =   true;
    }
    else
    {
        flag    =   false;
    }

    if(isLog)
    {
        this->xmin  =   std::log(xmin);
        this->xmax  =   std::log(xmax);
    }
    else
    {
        this->xmin  =   xmin;
        this->xmax  =   xmax;
    }

    this->romberg   =   romberg;
    this->rombergY  =   rombergY;


    iX  =   (double *)calloc(max,sizeof(double));
    iY  =   (double *)calloc(max,sizeof(double));

    step=   (this->xmax-this->xmin)/max;

    int VectorMax_ = std::max(romberg, rombergY);

    c   =   (double *)calloc(VectorMax_,sizeof(double));
    d   =   (double *)calloc(VectorMax_,sizeof(double));

    precision       =   0;
    this->isLog     =   isLog;
    this->logSubst  =   logSubst;
    this->rational  =   rational;
    this->relative  =   relative;
    this->rationalY =   rationalY;
    this->relativeY =   relativeY;
}

//----------------------------------------------------------------------------//

double Interpolant::Get2dFunctionFixedY(double x)
{
    if(isLog)
    {
        return function2d_(x, std::exp(iX[row]));
    }
    else
    {
        return function2d_(x, iX[row]);
    }
}

//----------------------------------------------------------------------------//

double Interpolant::interpolate(double x)
{
    int     start;
    double  result, aux;

    if(isLog)
    {
        x=slog(x);
    }

    reverse =   true;
    aux     =   (x - xmin)/step;
    starti  =   (int)aux;

    if(starti<0)
    {
        starti  =   0;
    }

    if(starti>=max)
    {
        starti  =   max-1;
    }

    start   =   (int)(aux - 0.5*(romberg-1));

    if(start<0)
    {
        start   =   0;
    }

    if(start+romberg>max || start>max)
    {
        start   =   max - romberg;
    }

    result  =   interpolate(x, start);

    if(logSubst)
    {
        if(self)
        {
            result  =   exp(result);
        }
    }

    return result;
}

//----------------------------------------------------------------------------//

double Interpolant::interpolate(double x1, double x2)
{
    int     i, start;
    double  aux, aux2=0, result;

    if(isLog)
    {
        x2  =   std::log(x2);
    }

    reverse =   true;

    aux     =   (x2 - xmin)/step;
    starti  =   (int)aux;

    if(starti<0)
    {
        starti  =   0;
    }

    if(starti>=max)
    {
        starti  =   max - 1;
    }

    start   =   (int)(aux-0.5*(romberg-1));

    if(start<0)
    {
        start   =   0;
    }

    if(start+romberg>max || start>max)
    {
        start   =   max - romberg;
    }

    for(i=start ; i<start+romberg ; i++)
    {
        iY[i]   =   Interpolant_[i].interpolate(x1);
    }

    if(!fast)
    {
        aux =   0;

        for(i=start ; i<start+romberg ; i++)
        {
            if(Interpolant_[i].precision>aux)
            {
                aux     =   Interpolant_[i].precision;
                aux2    =   Interpolant_[i].worstX;
            }
        }

        if(aux>precision2)
        {
            precision2  =   aux;
            worstX2     =   aux2;
        }
    }

    result  =   interpolate(x2, start);

    if(logSubst)
    {
        result  =   exp(result);
    }

    return result;
}

//----------------------------------------------------------------------------//

double Interpolant::interpolateArray(double x)
{
    int     i, j, m, start, auxdir;
    bool    dir;

    reverse =   false;
    i       =   0;
    j       =   max-1;
    dir     =   iX[max-1]>iX[0];

    while(j-i > 1)
    {
        m=(i+j)/2;

        if((x>iX[m])==dir)
        {
            i=m;
        }
        else
        {
            j=m;
        }
    }

    if(i+1 < max)
    {
        if(((x-iX[i])<(iX[i+1]-x))==dir)
        {
            auxdir=0;
        }
        else
        {
            auxdir=1;
        }
    }
    else
    {
        auxdir=0;
    }

    starti  =   i+auxdir;
    start   =   i-(int)(0.5*(romberg-1-auxdir));

    if(start<0)
    {
        start   =   0;
    }

    if(start+romberg>max || start>max)
    {
        start   =   max-romberg;
    }

    return interpolate(x, start);
}

//----------------------------------------------------------------------------//

double Interpolant::interpolate(double x, int start)
{
    int     num, i, k;
    bool    dd, doLog;
    double  error=0, result=0;
    double  aux, aux2, dx1, dx2;

    doLog   =   false;

    if(logSubst)
    {
        if(reverse)
        {
            for(i=0 ; i<romberg ; i++)
            {
                if(iY[start+i]==bigNumber)
                {
                    doLog   =   true;
                    break;
                }
            }
        }
    }

    if(fast)
    {
        num =   starti - start;

        if(x==iX[starti])
        {
            return iY[starti];
        }

        if(doLog)
        {
            for(i=0; i<romberg; i++)
            {
                c[i]    =   exp(iY[start+i]);
                d[i]    =   c[i];
            }
        }
        else
        {
            for(i=0; i<romberg; i++)
            {
                c[i]    =   iY[start+i];
                d[i]    =   c[i];
            }
        }
    }
    else
    {
        num =   0;
        aux =   fabs(x-iX[start+0]);

        for(i=0 ; i<romberg ; i++)
        {
            aux2    =   fabs(x-iX[start+i]);

            if(aux2==0)
            {
                return iY[start+i];
            }

            if(aux2<aux)
            {
                num =   i;
                aux =   aux2;
            }

            if(doLog)
            {
                c[i]    =   exp(iY[start+i]);
                d[i]    =   c[i];
            }
            else
            {
                c[i]    =   iY[start+i];
                d[i]    =   c[i];
            }
        }
    }

    if(num==0)
    {
        dd  =   true;
    }
    else if(num==romberg-1)
    {
        dd  =   false;
    }
    else
    {
        k       =   start+num;
        aux     =   iX[k-1];
        aux2    =   iX[k+1];

        if(fast)
        {
            if(((x-aux)>(aux2-x)) == (aux2>aux))
            {
                dd  =   true;
            }
            else
            {
                dd  =   false;
            }
        }
        else
        {
            if(fabs(x-aux) > fabs(x-aux2))
            {
                dd  =   true;
            }
            else
            {
                dd  =   false;
            }
        }
    }

    result  =   iY[start+num];

    if(doLog)
    {
        result  =   exp(result);
    }

    for(k=1 ; k<romberg ; k++)
    {
        for(i=0 ; i<romberg-k ; i++)
        {
            if(rational)
            {
                aux     =   c[i+1] - d[i];
                dx2     =   iX[start+i+k] - x;
                dx1     =   d[i]*(iX[start+i] - x)/dx2;
                aux2    =   dx1 - c[i+1];

                if(aux2!=0)
                {
                    aux     =   aux/aux2;
                    d[i]    =   c[i+1]*aux;
                    c[i]    =   dx1*aux;
                }
                else
                {
                    c[i]    =   0;
                    d[i]    =   0;
                }
            }
            else
            {
                dx1     =   iX[start+i] - x;
                dx2     =   iX[start+i+k] - x;
                aux     =   c[i+1] - d[i];
                aux2    =   dx1 - dx2;

                if(aux2!=0)
                {
                    aux     =   aux/aux2;
                    c[i]    =   dx1*aux;
                    d[i]    =   dx2*aux;
                }
                else
                {
                    c[i]    =   0;
                    d[i]    =   0;
                }
            }
        }

        if(num==0)
        {
            dd  =   true;
        }

        if(num==romberg-k)
        {
            dd  =   false;
        }

        if(dd)
        {
            error   =   c[num];
        }
        else
        {
            num--;
            error   =   d[num];
        }

        dd      =   !dd;
        result  +=  error;
    }

    if(!fast)
    {
        if(relative)
        {
            if(result!=0)
            {
                aux =   fabs(error/result);
            }
            else
            {
                aux =   0;
            }
        }
        else
        {
            aux =   fabs(error);
        }

        if(aux>precision)
        {
            precision   =   aux;
            worstX      =   x;
        }
    }

    if(doLog)
    {
        result  =   log(result);
    }

    return result;
}

//----------------------------------------------------------------------------//

double Interpolant::findLimit(double y)
{
    int         i, j, m, start, auxdir, auxR;
    bool        dir, rat, rel;
    double      result, aux;
    double*     ii;

    rat     =   false;
    reverse =   false;

    if(logSubst)
    {
        y   =   log(y);
    }

    i   =   0;
    j   =   max - 1;
    dir =   iY[max-1]>iY[0];

    while(j-i > 1)
    {
        m   =   (i+j)/2;

        if((y>iY[m])==dir)
        {
            i   =   m;
        }
        else
        {
            j   =   m;
        }
    }

    ii  =   iX;
    iX  =   iY;
    iY  =   ii;

    if(!fast)
    {
        aux         =   precision;
        precision   =   precisionY;
        precisionY  =   aux;

        aux         =   worstX;
        worstX      =   worstY;
        worstY      =   aux;

        rat         =   rational;
        rational    =   rationalY;
    }

    auxR        =   romberg;
    romberg     =   rombergY;
    rombergY    =   auxR;
    rel         =   relative;
    relative    =   relativeY;

    if(i+1 < max)
    {
        if(((y-iX[i])<(iX[i+1]-y))==dir)
        {
            auxdir  =   0;
        }
        else
        {
            auxdir  =   1;
        }
    }
    else
    {
        auxdir  =   0;
    }

    starti  =   i + auxdir;
    start   =   i - (int)(0.5*(romberg-1-auxdir));

    if(start<0)
    {
        start   =   0;
    }

    if(start+romberg>max || start>max)
    {
        start   =   max - romberg;
    }

    result  =   interpolate(y, start);
    ii      =   iX;
    iX      =   iY;
    iY      =   ii;

    if(!fast)
    {
        aux         =   precision;
        precision   =   precisionY;
        precisionY  =   aux;

        aux         =   worstX;
        worstX      =   worstY;
        worstY      =   aux;

        rational    =   rat;
    }

    auxR        =   romberg;
    romberg     =   rombergY;
    rombergY    =   auxR;
    relative    =   rel;

    if(result<xmin)
    {
        result  =   xmin;
    }
    else if(result>xmax)
    {
        result  =   xmax;
    }

    if(isLog)
    {
        result  =   exp(result);
    }

    return result;
}

//----------------------------------------------------------------------------//

double Interpolant::findLimit(double x1, double y)
{
    int         i, j, m, start, auxdir, auxR;
    bool        dir, rat, rel;
    double      result, aux, aux2=0;
    double*     ii;
    rat     =   false;
    reverse =   false;

    if(logSubst)
    {
        y   =   log(y);
    }

    if(!flag)
    {
        for(i=0 ; i<max ; i++)
        {
            iY[i]   =   Interpolant_[i].interpolate(x1);
        }
    }

    i   =   0;
    j   =   max-1;

    if(flag)
    {
        dir =   Interpolant_[max-1].interpolate(x1)>Interpolant_[0].interpolate(x1);
    }
    else
    {
        dir =   iY[max-1]>iY[0];
    }

    while(j-i > 1)
    {
        m=(i+j)/2;

        if(flag)
        {
            aux =   Interpolant_[m].interpolate(x1);
        }
        else
        {
            aux =   iY[m];
        }

        if((y>aux) == dir)
        {
            i   =   m;
        }
        else
        {
            j   =   m;
        }
    }

    ii  =   iX;
    iX  =   iY;
    iY  =   ii;

    if(!fast)
    {
        aux         =   precision;
        precision   =   precisionY;
        precisionY  =   aux;

        aux         =   worstX;
        worstX      =   worstY;
        worstY      =   aux;

        rat         =   rational;
        rational    =   rationalY;
    }

    auxR        =   romberg;
    romberg     =   rombergY;
    rombergY    =   auxR;
    rel         =   relative;
    relative    =   relativeY;

    if(i+1<max)
    {
        if(((y-iX[i])<(iX[i+1]-y))==dir)
        {
            auxdir  =   0;
        }
        else
        {
            auxdir  =   1;
        }
    }
    else
    {
        auxdir  =   0;
    }

    starti  =   i + auxdir;
    start   =   i - (int)(0.5*(romberg-1-auxdir));

    if(start<0)
    {
        start   =   0;
    }

    if(start+romberg>max || start>max)
    {
        start   =   max - romberg;
    }

    if(flag)
    {
        for(i=start ; i<start+romberg ; i++)
        {
            iX[i]   =   Interpolant_[i].interpolate(x1);
        }
    }

    result  =   interpolate(y, start);
    ii      =   iX;
    iX      =   iY;
    iY      =   ii;

    if(!fast)
    {
        aux         =   precision;
        precision   =   precisionY;
        precisionY  =   aux;

        aux         =   worstX;
        worstX      =   worstY;
        worstY      =   aux;

        rational    =   rat;
    }

    auxR        =   romberg;
    romberg     =   rombergY;
    rombergY    =   auxR;
    relative    =   rel;

    if(result<xmin)
    {
        result  =   xmin;
    }
    else if(result>xmax)
    {
        result  =   xmax;
    }

    if(!fast)
    {
        aux=0;

        for(i=start ; i<start+romberg ; i++)
        {
            if(Interpolant_[i].precision>aux)
            {
                aux     =   Interpolant_[i].precision;
                aux2    =   Interpolant_[i].worstX;
            }
        }

        if(aux>precision2)
        {
            precision2  =   aux;
            worstX2     =   aux2;
        }
    }

    if(isLog)
    {
        result  =   std::exp(result);
    }

    return result;
}

//----------------------------------------------------------------------------//

double Interpolant::exp(double x)
{
    if(x<=aBigNumber)
    {
        return 0;
    }
    else
    {
        return std::exp(x);
    }
}

//----------------------------------------------------------------------------//

double Interpolant::log(double x)
{
    if(x<=0)
    {
        return bigNumber;
    }
    else
    {
        return std::log(x);
    }
}

//----------------------------------------------------------------------------//

double Interpolant::slog(double x)
{
    if(x == x_save)
    {
        return y_save;
    }
    else
    {
        x_save  =   x;
        y_save  =   log(x_save);

        return y_save;
    }
}
