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
#include "PROPOSAL/methods.h"
#include <algorithm>


using namespace std;

const double Interpolant::bigNumber_   =   -300;
const double Interpolant::aBigNumber_  =   -299;

//----------------------------------------------------------------------------//

Interpolant::Interpolant()
:iX_()
,iY_()
,c_()
,d_()
,function1d_(NULL)
,function2d_(NULL)
,self_   (true)
,fast_   (true)
,x_save_ (1)
,y_save_ (0) { }

//----------------------------------------------------------------------------//

Interpolant::Interpolant(int max, double xmin, double xmax, boost::function<double (double)> function1d,
                         int romberg, bool rational, bool relative, bool isLog,
                         int rombergY, bool rationalY, bool relativeY, bool logSubst)
:iX_()
,iY_()
,c_()
,d_()
,function1d_(NULL)
,function2d_(NULL)
,self_   (true)
,fast_   (true)
,x_save_ (1)
,y_save_ (0)
{
    InitInterpolant(max, xmin, xmax, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);

    int     i;
    double  aux, xaux;


    function1d_     = function1d;

    aux     =   this->xmin_+step_/2;

    for(i=0; i<max_; i++)
    {
        iX_.at(i)   =   aux;

        if(Output::inf)
        {
            try
            {
                iY_.at(i)   =   Output::read();
            }
            catch(int a)
            {
                throw 0;
            }
        }
        else
        {

            if(isLog_)
            {
                xaux    =   std::exp(aux);
            }
            else
            {
                xaux    =   aux;
            }

            iY_.at(i)=function1d_(xaux);

            if(logSubst_)
            {
                iY_.at(i)=log(iY_.at(i));
            }

        }

        if(output_->outf)
        {
            output_->write(iY_.at(i));
        }

        aux +=  step_;
    }

}

//----------------------------------------------------------------------------//

Interpolant::Interpolant(int max1, double x1min, double x1max, int max2, double x2min, double x2max, boost::function<double (double,double)> function2d,
                         int romberg1, bool rational1, bool relative1, bool isLog1,
                         int romberg2, bool rational2, bool relative2, bool isLog2,
                         int rombergY, bool rationalY, bool relativeY, bool logSubst)
:iX_()
,iY_()
,c_()
,d_()
,function1d_(NULL)
,function2d_(NULL)
,self_   (true)
,fast_   (true)
,x_save_ (1)
,y_save_ (0)
{
    InitInterpolant(max2, x2min, x2max, romberg2, rational2, relative2, isLog2, rombergY, rationalY, relativeY, logSubst);

    int     i;
    double  aux;


    function2d_     = function2d;
    function1d_     = boost::bind(&Interpolant::Get2dFunctionFixedY, this, _1);

    Interpolant_.resize(max_);

    for(i=0,aux=xmin_+step_/2 ; i<max_ ; i++,aux+=step_)
    {
        iX_.at(i)   =   aux;
        row_     =   i;

        Interpolant_.at(i)     =  new Interpolant(max1, x1min, x1max, function1d_,
                                       romberg1, rational1, relative1, isLog1,
                                       rombergY, rationalY, relativeY, logSubst_);

        Interpolant_.at(i)->self_=false;
    }

    precision2_  =   0;
}

//----------------------------------------------------------------------------//

Interpolant::Interpolant(vector<double> x, vector<double> y, int romberg, bool rational, bool relative)
:iX_()
,iY_()
,c_()
,d_()
,function1d_(NULL)
,function2d_(NULL)
,self_   (true)
,fast_   (true)
,x_save_ (1)
,y_save_ (0)
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
        iX_.at(i)=x.at(i);
    }

    for(int i =0; i<(int)y.size(); i++)
    {
        iY_.at(i)=y.at(i);
    }
}

//----------------------------------------------------------------------------//

void Interpolant::InitInterpolant(int max, double xmin, double xmax,
                                  int romberg, bool rational, bool relative, bool isLog,
                                  int rombergY, bool rationalY, bool relativeY, bool logSubst)
{
    //double aux;

    self_   =    true;
    fast_   =    true;
    x_save_ =    1;
    y_save_ =    0;

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
        SWAP(xmin,xmax,double);
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

    this->max_   =   max;

    if((std::log(max)/std::log(2)+romberg)<max)
    {
        flag_    =   true;
    }
    else
    {
        flag_    =   false;
    }

    if(isLog)
    {
        this->xmin_  =   std::log(xmin);
        this->xmax_  =   std::log(xmax);
    }
    else
    {
        this->xmin_  =   xmin;
        this->xmax_  =   xmax;
    }

    this->romberg_   =   romberg;
    this->rombergY_  =   rombergY;


    iX_.resize(max);
    iY_.resize(max);

    step_=   (this->xmax_-this->xmin_)/max;

    int VectorMax_ = std::max(romberg, rombergY);

    c_.resize(VectorMax_);
    d_.resize(VectorMax_);

    precision_       =   0;
    this->isLog_     =   isLog;
    this->logSubst_  =   logSubst;
    this->rational_  =   rational;
    this->relative_  =   relative;
    this->rationalY_ =   rationalY;
    this->relativeY_ =   relativeY;
}

//----------------------------------------------------------------------------//

double Interpolant::Get2dFunctionFixedY(double x)
{
    if(isLog_)
    {
        return function2d_(x, std::exp(iX_.at(row_)));
    }
    else
    {
        return function2d_(x, iX_.at(row_));
    }
}

//----------------------------------------------------------------------------//

double Interpolant::interpolate(double x)
{
    int     start;
    double  result, aux;

    if(isLog_)
    {
        x=slog(x);
    }

    reverse_ =   true;
    aux     =   (x - xmin_)/step_;
    starti_  =   (int)aux;

    if(starti_<0)
    {
        starti_  =   0;
    }

    if(starti_>=max_)
    {
        starti_  =   max_-1;
    }

    start   =   (int)(aux - 0.5*(romberg_-1));

    if(start<0)
    {
        start   =   0;
    }

    if(start+romberg_>max_ || start>max_)
    {
        start   =   max_ - romberg_;
    }
    result  =   interpolate(x, start);

    if(logSubst_)
    {
        if(self_)
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

    if(isLog_)
    {
        x2  =   std::log(x2);
    }

    reverse_ =   true;

    aux     =   (x2 - xmin_)/step_;
    starti_  =   (int)aux;

    if(starti_<0)
    {
        starti_  =   0;
    }

    if(starti_>=max_)
    {
        starti_  =   max_ - 1;
    }

    start   =   (int)(aux-0.5*(romberg_-1));

    if(start<0)
    {
        start   =   0;
    }

    if(start+romberg_>max_ || start>max_)
    {
        start   =   max_ - romberg_;
    }

    for(i=start ; i<start+romberg_ ; i++)
    {
        iY_.at(i)   =   Interpolant_.at(i)->interpolate(x1);
    }

    if(!fast_)
    {
        aux =   0;

        for(i=start ; i<start+romberg_ ; i++)
        {
            if(Interpolant_.at(i)->precision_>aux)
            {
                aux     =   Interpolant_.at(i)->precision_;
                aux2    =   Interpolant_.at(i)->worstX_;
            }
        }

        if(aux>precision2_)
        {
            precision2_  =   aux;
            worstX2_     =   aux2;
        }
    }

    result  =   interpolate(x2, start);

    if(logSubst_)
    {
        result  =   exp(result);
    }

    return result;
}

//----------------------------------------------------------------------------//

double Interpolant::InterpolateArray(double x)
{
    int     i, j, m, start, auxdir;
    bool    dir;

    reverse_ =   false;
    i       =   0;
    j       =   max_-1;
    dir     =   iX_.at(max_-1)>iX_.at(0);

    while(j-i > 1)
    {
        m=(i+j)/2;

        if((x>iX_.at(m))==dir)
        {
            i=m;
        }
        else
        {
            j=m;
        }
    }

    if(i+1 < max_)
    {
        if(((x-iX_.at(i))<(iX_.at(i+1)-x))==dir)
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

    starti_  =   i+auxdir;
    start   =   i-(int)(0.5*(romberg_-1-auxdir));

    if(start<0)
    {
        start   =   0;
    }

    if(start+romberg_>max_ || start>max_)
    {
        start   =   max_-romberg_;
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

    if(logSubst_)
    {
        if(reverse_)
        {
            for(i=0 ; i<romberg_ ; i++)
            {
                if(iY_.at(start+i)==bigNumber_)
                {
                    doLog   =   true;
                    break;
                }
            }
        }
    }

    if(fast_)
    {
        num =   starti_ - start;

        if(x==iX_.at(starti_))
        {
            return iY_.at(starti_);
        }

        if(doLog)
        {
            for(i=0; i<romberg_; i++)
            {
                c_.at(i)    =   exp(iY_.at(start+i));
                d_.at(i)    =   c_.at(i);
            }
        }
        else
        {
            for(i=0; i<romberg_; i++)
            {
                c_.at(i)    =   iY_.at(start+i);
                d_.at(i)    =   c_.at(i);
            }
        }
    }
    else
    {
        num =   0;
        aux =   fabs(x-iX_.at(start+0));

        for(i=0 ; i<romberg_ ; i++)
        {
            aux2    =   fabs(x-iX_.at(start+i));

            if(aux2==0)
            {
                return iY_.at(start+i);
            }

            if(aux2<aux)
            {
                num =   i;
                aux =   aux2;
            }

            if(doLog)
            {
                c_.at(i)    =   exp(iY_.at(start+i));
                d_.at(i)    =   c_.at(i);
            }
            else
            {
                c_.at(i)    =   iY_.at(start+i);
                d_.at(i)    =   c_.at(i);
            }
        }
    }

    if(num==0)
    {
        dd  =   true;
    }
    else if(num==romberg_-1)
    {
        dd  =   false;
    }
    else
    {
        k       =   start+num;
        aux     =   iX_.at(k-1);
        aux2    =   iX_.at(k+1);

        if(fast_)
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

    result  =   iY_.at(start+num);

    if(doLog)
    {
        result  =   exp(result);
    }

    for(k=1 ; k<romberg_ ; k++)
    {
        for(i=0 ; i<romberg_-k ; i++)
        {
            if(rational_)
            {
                aux     =   c_.at(i+1) - d_.at(i);
                dx2     =   iX_.at(start+i+k) - x;
                dx1     =   d_.at(i)*(iX_.at(start+i) - x)/dx2;
                aux2    =   dx1 - c_.at(i+1);

                if(aux2!=0)
                {
                    aux     =   aux/aux2;
                    d_.at(i)    =   c_.at(i+1)*aux;
                    c_.at(i)    =   dx1*aux;
                }
                else
                {
                    c_.at(i)    =   0;
                    d_.at(i)    =   0;
                }
            }
            else
            {
                dx1     =   iX_.at(start+i) - x;
                dx2     =   iX_.at(start+i+k) - x;
                aux     =   c_.at(i+1) - d_.at(i);
                aux2    =   dx1 - dx2;

                if(aux2!=0)
                {
                    aux     =   aux/aux2;
                    c_.at(i)    =   dx1*aux;
                    d_.at(i)    =   dx2*aux;
                }
                else
                {
                    c_.at(i)    =   0;
                    d_.at(i)    =   0;
                }
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
            error   =   c_.at(num);
        }
        else
        {
            num--;
            error   =   d_.at(num);
        }

        dd      =   !dd;
        result  +=  error;
    }

    if(!fast_)
    {
        if(relative_)
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

        if(aux>precision_)
        {
            precision_   =   aux;
            worstX_      =   x;
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
    int         i, j, m, start, auxdir;
    bool        dir, rat, rel;
    double      result;

    rat     =   false;
    reverse_ =   false;

    if(logSubst_)
    {
        y   =   log(y);
    }

    i   =   0;
    j   =   max_ - 1;
    dir =   iY_.at(max_-1)>iY_.at(0);

    while(j-i > 1)
    {
        m   =   (i+j)/2;

        if((y>iY_.at(m))==dir)
        {
            i   =   m;
        }
        else
        {
            j   =   m;
        }
    }

    SWAP(iX_,iY_, vector<double>);
    if(!fast_)
    {
        SWAP(precision_,precisionY_,double);
        SWAP(worstX_,worstY_,double);

        rat         =   rational_;
        rational_    =   rationalY_;
    }

    SWAP(romberg_,rombergY_,double);

    rel         =   relative_;
    relative_    =   relativeY_;

    if(i+1 < max_)
    {
        if(((y-iX_.at(i))<(iX_.at(i+1)-y))==dir)
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

    starti_  =   i + auxdir;
    start   =   i - (int)(0.5*(romberg_-1-auxdir));

    if(start<0)
    {
        start   =   0;
    }

    if(start+romberg_>max_ || start>max_)
    {
        start   =   max_ - romberg_;
    }

    result  =   interpolate(y, start);

    SWAP(iX_,iY_, vector<double>);

    if(!fast_)
    {
        SWAP(precision_,precisionY_,double);
        SWAP(worstX_,worstY_,double);

        rational_    =   rat;
    }

    SWAP(romberg_,rombergY_,double);

    relative_    =   rel;

    if(result<xmin_)
    {
        result  =   xmin_;
    }
    else if(result>xmax_)
    {
        result  =   xmax_;
    }

    if(isLog_)
    {
        result  =   exp(result);
    }

    return result;
}

//----------------------------------------------------------------------------//

double Interpolant::findLimit(double x1, double y)
{
    int         i, j, m, start, auxdir;
    bool        dir;
    double      result, aux, aux2=0;

    reverse_ =   false;

    if(logSubst_)
    {
        y   =   log(y);
    }

    if(!flag_)
    {
        for(i=0 ; i<max_ ; i++)
        {
            iY_.at(i)   =   Interpolant_.at(i)->interpolate(x1);
        }
    }

    i   =   0;
    j   =   max_-1;

    if(flag_)
    {
        dir =   Interpolant_.at(max_-1)->interpolate(x1) > Interpolant_.at(0)->interpolate(x1);
    }
    else
    {
        dir =   iY_.at(max_-1)>iY_.at(0);
    }

    while(j-i > 1)
    {
        m=(i+j)/2;

        if(flag_)
        {
            aux =   Interpolant_.at(m)->interpolate(x1);
        }
        else
        {
            aux =   iY_.at(m);
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

    SWAP(iX_,iY_, vector<double>);

    if(!fast_)
    {
        SWAP(precision_, precisionY_, double);
        SWAP(worstX_,worstY_,double);
        SWAP(rational_,rationalY_,double);
    }

    SWAP(romberg_, rombergY_,double);
    SWAP(relative_,relativeY_,double);

    if(i+1<max_)
    {
        if(((y-iX_.at(i))<(iX_.at(i+1)-y))==dir)
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

    starti_  =   i + auxdir;
    start   =   i - (int)(0.5*(romberg_-1-auxdir));

    if(start<0)
    {
        start   =   0;
    }

    if(start+romberg_>max_ || start>max_)
    {
        start   =   max_ - romberg_;
    }

    if(flag_)
    {
        for(i=start ; i<start+romberg_ ; i++)
        {
            iX_.at(i)   =   Interpolant_.at(i)->interpolate(x1);
        }
    }

    result  =   interpolate(y, start);

    SWAP(iX_,iY_, vector<double>);

    if(!fast_)
    {
        SWAP(precision_, precisionY_, double);
        SWAP(worstX_,worstY_,double);
        SWAP(rational_,rationalY_,double);
    }

    SWAP(romberg_, rombergY_,double);
    SWAP(relative_,relativeY_,double);

    if(result<xmin_)
    {
        result  =   xmin_;
    }
    else if(result>xmax_)
    {
        result  =   xmax_;
    }

    if(!fast_)
    {
        aux=0;

        for(i=start ; i<start+romberg_ ; i++)
        {
            if(Interpolant_.at(i)->precision_>aux)
            {
                aux     =   Interpolant_.at(i)->precision_;
                aux2    =   Interpolant_.at(i)->worstX_;
            }
        }

        if(aux>precision2_)
        {
            precision2_  =   aux;
            worstX2_     =   aux2;
        }
    }

    if(isLog_)
    {
        result  =   std::exp(result);
    }

    return result;
}

bool Interpolant::Save(string Path){
    ofstream out;
    out.open(Path.c_str());
    out.precision(16);

    Save(out);

    out.close();
    return 1;
}

bool Interpolant::Save(ofstream& out){
    bool D2 = false;
    if(function2d_ != NULL)D2 = true;

    out << D2 << endl;

    if(D2)
    {
        out << max_ << "\t" << xmin_ << "\t" << xmax_ << endl;
        out << romberg_ << "\t" << rational_ << "\t" << relative_ << "\t" << isLog_ << endl;
        out << rombergY_ << "\t" << rationalY_ << "\t" << relativeY_ << "\t" << logSubst_ << endl;

        for(int i =0; i<max_ ;i++){
            out << iX_.at(i) << endl;
            Interpolant_.at(i)->Save(out);
        }
    }
    else
    {
        out << max_ << "\t" << xmin_ << "\t" << xmax_ << endl;
        out << romberg_ << "\t" << rational_ << "\t" << relative_ << "\t" << isLog_ << endl;
        out << rombergY_ << "\t" << rationalY_ << "\t" << relativeY_ << "\t" << logSubst_ << endl;

        for(int i =0; i<max_ ;i++){
            out << iX_.at(i) << "\t" << iY_.at(i) << endl;
        }
    }
    return 1;
}

bool Interpolant::Load(std::string Path){
    ifstream in;
    in.open(Path.c_str());

    Load(in);

    in.close();
    return 1;
}

bool Interpolant::Load(ifstream& in){
    bool D2;
    in >> D2 ;

    int max;
    double xmin, xmax;
    int romberg, rombergY;
    bool rational,rationalY,relative,relativeY,isLog, logSubst;

    if(D2)
    {
        in >> max >> xmin >> xmax;
        in >> romberg >> rational >> relative >> isLog ;
        in >> rombergY >> rationalY >> relativeY >> logSubst ;

        InitInterpolant(max, xmin, xmax,
                        romberg, rational, relative, isLog,
                        rombergY,rationalY, relativeY, logSubst);

        Interpolant_.resize(max_);

        for(int i =0; i<max_ ;i++){
            in >> iX_.at(i);

            Interpolant_.at(i) = new Interpolant();
            Interpolant_.at(i)->Load(in);
            Interpolant_.at(i)->self_ =false;
        }

    }
    else
    {
        in >> max >> xmin >> xmax;
        in >> romberg >> rational >> relative >> isLog ;
        in >> rombergY >> rationalY >> relativeY >> logSubst ;

        InitInterpolant(max, xmin, xmax,
                        romberg, rational, relative, isLog,
                        rombergY,rationalY, relativeY, logSubst);

        for(int i =0; i<max_ ;i++){
            in >> iX_.at(i) >> iY_.at(i) ;
        }
    }
    return 1;
}

//----------------------------------------------------------------------------//

double Interpolant::exp(double x)
{
    if(x<=aBigNumber_)
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
        return bigNumber_;
    }
    else
    {
        return std::log(x);
    }
}

//----------------------------------------------------------------------------//

double Interpolant::slog(double x)
{
    if(x == x_save_)
    {
        return y_save_;
    }
    else
    {
        x_save_  =   x;
        y_save_  =   log(x_save_);

        return y_save_;
    }
}

Interpolant::~Interpolant(){
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
