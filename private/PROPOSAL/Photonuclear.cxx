/*! \file   PhotoStochastic.cxx
*   \brief  Source file for the photonuclear routines.
*
*   For more details see the class documentation.
*
*   \date   29.06.2010
*   \author Jan-Hendrik Koehne
*/


#include "PROPOSAL/Photonuclear.h"
#include <algorithm>
#include <cmath>
#include <string>
#include "PROPOSAL/PhotoStochastic.h"
#include "PROPOSAL/PhotoContinuous.h"
#include "PROPOSAL/Medium.h"

using namespace std;

Photonuclear::Photonuclear(Photonuclear *cros_)
{

    bb              =   1;
    form            =   1;
    initM           =   true;
    shadow          =   1;
    jt_             =   false;
    initH           =   true;
    this->particle_ =   cros_->particle_;
    this->medium_   =   cros_->medium_;
    this->cros      =   cros_->cros;
    this->photo     =   cros_;
}

//----------------------------------------------------------------------------//

Photonuclear::Photonuclear(CrossSections *cros)
:CrossSections(*cros)
{

    bb          =   1;
    form        =   1;
    initM       =   true;
    shadow      =   1;
    jt_         =   false;
    initH       =   true;
    photo       =   this;
    continuous_ =   new PhotoContinuous(this);
    stochastic_ =   new PhotoStochastic(this);
    integral_   =   new Integral(IROMB, IMAXS, IPREC);
}

//----------------------------------------------------------------------------//

void Photonuclear::setEnergy(int i)
{
    double aux;

    cros->set_component(i);
    vMin    =   (MPI + (MPI*MPI)/(2*medium_->get_M().at(i)))/particle_->e;

    if(particle_->m<MPI)
    {
        aux     =   particle_->m/medium_->get_M().at(i);
        vMax    =   1 - medium_->get_M().at(i)*(1 + aux*aux)/(2*particle_->e);
    }
    else
    {
        vMax    =   1;
    }

    vMax    =   min(vMax, 1-particle_->m/particle_->e);

    if(vMax<vMin)
    {
        vMax    =   vMin;
    }

    vUp     =   min(vMax, medium_->vCut(particle_->e));

    if(vUp<vMin)
    {
        vUp =   vMin;
    }

    if(photo!=this)
    {
        photo->vMin =   vMin;
    }
}

//----------------------------------------------------------------------------//

double Photonuclear::photoN(double v, int i)
{
    switch(form)
    {
    case 1:
    case 2:
    {

        double nu, sgn, aux, aum, k, G, t;

        const double m1 =   0.54;
        const double m2 =   1.80;
        nu              =   v*particle_->e;
        nu              *=  1.e-3;


        switch(bb)
        {
        case 1:
        case 2:

            if(nu<=200.)
            {
                if(bb==2)
                {

                    sgn =   measuredSgN(nu);

                }
                else
                {

                    if(nu<=17.)
                    {
                        sgn =   96.1+82./sqrt(nu);
                    }
                    else
                    {
                        aux =   log(0.0213*nu);
                        sgn =   114.3 + 1.647*aux*aux;
                    }
                }
            }
            else
            {
                sgn =   49.2 + 11.1*log(nu) + 151.8/sqrt(nu);
            }

            break;

        case 3:  // Bezrukov and Bugaev (BB) parametrization
            aux =   log(0.0213*nu);
            sgn =   114.3 + 1.647*aux*aux;
            break;

        case 4:  // ZEUS parametrization

        default:
            aux =   nu*2.e-3*medium_->get_M().at(i);
            sgn =   63.5*pow(aux, 0.097) + 145*pow(aux , -0.5);
        }

        k   =   1 - 2/v + 2/(v*v);

        if(medium_->get_NucCharge().at(i)==1)
        {
            G   =   1;
        }
        else
        {
            aux =   0.00282*pow(medium_->get_atomicNum().at(i), 1./3)*sgn;
            G   =   (3/aux)*(0.5 + ((1 + aux)*exp(-aux) - 1)/(aux*aux));
        }

        G       *=  3;
        aux     =   v*particle_->m*1.e-3;
        t       =   aux*aux/(1-v);
        sgn     *=  1.e-30;
        aum     =   particle_->m*1.e-3;
        aum     *=  aum;
        aux     =   2*aum/t;
        aux     =   G*((k + 4*aum/m1)*log(1 + m1/t) - (k*m1)/(m1 + t) - aux) +
                    ((k + 2*aum/m2)*log(1 + m2/t) - aux) + 0.25*aux*((G*(m1 - 4*t))/(m1 + t) + (m2/t)*log(1 + t/m2));

        aux     *=  ALPHA/(8*PI)*medium_->get_atomicNum().at(i)*sgn*v;

        if(form==2)
        {
            if(particle_->type==1 || particle_->type==2)
            {
                aux +=  medium_->get_atomicNum().at(i)*1.e-30*hardBB(particle_->e, v);
            }
        }

        return medium_->get_molDensity()*medium_->get_atomInMolecule().at(i)*particle_->c*particle_->c*aux;
    }
    case 3:
    case 4:
    default:
    {
        if(jt_)
        {
            setEnergy(i);

            if(v>=vUp)
            {
                return max(interpolateJ_[i].interpolate(particle_->e, log(v/vUp)/log(vMax/vUp)), 0.0);
            }
        }

        double aux, min, max;

        this->i =   i;
        this->v =   v;
        min     =   particle_->m*v;
        min     *=  min/(1-v);

        if(particle_->m < MPI)
        {
            aux     =   particle_->m*particle_->m/particle_->e;
            min     -=  (aux*aux)/(2*(1-v));
        }

        max =   2*medium_->get_M().at(i)*particle_->e*(v-vMin);

        //  if(form==4) max=Math.min(max, 5.5e6);  // as requested in Butkevich and Mikheyev
        if(min>max)
        {
            return 0;
        }

        return (1 + photoerror)*medium_->get_molDensity()*medium_->get_atomInMolecule().at(i)*particle_->c*particle_->c*integral_->integrateWithLog(min, max, this);
    }
    }
}



//----------------------------------------------------------------------------//

void Photonuclear::setMeasured()
{
    if(initM)
    {
        initM=false;

        double x_aux[]  =   {0, 0.1, 0.144544, 0.20893, 0.301995,
                            0.436516, 0.630957, 0.912011, 1.31826, 1.90546,
                            2.75423, 3.98107, 5.7544, 8.31764, 12.0226,
                            17.378, 25.1189, 36.3078, 52.4807, 75.8577,
                            109.648, 158.489, 229.087, 331.131, 478.63,
                            691.831, 1000, 1445.44, 2089.3, 3019.95,
                            4365.16, 6309.58, 9120.12, 13182.6, 19054.6,
                            27542.3, 39810.8, 57544, 83176.4, 120226,
                            173780, 251188, 363078, 524807, 758576,
                            1.09648e+06, 1.58489e+06, 2.29086e+06, 3.3113e+06, 4.78628e+06,
                            6.91828e+06, 9.99996e+06};

        double y_aux[]  =   {0, 0.0666667, 0.0963626, 159.74, 508.103,
                            215.77, 236.403, 201.919, 151.381, 145.407,
                            132.096, 128.546, 125.046, 121.863, 119.16,
                            117.022, 115.496, 114.607, 114.368, 114.786,
                            115.864, 117.606, 120.011, 123.08, 126.815,
                            131.214, 136.278, 142.007, 148.401, 155.46,
                            163.185, 171.574, 180.628, 190.348, 200.732,
                            211.782, 223.497, 235.876, 248.921, 262.631,
                            277.006, 292.046, 307.751, 324.121, 341.157,
                            358.857, 377.222, 396.253, 415.948, 436.309,
                            457.334, 479.025};

        vector<double> x(x_aux, x_aux + sizeof(x_aux) / sizeof(double) );
        vector<double> y(y_aux, y_aux + sizeof(y_aux) / sizeof(double) );

        interpolateM_   =   new Interpolate(x, y, 4, false, false);
    }
}

//----------------------------------------------------------------------------//

double Photonuclear::measuredSgN(double e)
{
    setMeasured();
    return interpolateM_->interpolateArray(e);
}

//----------------------------------------------------------------------------//
void Photonuclear::enableHardBB()
{
    if(initH)
    {
        initH           =   false;
        double x_aux[]  =   {3, 4, 5, 6, 7, 8, 9};


        double y_aux[][56][7]=
        {
            {
                {7.174409e-4, 1.7132e-3, 4.082304e-3, 8.628455e-3, 0.01244159, 0.02204591, 0.03228755},
                {-0.2436045, -0.5756682, -1.553973, -3.251305, -5.976818, -9.495636, -13.92918},
                {-0.2942209, -0.68615, -2.004218, -3.999623, -6.855045, -10.05705, -14.37232},
                {-0.1658391, -0.3825223, -1.207777, -2.33175, -3.88775, -5.636636, -8.418409},
                {-0.05227727, -0.1196482, -0.4033373, -0.7614046, -1.270677, -1.883845, -2.948277},
                {-9.328318e-3, -0.02124577, -0.07555636, -0.1402496, -0.2370768, -0.3614146, -0.5819409},
                {-8.751909e-4, -1.987841e-3, -7.399682e-3, -0.01354059, -0.02325118, -0.03629659, -0.059275},
                {-3.343145e-5, -7.584046e-5, -2.943396e-4, -5.3155e-4, -9.265136e-4, -1.473118e-3, -2.419946e-3}
            },
            {
                {-1.269205e-4, -2.843877e-4, -5.761546e-4, -1.195445e-3, -1.317386e-3, -9.689228e-15, -6.4595e-15},
                {-0.01563032, -0.03589573, -0.07768545, -0.157375, -0.2720009, -0.4186136, -0.8045046},
                {0.04693954, 0.1162945, 0.3064255, 0.7041273, 1.440518, 2.533355, 3.217832},
                {0.05338546, 0.130975, 0.3410341, 0.7529364, 1.425927, 2.284968, 2.5487},
                {0.02240132, 0.05496, 0.144945, 0.3119032, 0.5576727, 0.8360727, 0.8085682},
                {4.658909e-3, 0.01146659, 0.03090286, 0.06514455, 0.1109868, 0.1589677, 0.1344223},
                {4.822364e-4, 1.193018e-3, 3.302773e-3, 6.843364e-3, 0.011191, 0.015614, 0.01173827},
                {1.9837e-5, 4.940182e-5, 1.409573e-4, 2.877909e-4, 4.544877e-4, 6.280818e-4, 4.281932e-4}
            }
        };

        interpolateH_   =   (Interpolate *)calloc(hmax,sizeof(Interpolate));

        for(int i=0; i<hmax; i++)
        {
            vector<double> x(x_aux, x_aux + sizeof(x_aux) / sizeof(double) );
            vector<double> y(y_aux[particle_->type-1][i], y_aux[particle_->type-1][i] + sizeof(y_aux[particle_->type-1][i]) / sizeof(double) );

            interpolateH_[i]    =   Interpolate(x, y, 4, false, false) ;
        }
    }
}

//----------------------------------------------------------------------------//

double Photonuclear::hardBB(double e, double v)
{
    enableHardBB();

    if(e<1.e5 || v<1.e-7)
    {
        return 0;
    }

    double aux, sum, lov, loe;

    sum =   0;
    aux =   1;
    lov =   log(v)/LOG10;
    loe =   log(e)/LOG10-3;

    for(int i=0; i<hmax; i++)
    {
        if(i>0)
        {
            aux *=  lov;
        }

        sum +=  aux*interpolateH_[i].interpolateArray(loe);
    }

    return sum/v;
}

//----------------------------------------------------------------------------------------------------//

double Photonuclear::function(double Q2)
{

    double x, aux, nu, G, F2, R2;

    nu  =   v*particle_->e;
    x   =   Q2/(2*medium_->get_M().at(i)*nu);

    if(medium_->get_NucCharge().at(i)==1)
    {
        G   =   1;
    }
    else switch(shadow)
        {
        case 1:
        {
            if(x<0.0014)
            {
                G   =   pow(medium_->get_atomicNum().at(i), -0.1);
            }
            else if(x<0.04)
            {
                G   =   pow(medium_->get_atomicNum().at(i), 0.069*log(x)/LOG10+0.097);
            }
            else
            {
                G   =   1;
            }

            break;

        }
        case 2:
        default:
        {
            if(x>0.3)
            {
                const double Mb =   0.437;
                const double la =   0.5;
                const double x2 =   0.278;

                double mb, Aosc, mu, au, ac;

                mb      =   Mb*medium_->get_mN().at(i);
                au      =   1/(1 - x);
                ac      =   1/(1 - x2);
                mu      =   MPI/medium_->get_M().at(i);
                Aosc    =   (1 - la*x)*((au - ac)-mu*(au*au - ac*ac));
                G       =   1 - mb*Aosc;
            }
            else
            {
                const double M1 =   0.129;
                const double M2 =   0.456;
                const double M3 =   0.553;

                double m1, m2, m3, x0, sgn;

                m1  =   M1*medium_->get_mN().at(i);
                m2  =   M2*medium_->get_mN().at(i);
                m3  =   M3*medium_->get_mN().at(i);
                nu  *=  1.e-3;
                sgn =   112.2*(0.609*pow(nu, 0.0988) + 1.037*pow(nu, -0.5944));
                aux =   0.00282*pow(medium_->get_atomicNum().at(i), 1./3)*sgn;
                G   =   (3/aux)*(0.5 + ((1 + aux)*exp(-aux)-1)/(aux*aux));
                G   =   0.75*G + 0.25;
                x0  =   pow(G/(1+m2), 1/m1);

                if(x>=x0)
                {
                    G   =   pow(x, m1)*(1+m2)*(1-m3*x);
                }
            }
        }
        }

    switch(form)
    {
    case 3:  // Abramowicz, Levin, Levy, and Maor parametrization
    {
        double P, W2;

        aux =   x*x;
        P   =   1 - 1.85*x + 2.45*aux - 2.35*aux*x + aux*aux;
        G   *=  (medium_->get_NucCharge().at(i) + (medium_->get_atomicNum().at(i) - medium_->get_NucCharge().at(i))*P);
        W2  =   medium_->get_M().at(i)*medium_->get_M().at(i) - Q2 + 2*medium_->get_M().at(i)*particle_->e*v;

        double cp1, cp2, cp3, cr1, cr2, cr3, ap1, ap2, ap3, ar1, ar2, ar3;
        double bp1, bp2, bp3, br1, br2, br3, m2o, m2r, L2, m2p, Q2o;

        switch(bb)
        {
        case 1:  // ALLM91
            cp1     =   0.26550;
            cp2     =   0.04856;
            cp3     =   1.04682;
            cr1     =   0.67639;
            cr2     =   0.49027;
            cr3     =   2.66275;
            ap1     =   -0.04503;
            ap2     =   -0.36407;
            ap3     =   8.17091;
            ar1     =   0.60408;
            ar2     =   0.17353;
            ar3     =   1.61812;
            bp1     =   0.49222;
            bp2     =   0.52116;
            bp3     =   3.55115;
            br1     =   1.26066;
            br2     =   1.83624;
            br3     =   0.81141;
            m2o     =   0.30508;
            m2r     =   0.20623;
            L2      =   0.06527;
            m2p     =   10.67564;
            Q2o     =   0.27799;
            break;
        case 2:  // ALLM97
        default:
            cp1     =   0.28067;
            cp2     =   0.22291;
            cp3     =   2.1979;
            cr1     =   0.80107;
            cr2     =   0.97307;
            cr3     =   3.4942;
            ap1     =   -0.0808;
            ap2     =   -0.44812;
            ap3     =   1.1709;
            ar1     =   0.58400;
            ar2     =   0.37888;
            ar3     =   2.6063;
            bp1     =   0.60243;
            bp2     =   1.3754;
            bp3     =   1.8439;
            br1     =   0.10711;
            br2     =   1.9386;
            br3     =   0.49338;
            m2o     =   0.31985;
            m2r     =   0.15052;
            L2      =   0.06527;
            m2p     =   49.457;
            Q2o     =   0.46017;
        }

        // GeV -> MeV conversion
        m2o     *=  1e6;
        m2r     *=  1e6;
        L2      *=  1e6;
        m2p     *=  1e6;
        Q2o     *=  1e6;

        // these values are corrected according to the file f2allm.f from Halina Abramowicz
        bp1     *=  bp1;
        bp2     *=  bp2;
        br1     *=  br1;
        br2     *=  br2;
        Q2o     +=  L2;

        const double R  =   0;

        double cr, ar, cp, ap, br, bp, t;

        t   =   log(log((Q2 + Q2o)/L2)/log(Q2o/L2));

        if(t<0)
        {
            t=0;
        }

        cr  =   cr1 + cr2*pow(t, cr3);
        ar  =   ar1 + ar2*pow(t, ar3);
        cp  =   cp1 + (cp1 - cp2)*(1/(1 + pow(t, cp3)) - 1);
        ap  =   ap1 + (ap1 - ap2)*(1/(1 + pow(t, ap3)) - 1);
        br  =   br1 + br2*pow(t, br3);
        bp  =   bp1 + bp2*pow(t, bp3);

        double xp, xr, F2p, F2r;

        xp  =   (Q2 + m2p)/(Q2 + m2p + W2 - medium_->get_M().at(i)*medium_->get_M().at(i));
        xr  =   (Q2 + m2r)/(Q2 + m2r + W2 - medium_->get_M().at(i)*medium_->get_M().at(i));
        F2p =   cp*pow(xp, ap)*pow(1 - x, bp);
        F2r =   cr*pow(xr, ar)*pow(1 - x, br);
        F2  =   (Q2/(Q2 + m2o))*(F2p + F2r)*G;
        R2  =   (2*(1 + R));

        break;
    }
    case 4:  // Butkevich and Mikheyev parametrization
    default:
    {
        const double a  =   0.2513e6;
        const double b  =   0.6186e6;
        const double c  =   3.0292e6;
        const double d  =   1.4817e6;
        const double d0 =   0.0988;
        const double ar =   0.4056;
        const double t  =   1.8152;
        const double As =   0.12;
        const double Bu =   1.2437;
        const double Bd =   0.1853;
        const double R  =   0.25;

        double F2p, F2n, FSp, FNp, FSn, FNn, n, dl, xUv, xDv;

        n   =   1.5*(1 + Q2/(Q2 + c));
        dl  =   d0*(1 + 2*Q2/(Q2 + d));
        aux =   As*pow(x, -dl)*pow(Q2/(Q2 + a), 1 + dl);
        FSp =   aux*pow(1 - x, n + 4);
        FSn =   aux*pow(1 - x, n + t);
        aux =   pow(x, 1 - ar)*pow(1 - x, n)*pow(Q2/(Q2 + b), ar);
        xUv =   Bu*aux;
        xDv =   Bd*aux*(1 - x);
        FNp =   xUv + xDv;
        FNn =   xUv/4 + xDv*4;
        F2p =   FSp + FNp;
        F2n =   FSn + FNn;
        F2  =   G*(medium_->get_NucCharge().at(i)*F2p + (medium_->get_atomicNum().at(i) - medium_->get_NucCharge().at(i))*F2n);
        R2  =   (2*(1 + R));
    }
    }

    aux =   ME*RE/Q2;
    aux *=  aux*(1 - v - medium_->get_M().at(i)*x*v/(2*particle_->e) + (1 - 2*particle_->m*particle_->m/Q2)*v*v*(1 + 4*medium_->get_M().at(i)*medium_->get_M().at(i)*x*x/Q2)/R2);

    return (4*PI*F2/v)*aux;
}

//----------------------------------------------------------------------------//

double Photonuclear::functionInt(double e, double v)
{
    particle_->setEnergy(e);
    setEnergy(cros->get_component());

    if(vUp==vMax)
    {
        return 0;
    }

    v   =   vUp*exp(v*log(vMax/vUp));

    return photoN(v, cros->get_component());
}

//----------------------------------------------------------------------------//
// Getter

PhotoContinuous  *Photonuclear::get_Continuous()
{
    return continuous_;
}

PhotoStochastic *Photonuclear::get_Stochastic()
{
    return stochastic_;
}
//----------------------------------------------------------------------------//
//Setter

void Photonuclear::set_jt(bool newjt)
{
    jt_ =   newjt;
}

void Photonuclear::set_form(int newForm)
{
    form    =   newForm;
}

void Photonuclear::set_bb(int newbb)
{
    bb  =   newbb;
}

void Photonuclear::set_shadow(int newshadow)
{
    shadow  =   newshadow;
}

void Photonuclear::SetRandomNumberGenerator(boost::function<double ()> &f)
{
	MathModel::SetRandomNumberGenerator(f);
	get_Continuous()->SetRandomNumberGenerator(f);
	get_Stochastic()->SetRandomNumberGenerator(f);
}
