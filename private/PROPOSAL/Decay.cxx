/*
* Decay.cxx
*
*  Created on: 24.06.2010
*      Author: koehne
*/

#include "PROPOSAL/Decay.h"
#include "algorithm"
#include <cmath>
#include "PROPOSAL/methods.h"
#include "PROPOSAL/PROPOSALParticle.h"

bool Decay::flag=false;

using namespace std;

Decay::Decay(CrossSections *cros) : CrossSections(cros)
{
    f       =   new FindRoot(IMAXS, IPREC);
}

//----------------------------------------------------------------------------//

double Decay::decay()
{

    if(cros->get_cd()<=0 || particle_->l<0)
    {
        return 0;
    }

    return cros->get_cd()/max((particle_->p/particle_->m)*particle_->l*SPEED, XRES);
}

//----------------------------------------------------------------------------//

double Decay::e(double ernd, double arnd, double srnd, Output *o)
{
    if(particle_->l<0)
    {
        return 0;
    }
    double emax, x0, f0, el, lm, pl;
    string out1 =   "nu";
    string out2 =   "nu";

    if(particle_->type==2)
    {

        const double brmu   =   0.1737;
        const double brel   =   0.1783+brmu;
        const double brpi   =   0.1109+brel;
        const double br2p   =   0.2540+brpi;
        const double br3p   =   0.1826+br2p;

        if(srnd<brmu)
        {
            lm  =   MMU;

            if(endsWith(particle_->name,"+"))
            {
                out =   "mu+";
                if(flag)
                {
                    out1    =   "~nu_tau";
                    out2    =   "nu_mu";
                }
            }
            else if(endsWith(particle_->name,"-"))
            {
                out =   "mu-";
                if(flag)
                {
                    out1    =   "nu_tau";
                    out2    =   "~nu_mu";
                }
            }
            else
            {
                out =    "mu";
                if(flag)
                {
                    out1    =   "nu_tau";
                    out2    =   "~nu_mu";
                }
            }
        }
        else if(srnd<brel)
        {
            lm  =   ME;
            if(endsWith(particle_->name,"+"))
            {
                if(o->get_AMASIM())
                {
                    out =    "epair";
                }
                else
                {
                    out =   "e+";
                }

                if(flag)
                {
                    out1    =   "~nu_tau";
                    out2    =   "nu_e";
                }
            }
            else if(endsWith(particle_->name,"-"))
            {
                if(o->get_AMASIM())
                {
                    out =   "delta";
                }
                else
                {
                    out =   "e-";
                }

                if(flag)
                {
                    out1    =   "nu_tau";
                    out2    =   "~nu_e";
                }
            }
            else
            {
                if(o->get_AMASIM())
                {
                    out =   "delta";
                }
                else
                {
                    out =   "e";
                }

                if(flag)
                {
                    out1    =   "nu_tau";
                    out2    =   "~nu_e";
                }
            }
        }
        else
        {
            if(srnd<brpi)
            {
                lm  =   MPI;
            }
            else if(srnd<br2p)
            {
                lm  =   MRH;
            }
            else if(srnd<br3p)
            {
                lm  =   MA1;
            }
            else
            {
                lm  =   MRS;
            }
            el  =   (MTAU*MTAU + lm*lm) / (2*MTAU);

            if(o->get_AMASIM())
            {
                out =   "munu";
            }
            else
            {
                out =   "hadr";
            }

            el  =   el * (particle_->e/particle_->m) + sqrt(el*el - lm*lm) * (particle_->p/particle_->m) * (2*arnd-1);

            if(flag)
            {
                if(endsWith(particle_->name,"+"))
                {
                    out1    =   "~nu_tau";
                }
                else if(endsWith(particle_->name,"-"))
                {
                    out1    =   "nu_tau";
                }
                else
                {
                    out1    =   "nu_tau";
                }
                o->output(1, out1, particle_->e-el, 0);
            }
            return el;
        }
    }
    else
    {
        lm  =   ME;
        if(endsWith(particle_->name,"+"))
        {
            if(o->get_AMASIM())
            {
                out =   "epair";
            }
            else
            {
                out =   "e+";
            }
            if(flag)
            {
                out1    =   "~nu_mu";
                out2    =   "nu_e";
            }
        }
        else if(endsWith(particle_->name,"-"))
        {
            if(o->get_AMASIM())
            {
                out =   "delta";
            }
            else
            {
                out =   "e-";
            }
            if(flag)
            {
                out1    =   "nu_mu";
                out2    =   "~nu_e";
            }
        }
        else
        {
            if(o->get_AMASIM())
            {
                out =   "delta";
            }
            else
            {
                out =   "e";
            }
            if(flag)
            {
                out1    =   "nu_mu";
                out2    =   "~nu_e";
            }
        }
    }
    emax    =   (particle_->m*particle_->m + lm*lm) / (2*particle_->m);
    x0      =   lm/emax;
    f0      =   x0*x0;
    f0      =   f0*x0 - f0*f0/2;
    el      =   max(f->findRoot(x0, 1.0, 0.5, this, f0+(0.5-f0)*ernd)*emax , lm);;
    pl      =   sqrt(el*el - lm*lm);

    if(flag)
    {
        int sign;
        double cp, sp, ct, st, en, En1, En2;
        cp  =   pl/(particle_->m - el);
        sp  =   sqrt(1 - cp*cp);

        if(arnd<0.5)
        {
            sign    =   1;
            arnd    =   2*arnd;
        }
        else
        {
            sign    =   -1;
            arnd    =   2*arnd-1;
        }
        ct  =   2*arnd - 1;
        st  =   sqrt(1 - ct*ct);
        en  =   (particle_->m - el)/2;
        En1 =   -cp*ct - sign*sp*st;
        En2 =   -cp*ct + sign*sp*st;
        En1 =   en * ((particle_->e/particle_->m) + (particle_->p/particle_->m) * En1);
        En2 =   en * ((particle_->e/particle_->m) + (particle_->p/particle_->m) * En2);
        o->output(1, out1, En1, 0);
        o->output(1, out2, En2, 0);
    }

    return el * (particle_->e/particle_->m) + pl * (particle_->p/particle_->m) * (2*arnd - 1);
}

//----------------------------------------------------------------------------//

double Decay::function(double x)
{
    double x2;
    x2  =   x*x;
    return  x*x2 - x2*x2/2;
}

//----------------------------------------------------------------------------//

double Decay::dFunction(double x)
{
    return (3 - 2*x) * x*x;
}


