/*! \file   Propagate.cxx
*   \brief  Source file for the propagate routines.
*
*   For more details see the class documentation.
*
*   \date   05.07.2010
*   \author Jan-Hendrik Koehne
*/

#include "PROPOSAL/Propagate.h"
#include <math.h>
#include "PROPOSAL/methods.h"
#include "algorithm"
#include "stdlib.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/Energy2Loss.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/StandardNormal.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/StandardNormal.h"
#include "PROPOSAL/IonizContinuous.h"
#include "PROPOSAL/IonizStochastic.h"
#include "PROPOSAL/EpairStochastic.h"
#include "PROPOSAL/EpairContinuous.h"
#include "PROPOSAL/BremsContinuous.h"
#include "PROPOSAL/BremsStochastic.h"
#include "PROPOSAL/Decay.h"
#include "PROPOSAL/Energy2LossX.h"
#include "PROPOSAL/Energy2LossE.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/PhotoStochastic.h"
#include "PROPOSAL/PhotoContinuous.h"
#include "PROPOSAL/Interpolate.h"
#include <sstream>
#include <fstream>

int Propagate::g    =   5;

using namespace std;

//----------------------------------------------------------------------------//

Propagate::Propagate(){}

//----------------------------------------------------------------------------//

Propagate::Propagate(string w, double ecut, double vcut)
{

    init(w, ecut, vcut, "mu", 1.);
}

//----------------------------------------------------------------------------//

Propagate::Propagate(string w, double ecut, double vcut, string type)
{

    init(w, ecut, vcut, type, 1.);
}

//----------------------------------------------------------------------------//


Propagate::Propagate(string w, double ecut, double vcut, string type, double rho)
{

    init(w, ecut, vcut, type, rho);
}

//----------------------------------------------------------------------------//



void Propagate::init(std::string w, double ecut, double vcut, std::string type, double rho)
{
    this->rho   =   1.;
    sdec        =   false;
    recc        =   false;
    exactTime   =   false;
    contiCorr   =   false;
    molieScat   =   false;
    dw          =   false;
    rw          =   0;
    hw          =   0;
    df          =   false;
    jt          =   false;
    up          =   true;

    particle_   =   new PROPOSALParticle(this, type);
    medium_     =   new Medium(w, ecut, vcut, rho);
    cros        =   new CrossSections(particle_, medium_);
    E2Loss_     =   new Energy2Loss(cros);

    integral_.resize(2);

    for(int i=0; i<2; i++)
    {
        integral_.at(i) =   new Integral(IROMB, IMAXS, IPREC2);
    }

    StandardN   =   new StandardNormal(IROMB, IMAXS, IPREC);
    o           =   new Output(particle_);

}

//----------------------------------------------------------------------------//


double Propagate::propagateTo(double r, double e)
{
    int wint;
    bool DEBUG, flag;
    double ei, ef=0, efd, efi, aux=0, dr;
    double rndd, rndi, rnd1, rnd2, rnd3, rnddMin, rndiMin, rndTot;

    DEBUG=o->DEBUG;

    if(o->HIST==-1)
    {
        o->init(particle_->name);
    }

    ei  =   e;
    ef  =   ei;

    if(r<0)
    {
        r   =   0;
    }


    if(e<=particle_->low || r==0)
    {
        flag    =   false;
    }
    else
    {
        flag    =   true;
    }

    if(DEBUG)
    {

        cerr<<"\nPropagating "<<particle_->name<<" of energy "<<o->f(ei)<<" MeV to a distance of "<<o->f(r)<<" cm"<<endl;
    }

    while(flag)
    {

        rndd    = -log(RandomDouble());
        rndi    = -log(RandomDouble());

        if(DEBUG)
        {
            cerr<<"1. solving the tracking integral ...  rndd = "<<o->f(rndd)<<"  rndi = "<<o->f(rndi)<<" ...  "<<endl;
        }

        if(particle_->l<0)
        {
            rnddMin =   0;
        }
        else
        {

            rnddMin =   getpr(ei, rndd, false)/rho;
        }

        if(DEBUG)
        {
            cerr<<" \t \t \t rnddMin = "<<o->f(rnddMin)<<" (d)  "<<endl;
        }

        rndiMin =   getpr(ei, rndi, true);

        if(DEBUG)
        {
            cerr<<"rndiMin = "<<o->f(rndiMin)<<" (i)"<<endl;
        }

        if(DEBUG)
        {
            cerr<<"2. evaluating the energy loss ...  "<<endl;
        }

        if(rndd>=rnddMin || rnddMin<=0)
        {
            efd =   particle_->low;
        }
        else
        {
            efd =   getef(ei, rndd*rho, false);
        }

        if(DEBUG)
        {
            cerr<<"efd = "<<o->f(efd)<<" MeV  "<<endl;
        }

        if(rndi>=rndiMin || rndiMin<=0)
        {
            efi =   particle_->low;
        }
        else
        {
            efi =   getef(ei, rndi, true);
        }

        if(DEBUG)
        {
            cerr<<"efi = "<<o->f(efi)<<" MeV ...  "<<endl;
        }

        pint    =   (efi>efd);

        if(pint)
        {
            ef  =   efi;
        }
        else
        {
            ef  =   efd;
        }

        if(DEBUG)
        {
            cerr<<" \t \t \t lost "<<o->f(ei-ef)<<" MeV  ef = "+o->f(ef)+" MeV"<<endl;
        }

        if(DEBUG)
        {
            cerr<<"3. calculating the displacement ...  "<<endl;
        }

        dr  =   cros->getdx(ei, ef, rho*(r-particle_->r))/rho;

        if(DEBUG)
        {
            cerr<<"dr = "<<o->f(dr)<<" cm"<<endl;
        }

        if(dr<r-particle_->r)
        {
            if(DEBUG)
            {
                cerr<<"4. calculating the local time ...  "<<endl;
            }
        }
        else
        {
            dr  =   r - particle_->r;

            if(DEBUG)
            {
                cerr<<"4. getting the const energy ...  "<<endl;
            }

            ef  =   cros->getef(ei, rho*dr);

            if(DEBUG)
            {
                cerr<<"lost "<<o->f(ei-ef)<<" MeV  ef = "<<o->f(ef)+" MeV"<<endl;
            }

            if(DEBUG)
            {
                cerr<<"5. calculating the local time ...  "<<endl;
            }
        }

        if(recc)
        {
            o->output(0, "a"+particle_->name, ei, dr);
        }

        particle_->advance(dr, ei, ef);

        if(DEBUG)
        {
            cerr<<"t = "<<o->f(particle_->t)<<" s"<<endl;
        }

        if(abs(r-particle_->r)<abs(r)*COMPUTER_PRECISION)
        {
            particle_->r=r;  // computer precision control
        }

        if(contiCorr)
        {
            if(ef!= particle_->low)
            {
                ef  =   StandardN->sndrn(RandomDouble(), ef, sqrt(E2Loss_->e2le->dE2de(ei, ef)), particle_->low, ei, false);
            }

        }

        if(recc)
        {
            o->output(0, "conti", ei-ef, -dr);
        }

        if(ef==particle_->low || particle_->r==r)
        {
            break;
        }

        if(DEBUG)
        {
            cerr<<"5. choosing the cross section ..."<<endl;
        }

        rnd2    =   RandomDouble();
        rnd3    =   RandomDouble();

        particle_->setEnergy(ef);

        if(pint)
        {
            rnd1    =   RandomDouble();

            if(dw)
            {
                if(particle_->r>hw)
                {
                    double exp  =   abs(rw);
                    double pow_ =   pow(rnd2, exp);

                    if(rw>0)
                    {
                        rnd2    =   1 - pow_*rnd2;
                    }
                    else
                    {
                        rnd2    =   pow_*rnd2;
                    }

                    rw  =   (1 + exp)*pow_;
                    hw  =   particle_->r;
                    dw  =   false;
                }
            }

            if(DEBUG)
            {
                decayS  =   cros->get_decay()->decay();
            }

            ionizS  =   cros->get_ionization()->get_Stochastic()->dNdx(rnd2);
            bremsS  =   cros->get_bremsstrahlung()->get_Stochastic()->dNdx(rnd2);
            photoS  =   cros->get_photonuclear()->get_Stochastic()->dNdx(rnd2);
            epairS  =   cros->get_epairproduction()->stochastic_->dNdx(rnd2);
            totalS  =   ionizS+bremsS+photoS+epairS;
            rndTot  =   rnd1*totalS;

            if(DEBUG)
            {
                cerr<<" . rnd1 = "<<o->f(rnd1)<<" rnd2 = "<<o->f(rnd2)<<
                    " rnd3 = "<<o->f(rnd3)<<" decay = "<<o->f(decayS)<<endl;
            }

            if(DEBUG)
            {
                cerr<<" . ioniz = "<<o->f(ionizS)<<" brems = "<<o->f(bremsS)<<
                    " photo = "<<o->f(photoS)<<" epair = "<<o->f(epairS);
            }

            if(ionizS>rndTot)
            {

                aux     =   cros->get_ionization()->get_Stochastic()->e(rnd3);
                ef      -=  aux;
                wint    =   2;
            }
            else if(ionizS+bremsS>rndTot)
            {
                aux     =   cros->get_bremsstrahlung()->get_Stochastic()->e(rnd3);
                ef      -=  aux;
                wint    =   3;
            }
            else if(ionizS+bremsS+photoS>rndTot)
            {
                aux     =   cros->get_photonuclear()->get_Stochastic()->e(rnd3);
                ef      -=  aux;
                wint    =   4;
            }
            else if(ionizS+bremsS+photoS+epairS>rndTot)
            {
                aux     =   cros->get_epairproduction()->stochastic_->e(rnd3);
                ef      -=  aux;
                wint    =   5;
            }
            else  // due to the parameterization of the cross section cutoffs
            {
                ei  =   ef;
                continue;
            }

        }

        else
        {
            if(particle_->type==2)
            {
                aux     =   cros->get_decay()->e(rnd2, rnd3, RandomDouble(), o);
                ef      =   0;
                wint    =   1;
            }
            else
            {
                aux     =   cros->get_decay()->e(rnd2, rnd3, 0.5, o);
                ef      =   0;
                wint    =   1;
            }
        }

        if(wint==1)
        {
            o->output(wint, cros->get_decay()->get_out(), aux, ef);
        }
        else
        {
            o->output(wint, medium_->get_E(cros->get_component()), aux, ef);
        }



        if(ef<=particle_->low)
        {

            break;
        }

        ei  =   ef;

    }

    if(sdec)
    {
        if(particle_->r!=r && ef!=0 && particle_->l>=0)
        {
            particle_->setEnergy(particle_->m);

            particle_->t    +=  -particle_->l*log(RandomDouble());

            if(particle_->type==2)
            {
                aux =   cros->get_decay()->e(RandomDouble(), 0.5, RandomDouble(), o);
            }
            else
            {
                aux =   cros->get_decay()->e(RandomDouble(), 0.5, 0.5, o);
            }

            ef  =   0;

            o->output(1, cros->get_decay()->get_out(), aux, ef);
        }
    }

    particle_->setEnergy(ef);  // to remember const state of the particle
    o->HIST =   -1;            // to make sure user resets particle properties

    if(particle_->r==r)
    {
        if(DEBUG)
        {
            cerr<<"PROPOSALParticle reached the border with energy ef = "<<o->f(ef)<<" MeV";
        }

        return ef;
    }
    else
    {
        if(DEBUG)
        {
            if(particle_->l<0)
            {
                cerr<<"PROPOSALParticle stopped at rf = "<<o->f(particle_->r)<<" cm";
            }
            else
            {
                cerr<<"PROPOSALParticle disappeared at rf = "<<o->f(particle_->r)<<" cm";
            }
        }

        return -particle_->r;
    }
}

//----------------------------------------------------------------------------//


double Propagate::propagateTo(double r, double e, double rc)
{
    int HIST;
    double result;

    HIST    =   o->HIST;
    ec      =   propagateTo(rc, e);
    tc      =   particle_->t;

    if(ec>0)
    {
        o->HIST =   HIST;
        result  =   propagateTo(r, ec);
    }
    else
    {
        result  =   ec;
    }

    return result;
}

//----------------------------------------------------------------------------//


double Propagate::getPropEc()
{
    return ec;
}

//----------------------------------------------------------------------------//



double Propagate::getPropTc()
{
    return tc;
}

//----------------------------------------------------------------------------//


double Propagate::function(double E)
{
    const bool DEBUG    =   false;
    double aux;

    aux =   cros->function(E);

    if(!pint)
    {
        decayS  =   cros->get_decay()->decay();

        if(DEBUG)
        {
            o->err<<" + "+o->f(particle_->e);
        }

        return aux*decayS;
    }
    else
    {
        ionizS  =   cros->get_ionization()->get_Stochastic()->dNdx();

        if(DEBUG)
        {
            o->err<<" \t "+o->f(ionizS);
        }

        bremsS  =   cros->get_bremsstrahlung()->get_Stochastic()->dNdx();

        if(DEBUG)
        {
            o->err<<" \t "+o->f(bremsS);
        }

        photoS  =   cros->get_photonuclear()->get_Stochastic()->dNdx();

        if(DEBUG)
        {
            o->err<<" \t "+o->f(photoS);
        }

        epairS  =   cros->get_epairproduction()->stochastic_->dNdx();

        if(DEBUG)
        {
            o->err<<" \t "+o->f(epairS);
        }

        totalS  =   ionizS + bremsS + photoS + epairS;

        return aux*totalS;
    }
}

//----------------------------------------------------------------------------//


void Propagate::interpolate(string w)
{


    int i;
    stringstream tmpStr;

    double e_hi     =   PhysicsModel::get_ebig();
    double e_low    =   particle_->low;

    cout<<"Parameterizations apply in the energy range from "<<e_low<<" MeV to "<<e_hi<<" MeV"<<endl;

    if(w.find("all")!=string::npos)
    {
        w   +=  " crs trackE trackX";

        if(exactTime)
        {
            w   +=  " trackT";
        }

        if(contiCorr)
        {
            w   +=  " contiR";
        }

        if(molieScat)
        {
            w   +=  " molieS gaussR";
        }
    }

    if(w.find("crs")!=string::npos)
    {
        w   +=  " ionizE ionizS bremsE bremsS photoE photoS epairE epairP epairS";

        if(cros->get_photonuclear()->get_form()>2)
        {
            w   +=  " photoP";
        }
    }

    if(w.find("contiR")!=string::npos)
    {
        w   +=  " gaussR en2ldX en2ldE";
    }

    w=toLowerCase(w);

    try
    {
        cros->get_ionization()->get_Continuous()->set_jt(false);

        if(w.find("ionize")!=string::npos && cros->get_ci()>0)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "ionize" << "-" << NUM1;
            cout<<"Parameterizing ionizE ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->get_ionization()->get_Continuous()->interpolateJ_ = new Interpolate(NUM1, e_low, e_hi, cros->get_ionization()->get_Continuous(), g, true, false, true, g, false, false, true);
            cros->get_ionization()->get_Continuous()->set_jt(true);
            cout<<"done \n";
        }

        cros->get_ionization()->get_Stochastic()->set_jt(false);

        if(w.find("ionizs")!=string::npos && cros->get_ci()>0)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "ionizs" << "-" << NUM1;
            cout<<"Parameterizing ionizS ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->get_ionization()->get_Stochastic()->interpolateJ_     =   new Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, cros->get_ionization()->get_Stochastic(), g, false, false, true, g, false, false, false, g, true, false, false);
            cros->get_ionization()->get_Stochastic()->interpolateJo_    =   new Interpolate(NUM1, e_low, e_hi, cros->get_ionization()->get_Stochastic(), g, false, false, true, g, true, false, false);
            cros->get_ionization()->get_Stochastic()->set_jt(true);
            cout<<"done "<<endl;

        }

        cros->get_bremsstrahlung()->get_Continuous()->set_jt(false);

        if(w.find("bremse")!=string::npos && cros->get_cb()>0)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "bremse" << "-" << NUM1;
            cout<<"Parameterizing bremsE ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->get_bremsstrahlung()->get_Continuous()->interpolateJ_ =   new Interpolate(NUM1, e_low, e_hi, cros->get_bremsstrahlung()->get_Continuous(), g, true, false, true, g, false, false, false);
            cros->get_bremsstrahlung()->get_Continuous()->set_jt(true);
            cout<<"done \n";
        }

        cros->get_bremsstrahlung()->get_Stochastic()->set_jt(false);

        if(w.find("bremss")!=string::npos && cros->get_cb()>0)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "bremss" << "-" << NUM1;
            cout<<"Parameterizing bremsS ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->get_bremsstrahlung()->get_Stochastic()->interpolateJ_     =   (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));
            cros->get_bremsstrahlung()->get_Stochastic()->interpolateJo_    =   (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));

            for(i=0; i<medium_->get_numCompontents(); i++)
            {
                cros->set_component(i);
                cros->get_bremsstrahlung()->get_Stochastic()->interpolateJ_[i]  =   Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, cros->get_bremsstrahlung()->get_Stochastic(), g, false, false, true, g, false, false, false, g, true, false, false);
                cros->get_bremsstrahlung()->get_Stochastic()->interpolateJo_[i] =   Interpolate(NUM1, e_low, e_hi, cros->get_bremsstrahlung()->get_Stochastic(), g, false, false, true, g, true, false, false);

            }

            cros->get_bremsstrahlung()->get_Stochastic()->set_jt(true);
            cout<<"done \n";
        }

        cros->get_photonuclear()->set_jt(false);

        if(w.find("photop")!=string::npos && cros->get_cp()>0)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "photop" << "-" << NUM1;
            cout<<"Parameterizing photoP ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->get_photonuclear()->interpolateJ_ =   (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));

            for(i=0; i<medium_->get_numCompontents(); i++)
            {
                cros->set_component(i);
                cros->get_photonuclear()->interpolateJ_[i]  =   Interpolate(NUM1, e_low, e_hi, NUM1, 0., 1., cros->get_photonuclear(), g, false, false, true, g, false, false, false, g, false, false, false);

            }

            cros->get_photonuclear()->set_jt(true);
            cout<<"done \n";
        }

        cros->get_photonuclear()->get_Continuous()->set_jt(false);

        if(w.find("photoe")!=string::npos && cros->get_cp()>0)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "photoe" << "-" << NUM1;
            cout<<"Parameterizing photoE ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->get_photonuclear()->get_Continuous()->interpolateJ_   =   new Interpolate(NUM1, e_low, e_hi, cros->get_photonuclear()->get_Continuous(), g, true, false, true, g, false, false, false);
            cros->get_photonuclear()->get_Continuous()->set_jt(true);
            cout<<"done \n";
        }

        cros->get_photonuclear()->get_Stochastic()->set_jt(false);

        if(w.find("photos")!=string::npos && cros->get_cp()>0)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "photos" << "-" << NUM1;
            cout<<"Parameterizing photoS ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->get_photonuclear()->get_Stochastic()->interpolateJ_   =   (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));
            cros->get_photonuclear()->get_Stochastic()->interpolateJo_  =   (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));

            for(i=0; i<medium_->get_numCompontents(); i++)
            {
                cros->set_component(i);
                cros->get_photonuclear()->get_Stochastic()->interpolateJ_[i]    =   Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, cros->get_photonuclear()->get_Stochastic(), g, false, false, true, g, false, false, false, g, true, false, false);
                cros->get_photonuclear()->get_Stochastic()->interpolateJo_[i]   =   Interpolate(NUM1, e_low, e_hi, cros->get_photonuclear()->get_Stochastic(), g, false, false, true, g, true, false, false);
            }

            cros->get_photonuclear()->get_Stochastic()->set_jt(true);
            cout<<"done \n";
        }

        cros->get_epairproduction()->set_jt(false);

        if(w.find("epairp")!=string::npos && cros->get_ce()>0)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "epairp" << "-" << NUM1;
            cout<<"Parameterizing epairP ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->get_epairproduction()->interpolateJ_  =   (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));

            for(i=0; i<medium_->get_numCompontents(); i++)
            {
                cros->set_component(i);
                cros->get_epairproduction()->interpolateJ_[i]   =   Interpolate(NUM1, e_low, e_hi, NUM1, 0., 1., cros->get_epairproduction(), g, false, false, true, g, false, false, false, g, false, false, false);
            }

            cros->get_epairproduction()->set_jt(true);
            cout<<"done \n";
        }

        cros->get_epairproduction()->get_Continuous()->set_jt(false);

        if(w.find("epaire")!=string::npos && cros->get_ce()>0)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "epaire" << "-" << NUM1;
            cout<<"Parameterizing epairE ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->get_epairproduction()->get_Continuous()->interpolateJ_    =   new Interpolate(NUM1, e_low, e_hi, cros->get_epairproduction()->continuous_, g, true, false, true, g, false, false, false);
            cros->get_epairproduction()->get_Continuous()->set_jt(true);
            cout<<"done \n";

        }

        cros->get_epairproduction()->get_Stochastic()->set_jt(false);

        if(w.find("epairs")!=string::npos && cros->get_ce()>0)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "epairs" << "-" << NUM1;
            cout<<"Parameterizing epairS ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->get_epairproduction()->get_Stochastic()->interpolateJ_    =   (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));
            cros->get_epairproduction()->get_Stochastic()->interpolateJo_   =   (Interpolate *)calloc(medium_->get_numCompontents(),sizeof(Interpolate));

            for(i=0; i<medium_->get_numCompontents(); i++)
            {
                cros->set_component(i);
                cros->get_epairproduction()->get_Stochastic()->interpolateJ_[i]     =   Interpolate(NUM1, e_low, e_hi, NUM1, 0, 1, cros->get_epairproduction()->stochastic_, g, false, false, true, g, false, false, false, g, true, false, false);
                cros->get_epairproduction()->get_Stochastic()->interpolateJo_[i]    =   Interpolate(NUM1, e_low, e_hi, cros->get_epairproduction()->stochastic_, g, false, false, true, g, true, false, false);

            }

            cros->get_epairproduction()->get_Stochastic()->set_jt(true);
            cout<<"done \n";
        }

        this->jt=false;

        if(w.find("tracke")!=string::npos)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "tracke" << "-" << NUM3;
            cout<<"Parameterizing trackE ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            this->interpolateJ_     =   (Interpolate *)calloc(2,sizeof(Interpolate));
            this->interpolateJdf_   =   (Interpolate *)calloc(2,sizeof(Interpolate));
            pint=true;

            if(abs(-integral_.at(1)->integrateWithLog(particle_->low, particle_->low*10, this))<abs(-integral_.at(1)->integrateWithLog(PhysicsModel::get_ebig(), PhysicsModel::get_ebig()/10, this)))
            {
                up  =   true;
            }
            else
            {
                up  =   false;
            }

            for(i=0; i<2; i++)
            {
                pint                        =   (i==1);
                this->df                    =   false;
                this->interpolateJ_[i]      =   Interpolate(NUM3, e_low, e_hi, this, g, false, false, true, g, false, false, false);
                this->df                    =   true;
                this->interpolateJdf_[i]    =   Interpolate(NUM3, e_low, e_hi, this, g, false, false, true, g, false, false, false);

            }

            this->jt    =   true;

            for(i=0; i<2; i++)
            {
                if(!(up&&(i==1)))
                {
                    bigLow[i]=interpolateJ_[i].interpolate(particle_->low);
                }
            }

            if(up)
            {
                cout<<"done \n";
            }
            else
            {
                cout<<"down \n";
            }
        }

        cros->set_jt(false);

        if(w.find("trackx")!=string::npos)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "trackx" << "-" << NUM3;
            cout<<"Parameterizing trackX ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cros->set_df(false);
            cros->interpolateJ_     =   new Interpolate(NUM3, e_low, e_hi, cros, g, false, false, true, g, false, false, false);
            cros->set_df(true);
            cros->interpolateJdf_   =   new Interpolate(NUM3, e_low, e_hi, cros, g, false, false, true, g, false, false, false);
            cros->set_jt(true);
            cout<<"done \n";

        }

        particle_->jt   =   false;

        if(w.find("trackt")!=string::npos)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "trackt" << "-" << NUM3;
            cout<<"Parameterizing trackT ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            particle_->df               =   false;
            particle_->interpolateJ_    =   new Interpolate(NUM3, e_low, e_hi, particle_, g, false, false, true, g, false, false, false);
            particle_->df               =   true;
            particle_->interpolateJdf_  =   new Interpolate(NUM3, e_low, e_hi, particle_, g, false, false, true, g, false, false, false);
            particle_->jt               =   true;
            cout<<"done \n";
        }

        particle_->get_scattering()->set_jt(false);

        if(w.find("molies")!=string::npos)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "molies" << "-" << NUM2;
            cout<<"Parameterizing molieS ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            particle_->get_scattering()->set_df(false);
            particle_->get_scattering()->interpolateJ_      =   new Interpolate(NUM2, e_low, e_hi, particle_->get_scattering(), g, false, false, true, g, false, false, false);
            particle_->get_scattering()->set_df(true);
            particle_->get_scattering()->interpolateJdf_    =   new Interpolate(NUM2, e_low, e_hi, particle_->get_scattering(), g, false, false, true, g, false, false, false);
            particle_->get_scattering()->set_jt(true);
            cout<<"done"<<endl;
            ;
        }

        ;

        StandardN->set_jt(false);

        if(w.find("gaussr")!=string::npos)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "gaussr" << "-" << NUM2;
            cout<<"Parameterizing gaussR ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            StandardN->interpolateJ_    =   new Interpolate(NUM2, -5, 5, StandardN, g, true, false, false, g, true, false, false);
            StandardN->set_jt(true);
            cout<<"done \n";
            ;
        }

        E2Loss_->e2lx->set_jt(false);

        if(w.find("en2ldx")!=string::npos)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "en2ldx" << "-" << NUM2;
            cout<<"Parameterizing en2ldx ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            cout<<"Parameterizing en2ldX ... ";
            E2Loss_->e2lx->interpolateJ_    =   new Interpolate(NUM2, e_low, e_hi, E2Loss_->e2lx, g, false, false, true, g, false, false, false);
            E2Loss_->e2lx->set_jt(true);
            cout<<"done \n";
        }

        E2Loss_->e2le->set_jt(false);

        if(w.find("en2lde")!=string::npos)
        {
            tmpStr.clear();
            tmpStr.str("");
            tmpStr << "en2lde" << "-" << NUM2;
            cout<<"Parameterizing en2ldE ... ";

            if(Output::inf&&!Output::raw)
            {
                if((tmpStr.str().compare(Output::readStr()) == 0))
                {
                    cout << "reading table" << " ... ";
                }
                else
                {
                    throw 0;
                }
            }

            if(Output::outf)
            {
                Output::write(tmpStr.str());
                cout << "writing table" << " ... ";
            }

            E2Loss_->e2le->set_df(false);
            E2Loss_->e2le->interpolateJ_    =   new Interpolate(NUM2, e_low, e_hi, E2Loss_->e2le, g, false, false, true, g, false, false, false);
            E2Loss_->e2le->set_df(true);
            E2Loss_->e2le->interpolateJdf_  =   new Interpolate(NUM2, e_low, e_hi, E2Loss_->e2le, g, false, false, true, g, false, false, false);
            E2Loss_->e2le->set_jt(true);
            cout<<"done \n";
        }

        cout<<"Finished parameterizations"<<endl;
    }
    catch(int a)
    {
        throw 0;
    }

}

//----------------------------------------------------------------------------//

void Propagate::interpolate(string w, string filename)
{
    stringstream name;
    bool flag;

    name<<filename<<".";
    name<<particle_->name<<"_"<<replaceAll((toLowerCase((medium_->get_name()+"_"+w))),' ', '-');
    name<<"_"<<Output::f(medium_->get_ecut())<<"_"<<Output::f(medium_->get_vcut())<<"_";
    name<<cros->get_bremsstrahlung()->get_form();
    name<<cros->get_photonuclear()->get_form();
    name<<cros->get_photonuclear()->get_bb();
    name<<cros->get_photonuclear()->get_shadow();

    if(cros->get_lpm())
    {
        name<<"t";
    }
    else
    {
        name<<"f";
    }

    if(exactTime)
    {
        name<<"t";
    }
    else
    {
        name<<"f";
    }

    if(contiCorr)
    {
        name<<"t";
    }
    else
    {
        name<<"f";
    }

    if(molieScat)
    {
        name<<"t";
    }
    else
    {
        name<<"f";
    }

    if(PhysicsModel::get_elow()==particle_->low)
    {
        name<<"_l"<<Output::f(PhysicsModel::get_elow());
    }

    if(PhysicsModel::get_ebig()!=BIGENERGY)
    {
        name<<"_b"<<Output::f(PhysicsModel::get_ebig());
    }

    if(cros->get_ci()!=1 || cros->get_cb()!=1 || cros->get_cp()!=1 || cros->get_ce()!=1 || cros->get_cd()!=1 || medium_->get_rho()!=1)
    {
        name<<"_"<<Output::f(cros->get_ci());
        name<<","<<Output::f(cros->get_cb());
        name<<","<<Output::f(cros->get_cp());
        name<<","<<Output::f(cros->get_ce());
        name<<","<<Output::f(cros->get_cd());
        name<<","<<Output::f(medium_->get_rho());
    }

    if(Output::raw)
    {
        name<<"_raw";
    }
    else
    {
        name<<"_ascii";
    }

    name<<".data";

    do
    {

        if(Output::texi)
        {
            return;
        }

        flag    =   false;

        try
        {
            Output::open(name.str());
            interpolate(w);
            Output::close();
        }
        catch (int a)
        {
            cout<<"EXCEPTION"<<endl;
            flag    =   true;
            Output::Delete(name.str());
        }
    }
    while(flag);
}

//----------------------------------------------------------------------------//



double Propagate::getpr(double ei, double rnd, bool pint)
{
    if(jt)
    {
        if(pint)
        {
            storeDif[1] =   interpolateJ_[1].interpolate(ei);
        }
        else
        {
            storeDif[0] =   interpolateJ_[0].interpolate(ei);
        }


        if(up&&pint)
        {
            if(pint)
            {
                return max(storeDif[1], 0.0);
            }
            else
            {
                return max(storeDif[0], 0.0);
            }
        }
        else
        {
            if(pint)
            {
                return max(bigLow[1]-storeDif[1], 0.0);
            }
            else
            {
                return max(bigLow[0]-storeDif[0], 0.0);
            }
        }
    }
    else
    {

        this->pint  =   pint;

        if(pint)
        {
            return integral_.at(1)->integrateWithLog(ei, particle_->low, this, -rnd);
        }
        else
        {
            return integral_.at(0)->integrateWithLog(ei, particle_->low, this, -rnd);
        }
    }
}

//----------------------------------------------------------------------------------------------------//


double Propagate::getef(double ei, double rnd, bool pint)
{
    if(jt)
    {
        if(pint)
        {
            if(abs(rnd)>abs(storeDif[1])*HALF_PRECISION)
            {
                double aux;

                if(up&&pint)
                {
                    if(pint)
                    {
                        aux =   interpolateJ_[1].findLimit(storeDif[1]-rnd);
                    }
                    else
                    {
                        aux =   interpolateJ_[0].findLimit(storeDif[0]-rnd);
                    }
                }
                else
                {
                    if(pint)
                    {
                        aux =   interpolateJ_[1].findLimit(storeDif[1]+rnd);
                    }
                    else
                    {
                        aux =   interpolateJ_[0].findLimit(storeDif[0]+rnd);
                    }
                }

                if(abs(ei-aux)>abs(ei)*HALF_PRECISION)
                {
                    return min(max(aux, particle_->low), ei);
                }
            }
        }
        else
        {
            if(abs(rnd)>abs(storeDif[0])*HALF_PRECISION)
            {
                double aux;

                if(up&&pint)
                {
                    if(pint)
                    {
                        aux =   interpolateJ_[1].findLimit(storeDif[1]-rnd);
                    }
                    else
                    {
                        aux =   interpolateJ_[0].findLimit(storeDif[0]-rnd);
                    }
                }
                else
                {
                    if(pint)
                    {
                        aux =   interpolateJ_[1].findLimit(storeDif[1]+rnd);
                    }
                    else
                    {
                        aux =   interpolateJ_[0].findLimit(storeDif[0]+rnd);
                    }
                }

                if(abs(ei-aux)>abs(ei)*HALF_PRECISION)
                {
                    return min(max(aux, particle_->low), ei);
                }
            }
        }

        if(pint)
        {
            return min(max(ei + rnd/interpolateJdf_[1].interpolate(ei + rnd/(2*interpolateJdf_[1].interpolate(ei))), particle_->low), ei);
        }
        else
        {
            return min(max(ei + rnd/interpolateJdf_[0].interpolate(ei + rnd/(2*interpolateJdf_[0].interpolate(ei))), particle_->low), ei);
        }
    }
    else
    {
        this->pint  =   pint;

        if(pint)
        {
            return integral_.at(1)->getUpperLimit();
        }
        else
        {
            return integral_.at(0)->getUpperLimit();
        }

    }
}

//----------------------------------------------------------------------------//


double Propagate::functionInt(double e)
{
    if(df)
    {
        return function(e);
    }
    else
    {
        if(up&&pint)
        {
            if(pint)
            {
                return integral_.at(1)->integrateWithLog(e, particle_->low, this);
            }
            else
            {
                return integral_.at(0)->integrateWithLog(e, particle_->low, this);
            }
        }
        else
        {
            if(pint)
            {
                return -integral_.at(1)->integrateWithLog(e, PhysicsModel::get_ebig(), this);
            }
            else
            {
                return -integral_.at(0)->integrateWithLog(e, PhysicsModel::get_ebig(), this);
            }
        }
    }
}

//----------------------------------------------------------------------------//

void Propagate::setCSform(int bspar, int pncrs, int pncbb, int pncsh, string muta)
{

    if(bspar<1 || bspar>4)
    {
        cout<<"Warning: bs is not a valid number"<<endl;
        bspar   =   1;
    }

    if(pncrs<1 || pncrs>4 || (pncrs==2 && !(muta.compare("mu")==0 || muta.compare("tau")==0)))
    {
        cout<<"Warning: ph is not a valid number"<<endl;
        pncrs   =   1;
    }

    if(((pncrs==1 || pncrs==2) && (pncbb<1 || pncbb>4)) ||
            (pncrs==3 && (pncbb<1 || pncbb>2)) || (pncrs==4 && pncbb!=1))
    {
        cout<<"Warning: bb is not a valid number"<<endl;
        pncbb   =   1;
    }

    if(((pncrs==1 || pncrs==2) && (pncsh!=1)) || ((pncrs>2) && (pncsh<1 || pncsh>2)))
    {
        cout<<"Warning: sh is not a valid number"<<endl;
        pncsh   =   1;
    }

    cout<<"bspar="<<bspar<<"\t";
    cout<<"pncrs="<<pncrs<<"\t";
    cout<<"pncbb="<<pncbb<<"\t";
    cout<<"pncsh="<<pncsh<<endl;

    cros->get_bremsstrahlung()->set_form(bspar);    // Choose parametrization of the bremsstrahlung cross section
    cros->get_photonuclear()->set_form(pncrs);      // Choose parametrization of the photon-nucleon cross section
    cros->get_photonuclear()->set_bb(pncbb);        // Choose parametrization of the photon-nucleon cross section
    cros->get_photonuclear()->set_shadow(pncsh);    // Choose parametrization of the photon-nucleon cross section

}

//----------------------------------------------------------------------------//

void Propagate::setCShigh()
{
    int bspar=2;
    int pncrs=1;
    int pncbb=2;
    int pncsh=2;

    cros->get_bremsstrahlung()->set_form(bspar);    // Choose parametrization of the bremsstrahlung cross section
    cros->get_photonuclear()->set_form(pncrs);      // Choose parametrization of the photon-nucleon cross section
    cros->get_photonuclear()->set_bb(pncbb);        // Choose parametrization of the photon-nucleon cross section
    cros->get_photonuclear()->set_shadow(pncsh);    // Choose parametrization of the photon-nucleon cross section

}

//----------------------------------------------------------------------------//

void Propagate::setCSlow()
{
    int bspar=3;
    int pncrs=3;
    int pncbb=2;
    int pncsh=2;

    cros->get_bremsstrahlung()->set_form(bspar);    // Choose parametrization of the bremsstrahlung cross section
    cros->get_photonuclear()->set_form(pncrs);      // Choose parametrization of the photon-nucleon cross section
    cros->get_photonuclear()->set_bb(pncbb);        // Choose parametrization of the photon-nucleon cross section
    cros->get_photonuclear()->set_shadow(pncsh);    // Choose parametrization of the photon-nucleon cross section
}


//----------------------------------------------------------------------------//
// Getter

StandardNormal *Propagate::get_Standard()
{
    return StandardN;
}

CrossSections *Propagate::get_cros()
{
    return cros;
}

Medium *Propagate::get_Medium()
{
    return medium_;
}
//----------------------------------------------------------------------------//
//Setter

void Propagate::set_exactTime(bool eTime)
{
    exactTime   =   eTime;
}

void Propagate::set_molieScat(bool mScat)
{
    molieScat   =   mScat;
}

void Propagate::set_contiCorr(bool cCorr)
{
    contiCorr   =   cCorr;
}

void Propagate::set_pint(bool newPint)
{
    pint        =   newPint;
}

void Propagate::set_rho(double newRho)
{
    rho         =   newRho;
}
