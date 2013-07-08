/*! \file   Scattering.cxx
*   \brief  Source file for the scattering routines.
*
*   For more details see the class documentation.
*
*   \date   2013.06.13
*   \author Tomasz Fuchs
*/


#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Scattering.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Photonuclear.h"

using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Scattering::Scattering( )
    :x0_(0)
    ,do_interpolation_(false)
    ,order_of_interpolation_(5)
    ,integral_(  new Integral(IROMB, IMAXS, IPREC2) )
    ,interpolant_(NULL)
    ,interpolant_diff_(NULL)
    ,particle_( new Particle("mu") )
{
    crosssections_.push_back(new Ionization());
    crosssections_.push_back(new Bremsstrahlung());
    crosssections_.push_back(new Photonuclear());
    crosssections_.push_back(new Epairproduction());

    SetParticle(particle_);

    for(unsigned int i =0;i<crosssections_.size();i++)
    {
        if( crosssections_.at(i)->GetName().compare("Bremsstrahlung") ==0){
            Bremsstrahlung* tmp = (Bremsstrahlung*)crosssections_.at(i);
            x0_ = tmp->CalculateScatteringX0();
        }
    }
    standard_normal_    =   new StandardNormal(IROMB, IMAXS, IPREC);
}

//double Scattering::cutoff   =   1;


Scattering::Scattering(std::vector<CrossSections*> crosssections)
    :x0_(0)
    ,do_interpolation_(false)
    ,order_of_interpolation_(5)
    ,integral_( new Integral(IROMB, IMAXS, IPREC2) )
    ,interpolant_(NULL)
    ,interpolant_diff_(NULL)
    ,particle_(NULL)
    ,crosssections_( crosssections)
{
    SetParticle(crosssections.at(0)->GetParticle());
    for(unsigned int i =0;i<crosssections_.size();i++)
    {
        if( crosssections_.at(i)->GetName().compare("Bremsstrahlung") ==0){
            Bremsstrahlung* tmp = (Bremsstrahlung*)crosssections_.at(i);
            x0_ = tmp->CalculateScatteringX0();
        }
    }

    standard_normal_    =   new StandardNormal(IROMB, IMAXS, IPREC);
}

Scattering::Scattering(const Scattering &scattering)
    :x0_(scattering.x0_)
    ,do_interpolation_(scattering.do_interpolation_)
    ,particle_(scattering.particle_)
    ,crosssections_(scattering.crosssections_)
{
    if(scattering.interpolant_ != NULL)
    {
        interpolant_ = new Interpolant(*scattering.interpolant_) ;
    }
    else
    {
        interpolant_ = NULL;
    }

    if(scattering.interpolant_diff_ != NULL)
    {
        interpolant_diff_ = new Interpolant(*scattering.interpolant_diff_) ;
    }
    else
    {
        interpolant_diff_ = NULL;
    }

    integral_ = new Integral(* scattering.integral_);
    standard_normal_ = new StandardNormal( *scattering.standard_normal_);
};

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Scattering& Scattering::operator=(const Scattering &scattering){
    if (this != &scattering)
    {
      Scattering tmp(scattering);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Scattering::operator==(const Scattering &scattering) const
{
    if(*particle_ != *(scattering.particle_))return false;

    for(unsigned int i =2; i<scattering.crosssections_.size(); i++)
    {
        if(scattering.crosssections_.at(i)->GetName().compare("Bremsstrahlung")==0)
        {
            if( *(Bremsstrahlung*)crosssections_.at(i) !=  *(Bremsstrahlung*)scattering.crosssections_.at(i) ) return false;
        }
        else if(scattering.crosssections_.at(i)->GetName().compare("Ionization")==0)
        {
            if( *(Ionization*)crosssections_.at(i) != *(Ionization*)scattering.crosssections_.at(i) ) return false;
        }
        else if(scattering.crosssections_.at(i)->GetName().compare("Epairproduction")==0)
        {
            if( *(Epairproduction*)crosssections_.at(i) !=  *(Epairproduction*)scattering.crosssections_.at(i) ) return false;
        }
        else if(scattering.crosssections_.at(i)->GetName().compare("Photonuclear")==0)
        {
            if( *(Photonuclear*)crosssections_.at(i) !=  *(Photonuclear*)scattering.crosssections_.at(i) )  return false;
        }
        else
        {
            cout<<"In copy constructor of Scattering: Error: Unknown crossSection"<<endl;
            exit(1);
        }
    }

    if( standard_normal_ != NULL && scattering.standard_normal_ != NULL)
    {
        if( *standard_normal_   != *scattering.standard_normal_)                                        return false;
    }

    if( interpolant_ != NULL && scattering.interpolant_ != NULL)
    {
        if( *interpolant_   != *scattering.interpolant_)                                        return false;
    }

    if( interpolant_diff_ != NULL && scattering.interpolant_diff_ != NULL)
    {
        if( *interpolant_diff_   != *scattering.interpolant_diff_)                                        return false;
    }

    if( integral_ != NULL && scattering.integral_ != NULL)
    {
        if( *integral_   != *scattering.integral_)                                        return false;
    }

    //else
    return true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Scattering::operator!=(const Scattering &scattering) const {
  return !(*this == scattering);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Scattering::swap(Scattering &scattering)
{
    using std::swap;

    Particle tmp_particle1(*scattering.particle_);
    Particle tmp_particle2(*particle_);

    vector<CrossSections*> tmp_cross1(scattering.crosssections_);
    vector<CrossSections*> tmp_cross2(crosssections_);

    SetCrosssections(  tmp_cross1 );
    scattering.SetCrosssections(  tmp_cross2 );

    SetParticle( new Particle(tmp_particle1) );
    scattering.SetParticle( new Particle(tmp_particle2) );

    swap(x0_,scattering.x0_);
    swap(do_interpolation_,scattering.do_interpolation_);

    if(scattering.interpolant_ != NULL)
    {
        interpolant_->swap(*scattering.interpolant_);
    }
    else
    {
        interpolant_ = NULL;
    }

    if(scattering.interpolant_diff_ != NULL)
    {
        interpolant_diff_->swap(*scattering.interpolant_diff_) ;
    }
    else
    {
        interpolant_diff_ = NULL;
    }

    integral_->swap(*scattering.integral_);

    if(scattering.standard_normal_ != NULL)
    {
        standard_normal_->swap(*scattering.standard_normal_) ;
    }
    else
    {
        standard_normal_ = NULL;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


/*! \file   Scattering.cxx
*   \brief  Source file for the scattering routines.
*
*   For more details see the class documentation.
*
*   \date   29.06.2010
*   \author Jan-Hendrik Koehne
*/



//----------------------------------------------------------------------------//


double Scattering::FunctionToIntegral(double energy)
{
    double aux, aux2;
    double result;

    //Same Implentation as in double ProcessCollection::FunctionToIntegral(double energy)
    //It was reimplemented to avoid building a ProcessCollecion Object to calculate or
    //test the scattering class.
    particle_->SetEnergy(energy);
    result  =    0;

    for(unsigned int i =0;i<crosssections_.size();i++)
    {
        aux     =   crosssections_.at(i)->CalculatedEdx();
        result  +=  aux;
    }

    aux = -1/result;
    //End of the reimplementation

    aux2    =   RY*particle_->GetEnergy() / (particle_->GetMomentum() *particle_->GetMomentum());
    aux     *=  aux2*aux2;

    return aux;
}

//----------------------------------------------------------------------------//

long double Scattering::CalculateTheta0(double dr, double ei, double ef)
{

    double aux=-1;
    double cutoff=1;
    if(do_interpolation_)
    {
        if(fabs(ei-ef)>fabs(ei)*HALF_PRECISION)
        {
            aux         =   interpolant_->Interpolate(ei);
            double aux2 =   aux - interpolant_->Interpolate(ef);

            if(fabs(aux2)>fabs(aux)*HALF_PRECISION)
            {
                aux =   aux2;
            }
            else
            {
                aux =   interpolant_diff_->Interpolate((ei+ef)/2)*(ef-ei);
            }
        }
        else
        {
            aux =   interpolant_diff_->Interpolate((ei+ef)/2)*(ef-ei);
        }
    }
    else
    {
        aux = integral_->Integrate(ei, ef, boost::bind(&Scattering::FunctionToIntegral, this, _1),4);
    }

    aux =   sqrt(max(aux, 0.0)/x0_) *particle_->GetCharge();

    //TOMASZ IMPLENETATION
    //This should be the correct result, but is way out of range
    /*
    double beta = sqrt(1 -  particle_->GetMass() * particle_->GetMass()/ (particle_->GetEnergy()*particle_->GetEnergy() ));
    beta = (RY/(beta*particle_->GetMomentum())) * particle_->GetCharge() * sqrt(dr/x0_) ;
    cout << "auxAlt: " << aux << endl;
    cout << "auxNeu: " << beta << endl;
    */
    aux *=  max(1 + 0.038*log(dr/x0_), 0.0);

    return min(aux, cutoff);
}

void Scattering::Scatter(double dr, double ei, double ef)
{
    //    Implement the Molie Scattering here see PROPOSALParticle::advance of old version
        long double Theta0, Theta_max,rnd1,rnd2,sx,tx,sy,ty,sz,tz,ax,ay,az;
        double x,y,z;
        Theta0     =   CalculateTheta0(dr, ei, ef);
        Theta_max    =   1./SQRT2;


        rnd1    =   (long double)standard_normal_-> StandardNormalRandomNumber(RandomDouble(), 0, Theta0, -Theta_max, Theta_max, false);
        rnd2    =   (long double)standard_normal_-> StandardNormalRandomNumber(RandomDouble(), 0, Theta0, -Theta_max, Theta_max, false);
        sx      =   (rnd1/SQRT3+rnd2)/2;
        tx      =   rnd2;


        rnd1    =   (long double)standard_normal_-> StandardNormalRandomNumber(RandomDouble(), 0, Theta0, -Theta_max, Theta_max, false);
        rnd2    =   (long double)standard_normal_-> StandardNormalRandomNumber(RandomDouble(), 0, Theta0, -Theta_max, Theta_max, false);
        sy      =   (rnd1/SQRT3+rnd2)/2;
        ty      =   rnd2;


        sz      =   sqrt(max(1.-(sx*sx+sy*sy), (long double)0.));
        tz      =   sqrt(max(1.-(tx*tx+ty*ty), (long double)0.));


        long double sinth, costh,sinph,cosph;
        long double theta,phi;
        sinth = (long double)particle_->GetSinTheta();
        costh = (long double)particle_->GetCosTheta();
        sinph = (long double)particle_->GetSinPhi();
        cosph = (long double)particle_->GetCosPhi();
        x   = particle_->GetX();
        y   = particle_->GetY();
        z   = particle_->GetZ();


        ax      =   sinth*cosph*sz+costh*cosph*sx-sinph*sy;
        ay      =   sinth*sinph*sz+costh*sinph*sx+cosph*sy;
        az      =   costh*sz-sinth*sx;

        x       +=  ax*dr;
        y       +=  ay*dr;
        z       +=  az*dr;


        ax      =   sinth*cosph*tz+costh*cosph*tx-sinph*ty;
        ay      =   sinth*sinph*tz+costh*sinph*tx+cosph*ty;
        az      =   costh*tz-sinth*tx;



        costh   =   az;
        sinth   =   sqrt(max(1-costh*costh, (long double)0));

        if(sinth!=0)
        {
            sinph   =   ay/sinth;
            cosph   =   ax/sinth;
        }

        if(costh>1)
        {
            theta   =   acos(1)*180./PI;
        }
        else if(costh<-1)
        {
            theta   =   acos(-1)*180./PI;
        }
        else
        {
            theta   =   acos(costh)*180./PI;
        }

        if(cosph>1)
        {
            phi =   acos(1)*180./PI;
        }
        else if(cosph<-1)
        {
            phi =   acos(-1)*180./PI;
        }
        else
        {
            phi =   acos(cosph)*180./PI;
        }

        if(sinph<0)
        {
            phi =   360.-phi;
        }

        if(phi>=360)
        {
            phi -=  360.;
        }

        particle_->SetX(x);
        particle_->SetY(y);
        particle_->SetZ(z);
        particle_->SetPhi(phi);
        particle_->SetTheta(theta);
}

//----------------------------------------------------------------------------//

double Scattering::FunctionToBuildInterpolant(double energy)
{
        return integral_->Integrate(energy, BIGENERGY, boost::bind(&Scattering::FunctionToIntegral, this, _1),4);
}

void Scattering::EnableInterpolation(string path)
{
    if(do_interpolation_)return;

    bool reading_worked=true, storing_failed=false;
    if(!path.empty())
    {
        for(unsigned int i = 0; i< crosssections_.size(); i++)
        {
            crosssections_.at(i)->EnableDEdxInterpolation(path);
            crosssections_.at(i)->EnableDNdxInterpolation(path);
        }

        stringstream filename;
        filename<<path<<"/Scattering_"<<particle_->GetName()
               <<"_"<<crosssections_.at(0)->GetMedium()->GetName()
               <<"_"<< crosssections_.at(0)->GetEnergyCutSettings()->GetEcut()
               <<"_"<<crosssections_.at(0)->GetEnergyCutSettings()->GetVcut();

        for(unsigned int i =0; i<crosssections_.size(); i++)
        {
            if(crosssections_.at(i)->GetName().compare("Bremsstrahlung")==0)
            {
                filename << "_b_"
                         << "_" << crosssections_.at(i)->GetParametrization()
                         << "_" << crosssections_.at(i)->GetMultiplier()
                         << "_" << crosssections_.at(i)->GetLpmEffectEnabled()
                         << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetEcut()
                         << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetVcut();

            }
            else if(crosssections_.at(i)->GetName().compare("Ionization")==0)
            {
                filename << "_i_"
                         << "_" << crosssections_.at(i)->GetMultiplier()
                         << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetEcut()
                         << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetVcut();
            }
            else if(crosssections_.at(i)->GetName().compare("Epairproduction")==0)
            {
                filename << "_e_"
                         << "_" << crosssections_.at(i)->GetMultiplier()
                         << "_" << crosssections_.at(i)->GetLpmEffectEnabled()
                         << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetEcut()
                         << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetVcut();
            }
            else if(crosssections_.at(i)->GetName().compare("Photonuclear")==0)
            {
                filename << "_p_"
                         << "_" << crosssections_.at(i)->GetParametrization()
                         << "_" << crosssections_.at(i)->GetMultiplier()
                         << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetEcut()
                         << "_" << crosssections_.at(i)->GetEnergyCutSettings()->GetVcut();
            }

        }


        if( FileExist(filename.str()) )
        {
            cerr<<"Info: Scattering tables will be read from file:"<<endl;
            cerr<<"\t"<<filename.str()<<endl;
            ifstream input;

            input.open(filename.str().c_str());

            interpolant_ = new Interpolant();
            interpolant_diff_ = new Interpolant();
            reading_worked = interpolant_->Load(input);
            reading_worked = reading_worked && interpolant_diff_->Load(input);

            input.close();
        }

        if(!FileExist(filename.str()) || !reading_worked )
        {
            if(!reading_worked)
            {
                cerr<<"Info: file "<<filename.str()<<" is corrupted! Write is again!"<< endl;
            }

            cerr<<"Info: Scattering tables will be saved to file:"<<endl;
            cerr<<"\t"<<filename.str()<<endl;

            double energy = particle_->GetEnergy();

            ofstream output;
            output.open(filename.str().c_str());

            if(output.good())
            {
                output.precision(16);

                interpolant_ = new Interpolant(NUM2,particle_->GetLow() , BIGENERGY ,boost::bind(&Scattering::FunctionToBuildInterpolant, this, _1), order_of_interpolation_ , false, false, true, order_of_interpolation_, false, false, false );
                interpolant_diff_ = new Interpolant(NUM2,particle_->GetLow() , BIGENERGY ,boost::bind(&Scattering::FunctionToIntegral, this, _1), order_of_interpolation_ , false, false, true, order_of_interpolation_, false, false, false );

                interpolant_->Save(output);
                interpolant_diff_->Save(output);
            }
            else
            {
                storing_failed  =   true;
                cerr<<"Warning: Can not open file "<<filename.str()<<" for writing!"<<endl;
                cerr<<"\t Table will not be stored!"<<endl;
            }
            particle_->SetEnergy(energy);

            output.close();
        }
    }
    if(path.empty() || storing_failed)
    {
        double energy = particle_->GetEnergy();

        for(unsigned int i = 0; i< crosssections_.size(); i++)
        {
            crosssections_.at(i)->EnableDEdxInterpolation(path);
            crosssections_.at(i)->EnableDNdxInterpolation(path);
        }

        interpolant_ = new Interpolant(NUM2,particle_->GetLow() , BIGENERGY ,boost::bind(&Scattering::FunctionToBuildInterpolant, this, _1), order_of_interpolation_ , false, false, true, order_of_interpolation_, false, false, false );
        interpolant_diff_ = new Interpolant(NUM2,particle_->GetLow() , BIGENERGY ,boost::bind(&Scattering::FunctionToIntegral, this, _1), order_of_interpolation_ , false, false, true, order_of_interpolation_, false, false, false );

        particle_->SetEnergy(energy);
    }


    do_interpolation_ = true;
}

void Scattering::DisableInterpolation()
{
    delete interpolant_;
    delete interpolant_diff_;

    interpolant_ = NULL;
    interpolant_diff_ = NULL;

    do_interpolation_ = false;
}

//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void Scattering::SetParticle(Particle* particle)
{
    if(particle == NULL || particle == particle_)return;

    for(unsigned int i = 0 ; i< crosssections_.size() ; i++)
    {
        crosssections_.at(i)->SetParticle(particle);
    }

    particle_ = particle;
}

void Scattering::SetCrosssections(std::vector<CrossSections*> crosssections)
{
    crosssections_ = crosssections;
}




