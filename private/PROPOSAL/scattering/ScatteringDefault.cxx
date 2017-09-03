
#include <boost/bind.hpp>

#include "PROPOSAL/math/MathModel.h"
#include "PROPOSAL/scattering/ScatteringDefault.h"
#include "PROPOSAL/crossection/CrossSections.h"
#include "PROPOSAL/crossection/Ionization.h"
#include "PROPOSAL/crossection/Bremsstrahlung.h"
#include "PROPOSAL/crossection/Epairproduction.h"
#include "PROPOSAL/crossection/Photonuclear.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;
using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ScatteringDefault::ScatteringDefault()
    :do_interpolation_(false)
    ,order_of_interpolation_(5)
    ,integral_(IROMB, IMAXS, IPREC2)
    ,interpolant_(NULL)
    ,interpolant_diff_(NULL)
{
}

ScatteringDefault::ScatteringDefault(const ScatteringDefault &scattering)
    :do_interpolation_(scattering.do_interpolation_)
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

    integral_ = Integral(scattering.integral_);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// ScatteringDefault& ScatteringDefault::operator=(const ScatteringDefault &scattering){
//     if (this != &scattering)
//     {
//       ScatteringDefault tmp(scattering);
//       swap(tmp);
//     }
//     return *this;
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// bool ScatteringDefault::operator==(const ScatteringDefault &scattering) const
// {
//     if( interpolant_ != NULL && scattering.interpolant_ != NULL)
//     {
//         if( *interpolant_   != *scattering.interpolant_) return false;
//     }
//
//     if( interpolant_diff_ != NULL && scattering.interpolant_diff_ != NULL)
//     {
//         if( *interpolant_diff_   != *scattering.interpolant_diff_) return false;
//     }
//
//     if( integral_   != scattering.integral_) return false;
//
//     //else
//     return true;
//
// }
//
//
// //----------------------------------------------------------------------------//
// //----------------------------------------------------------------------------//
//
//
// bool ScatteringDefault::operator!=(const ScatteringDefault &scattering) const {
//   return !(*this == scattering);
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void ScatteringDefault::swap(ScatteringDefault &scattering)
// {
//     using std::swap;
//
//     swap(do_interpolation_,scattering.do_interpolation_);
//
//     if(scattering.interpolant_ != NULL)
//     {
//         interpolant_->swap(*scattering.interpolant_);
//     }
//     else
//     {
//         interpolant_ = NULL;
//     }
//
//     if(scattering.interpolant_diff_ != NULL)
//     {
//         interpolant_diff_->swap(*scattering.interpolant_diff_) ;
//     }
//     else
//     {
//         interpolant_diff_ = NULL;
//     }
//
//     integral_.swap(scattering.integral_);
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//


double ScatteringDefault::FunctionToIntegral(const PROPOSALParticle& particle, const std::vector<CrossSections*>& cross_sections, double energy)
{
    double aux, aux2;
    double result;

    //Same Implentation as in double ProcessCollection::FunctionToIntegral(double energy)
    //It was reimplemented to avoid building a ProcessCollecion Object to calculate or
    //test the scattering class.
    PROPOSALParticle tmp_particle(particle);
    tmp_particle.SetEnergy(energy);
    result  =    0;

    for(unsigned int i =0;i<cross_sections.size();i++)
    {
        aux     =   cross_sections.at(i)->CalculatedEdx(particle.GetEnergy());
        result  +=  aux;
    }

    aux = -1/result;
    //End of the reimplementation

    aux2    =   RY*particle.GetEnergy() / (particle.GetMomentum() *particle.GetMomentum());
    aux     *=  aux2*aux2;

    return aux;
}

//----------------------------------------------------------------------------//

long double ScatteringDefault::CalculateTheta0(const PROPOSALParticle& particle, const std::vector<CrossSections*>& cross_sections, double dr, double ei, double ef)
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
        aux = integral_.Integrate(ei, ef, boost::bind(&ScatteringDefault::FunctionToIntegral, this, boost::cref(particle), boost::cref(cross_sections), _1),4);
    }
    // TODO: Decide if using RadiationLength calculated in Bremsstrahlung or Medium
    double radiation_lenght = cross_sections.at(0)->GetMedium()->GetRadiationLength();

    // TODO: check if one has to take the absolute value of the particle charge
    aux =   sqrt(max(aux, 0.0) / radiation_lenght) * particle.GetCharge();
    aux *=  max(1 + 0.038*log(dr / radiation_lenght), 0.0);

    return min(aux, cutoff);
}

//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringDefault::CalculateRandomAngle(const PROPOSALParticle& particle, const std::vector<CrossSections*>& cross_sections, double dr, double ei, double ef)
{
    double Theta0,rnd1,rnd2;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(particle, cross_sections, dr, ei, ef);

    rnd1 = SQRT2*Theta0*erfInv( 2.*(RandomGenerator::Get().RandomDouble()-0.5) );
    rnd2 = SQRT2*Theta0*erfInv( 2.*(RandomGenerator::Get().RandomDouble()-0.5) );

    random_angles.sx = (rnd1/SQRT3+rnd2)/2;
    random_angles.tx = rnd2;

    rnd1 = SQRT2*Theta0*erfInv(2*(RandomGenerator::Get().RandomDouble()-0.5));
    rnd2 = SQRT2*Theta0*erfInv(2*(RandomGenerator::Get().RandomDouble()-0.5));

    random_angles.sy = (rnd1/SQRT3+rnd2)/2;
    random_angles.ty = rnd2;

    return random_angles;
}

//----------------------------------------------------------------------------//

double ScatteringDefault::FunctionToBuildInterpolant(const PROPOSALParticle& particle, const std::vector<CrossSections*>& cross_sections, double energy)
{
        return integral_.Integrate(energy, BIGENERGY, boost::bind(&ScatteringDefault::FunctionToIntegral, this, boost::cref(particle), boost::cref(cross_sections), _1),4);
}

void ScatteringDefault::EnableInterpolation(const PROPOSALParticle& particle, const std::vector<CrossSections*>& cross_sections, std::string path)
{
    if(do_interpolation_)return;

    bool reading_worked=true, storing_failed=false;

    // charged anti leptons have the same cross sections like charged leptons
    // (except of diffractive Bremsstrahlung, where one can analyse the interference term if implemented)
    // so they use the same interpolation tables
    string particle_name = particle.GetName();

    if(!path.empty())
    {
        stringstream filename;
        filename<<path<<"/ScatteringDefault"
                <<"_"<<particle_name
                <<"_mass_"<<particle.GetMass()
                <<"_charge_"<<particle.GetCharge()
                <<"_lifetime_"<<particle.GetLifetime()
                <<"_"<< cross_sections.at(0)->GetMedium()->GetName()
                <<"_"<< cross_sections.at(0)->GetMedium()->GetMassDensity()
                <<"_"<< cross_sections.at(0)->GetEnergyCutSettings()->GetEcut()
                <<"_"<< cross_sections.at(0)->GetEnergyCutSettings()->GetVcut();

        for(unsigned int i =0; i<cross_sections.size(); i++)
        {
            switch (cross_sections.at(i)->GetType())
            {
                case ParticleType::Brems:
                    filename << "_b"
                        << "_" << cross_sections.at(i)->GetParametrization()
                        << "_" << cross_sections.at(i)->GetLpmEffectEnabled();
                    break;
                case ParticleType::DeltaE:
                    filename << "_i";
                    break;
                case ParticleType::EPair:
                    filename << "_e"
                        << "_" << cross_sections.at(i)->GetLpmEffectEnabled();
                    break;
                case ParticleType::NuclInt:
                    filename << "_p"
                        << "_" << cross_sections.at(i)->GetParametrization();
                    break;
                default:
                    log_fatal("Unknown cross section");
                    exit(1);
            }
            filename<< "_" << cross_sections.at(i)->GetMultiplier()
                    << "_" << cross_sections.at(i)->GetEnergyCutSettings()->GetEcut()
                    << "_" << cross_sections.at(i)->GetEnergyCutSettings()->GetVcut();

        }


        if( FileExist(filename.str()) )
        {
            log_debug( "ScatteringDefault tables will be read from file:\t%s",filename.str().c_str());
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
                log_info( "File %s is corrupted! Write it again!",filename.str().c_str());
            }

            log_info("ScatteringDefault tables will be saved to file:\t%s",filename.str().c_str());

            ofstream output;

            output.open(filename.str().c_str());

            if(output.good())
            {
                output.precision(16);

                interpolant_ = new Interpolant(NUM2,particle.GetLow() , BIGENERGY ,boost::bind(&ScatteringDefault::FunctionToBuildInterpolant, this, boost::cref(particle), boost::cref(cross_sections), _1), order_of_interpolation_ , false, false, true, order_of_interpolation_, false, false, false );
                interpolant_diff_ = new Interpolant(NUM2,particle.GetLow() , BIGENERGY ,boost::bind(&ScatteringDefault::FunctionToIntegral, this, boost::cref(particle), boost::cref(cross_sections), _1), order_of_interpolation_ , false, false, true, order_of_interpolation_, false, false, false );

                interpolant_->Save(output);
                interpolant_diff_->Save(output);
            }
            else
            {
                storing_failed  =   true;
                log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
            }

            output.close();
        }
    }
    if(path.empty() || storing_failed)
    {
        interpolant_ = new Interpolant(NUM2,particle.GetLow() , BIGENERGY ,boost::bind(&ScatteringDefault::FunctionToBuildInterpolant, this, boost::cref(particle), boost::cref(cross_sections), _1), order_of_interpolation_ , false, false, true, order_of_interpolation_, false, false, false );
        interpolant_diff_ = new Interpolant(NUM2,particle.GetLow() , BIGENERGY ,boost::bind(&ScatteringDefault::FunctionToIntegral, this, boost::cref(particle), boost::cref(cross_sections), _1), order_of_interpolation_ , false, false, true, order_of_interpolation_, false, false, false );
    }


    do_interpolation_ = true;
}

void ScatteringDefault::DisableInterpolation()
{
    delete interpolant_;
    delete interpolant_diff_;

    interpolant_ = NULL;
    interpolant_diff_ = NULL;

    do_interpolation_ = false;
}
