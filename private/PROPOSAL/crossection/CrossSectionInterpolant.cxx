
#include <boost/bind.hpp>
#include <cmath>

#include "PROPOSAL/crossection/CrossSectionInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

CrossSectionInterpolant::CrossSectionInterpolant(Parametrization& param)
    : CrossSection(param)
    , dedx_interpolant_(NULL)
    , dndx_interpolant_1d_(param.GetMedium().GetNumComponents(), NULL)
    , dndx_interpolant_2d_(param.GetMedium().GetNumComponents(), NULL)
{
    Parametrization::Definition param_def = parametrization_.GetDefinition();

    std::vector<Interpolant1DBuilder> builder1d(components_.size());
    std::vector<Interpolant2DBuilder> builder2d(components_.size());

    Helper::InterpolantBuilderContainer builder_container1d(components_.size());
    Helper::InterpolantBuilderContainer builder_container2d(components_.size());
    Helper::InterpolantBuilderContainer builder_return;

    for (unsigned int i = 0; i < components_.size(); ++i)
    {
        // !!! IMPORTANT !!!
        // Order of builder matter because the functions needed for 1d interpolation
        // needs the already intitialized 2d interpolants.
        builder2d[i].SetMax1(NUM1)
            .SetX1Min(param.GetParticleDef().low)
            .SetX1Max(BIGENERGY)
            .SetMax2(NUM1)
            .SetX2Min(0.0)
            .SetX2Max(1.0)
            .SetRomberg1(param_def.order_of_interpolation)
            .SetRational1(false)
            .SetRelative1(false)
            .SetIsLog1(true)
            .SetRomberg2(param_def.order_of_interpolation)
            .SetRational2(false)
            .SetRelative2(false)
            .SetIsLog2(false)
            .SetRombergY(param_def.order_of_interpolation)
            .SetRationalY(true)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction2D(boost::bind(&CrossSection::FunctionToBuildDNdxInterpolant2D, this, _1, _2, i));

        builder_container2d[i].first = &builder2d[i];
        builder_container2d[i].second = &dndx_interpolant_2d_[i];

        builder1d[i].SetMax(NUM1)
            .SetXMin(param.GetParticleDef().low)
            .SetXMax(BIGENERGY)
            .SetRomberg(param_def.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(param_def.order_of_interpolation)
            .SetRationalY(true)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(boost::bind(&CrossSection::FunctionToBuildDNdxInterpolant, this, _1, i));

        builder_container1d[i].first = &builder1d[i];
        builder_container1d[i].second = &dndx_interpolant_1d_[i];
    }

    builder_return.insert(builder_return.end(), builder_container2d.begin(), builder_container2d.end());
    builder_return.insert(builder_return.end(), builder_container1d.begin(), builder_container1d.end());
    // builder2d.insert(builder2d.end(), builder1d.begin(), builder1d.end());

    log_info("Initialize dNdx for %s", typeid(parametrization_).name());
    Helper::InitializeInterpolation("dNdx",
                                    builder_return,
                                    std::vector<Parametrization*>(1, &parametrization_));
    log_info("Initialize dNdx for %s done!", typeid(parametrization_).name());
}

CrossSectionInterpolant::CrossSectionInterpolant(const CrossSectionInterpolant& cross_section)
    : CrossSection(cross_section)
{
    if (cross_section.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*cross_section.dedx_interpolant_);
    }

    int num_components = cross_section.parametrization_.GetMedium().GetNumComponents();

    dndx_interpolant_1d_.resize(num_components);
    for(InterpolantVec::iterator iter = dndx_interpolant_1d_.begin(); iter != dndx_interpolant_1d_.end(); ++iter)
    {
        if (*iter != NULL)
        {
            dndx_interpolant_1d_.push_back(new Interpolant(**iter));
        }
    }

    dndx_interpolant_2d_.resize(num_components);
    for(InterpolantVec::iterator iter = dndx_interpolant_2d_.begin(); iter != dndx_interpolant_2d_.end(); ++iter)
    {
        if (*iter != NULL)
        {
            dndx_interpolant_2d_.push_back(new Interpolant(**iter));
        }
    }
}

CrossSectionInterpolant::~CrossSectionInterpolant()
{
}

// ------------------------------------------------------------------------- //
// Pulblic methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculatedNdx(double energy)
{
    if(parametrization_.GetMultiplier() <= 0)
    {
        return 0;
    }

    sum_of_rates_ = 0;

    const ComponentVec& components = parametrization_.GetMedium().GetComponents();
    for(size_t i = 0; i < components.size(); ++i)
    {
        prob_for_component_[i] = std::max(dndx_interpolant_1d_.at(i)->Interpolate(energy), 0.);
        sum_of_rates_ += prob_for_component_[i];
    }
    return sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculatedNdx(double energy, double rnd)
{
    if(parametrization_.GetMultiplier() <= 0)
    {
        return 0;
    }

    // The random number will be stored to be able
    // to check if dNdx is already calculated for this random number.
    // This avoids a second calculation in CalculateStochaticLoss
    rnd_ = rnd;

    sum_of_rates_ = 0;

    const ComponentVec& components = parametrization_.GetMedium().GetComponents();
    for(size_t i = 0; i < components.size(); ++i)
    {
        prob_for_component_[i] = std::max(dndx_interpolant_1d_.at(i)->Interpolate(energy), 0.);
        sum_of_rates_ += prob_for_component_.at(i);
    }

    return sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculateStochasticLoss(double energy,  double rnd1, double rnd2)
{
    if(rnd1 != rnd_ )
    {
        CalculatedNdx(energy, rnd1);
    }

    return CalculateStochasticLoss(energy, rnd2);
}


// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculateStochasticLoss(double energy, double rnd1)
{

    double rnd;
    double rsum;

    rnd    =   rnd1*sum_of_rates_;
    rsum    =   0;

    for(size_t i = 0; i < components_.size(); ++i)
    {
        rsum    += prob_for_component_[i];

        if(rsum > rnd)
        {
            parametrization_.SetCurrentComponent(i);
            Parametrization::IntegralLimits limits = parametrization_.GetIntegralLimits(energy);

            if(limits.vUp == limits.vMax)
            {
                return energy * limits.vUp;
            }

            return energy * (limits.vUp * exp(dndx_interpolant_2d_.at(i)->FindLimit(energy, rnd_ * prob_for_component_[i])*log(limits.vMax / limits.vUp)));
        }
    }

    //sometime everything is fine, just the probability for interaction is zero
    bool prob_for_all_comp_is_zero=true;
    for(size_t i = 0; i < components_.size(); ++i)
    {
        parametrization_.SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_.GetIntegralLimits(energy);

        if (limits.vUp != limits.vMax)
            prob_for_all_comp_is_zero = false;
    }

    if (prob_for_all_comp_is_zero)
        return 0;

    //TODO(mario): revert Mon 2017/09/04
    // log_fatal("sum was not initialized correctly");
    return 0;
}

// ------------------------------------------------------------------------- //
// Function needed for interpolation intitialization
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::FunctionToBuildDNdxInterpolant(double energy, int component)
{
    return dndx_interpolant_2d_[component]->Interpolate(energy, 1.);
}

//----------------------------------------------------------------------------//
double CrossSectionInterpolant::FunctionToBuildDNdxInterpolant2D(double energy, double v, int component)
{
    parametrization_.SetCurrentComponent(component);
    Parametrization::IntegralLimits limits = parametrization_.GetIntegralLimits(energy);

    if (limits.vUp == limits.vMax)
    {
        return 0;
    }

    v = limits.vUp * exp(v * log(limits.vMax / limits.vUp));

    return dndx_integral_[component].Integrate(
        limits.vUp, v, boost::bind(&Parametrization::FunctionToDNdxIntegral, &parametrization_, energy, _1), 4);
}

// ------------------------------------------------------------------------- //
// Initialize interpolation
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
// void CrossSectionInterpolant::InitInterpolation(std::string pathname, bool raw)
// {
//     using namespace std;
//
//     bool storing_failed =   false;
//     bool reading_worked =   true;
//
//     // charged anti leptons have the same cross sections like charged leptons
//     // (except of diffractive Bremsstrahlung, where one can analyse the interference term if implemented)
//     // so they use the same interpolation tables
//
//     const ParticleDef& particle_def = parametrization_.GetParticleDef();
//     const Medium& medium = parametrization_.GetMedium();
//     const EnergyCutSettings& cut_settings = parametrization_.GetEnergyCuts();
//     const Parametrization::Definition& param_def = parametrization_.GetDefinition();
//
//     string particlename = parametrization_.GetParticleDef().name;
//
//     if(!pathname.empty())
//     {
//         stringstream filename;
//         filename<<pathname<<"/Brems_dNdx"
//                 <<"_particle"<<particlename
//                 <<"_mass_"<<particle_def.mass
//                 <<"_charge_"<<particle_def.charge
//                 <<"_lifetime_"<<particle_def.lifetime
//                 // <<"_para_"<<parametrization_
//                 <<"_med_"<<medium.GetName()
//                 <<"_"<<medium.GetMassDensity()
//                 <<"_ecut_"<<cut_settings.GetEcut()
//                 <<"_vcut_"<<cut_settings.GetVcut()
//                 <<"_lpm_"<<param_def.lpm_effect_enabled
//                 <<"_multiplier_"<<param_def.multiplier;
//
//         if(!raw)
//             filename<<".txt";
//
//         dndx_interpolant_1d_.resize( medium.GetNumComponents() );
//         dndx_interpolant_2d_.resize( medium.GetNumComponents() );
//
//         if( FileExist(filename.str()) )
//         {
//             log_debug("Bremsstrahlungs parametrisation tables (dNdx) will be read from file:\t%s",filename.str().c_str());
//
//             ifstream input;
//             if(raw)
//             {
//                 input.open(filename.str().c_str(), ios::binary);
//             }
//             else
//             {
//                 input.open(filename.str().c_str());
//             }
//             for(int i=0; i<(medium.GetNumComponents()); i++)
//             {
//                 // component_ = i;
//                 dndx_interpolant_2d_.at(i) = new Interpolant();
//                 dndx_interpolant_1d_.at(i) = new Interpolant();
//                 reading_worked = dndx_interpolant_2d_.at(i)->Load(input,raw);
//                 reading_worked = dndx_interpolant_1d_.at(i)->Load(input,raw);
//
//             }
//             input.close();
//         }
//         if(!FileExist(filename.str()) || !reading_worked )
//         {
//             if(!reading_worked)
//             {
//                 log_info("file %s is corrupted! Write it again!",filename.str().c_str());
//             }
//
//             log_info("Info: Bremsstrahlungs parametrisation tables (dNdx) will be saved to file:\t%s",filename.str().c_str());
//
//             ofstream output;
//
//             if(raw)
//             {
//                 output.open(filename.str().c_str(), ios::binary);
//             }
//             else
//             {
//                 output.open(filename.str().c_str());
//             }
//
//             if(output.good())
//             {
//                 output.precision(16);
//
//                 for(int i=0; i<(medium.GetNumComponents()); i++)
//                 {
//                     // component_ = i;
//
//                     dndx_interpolant_2d_[i] =
//                         new Interpolant(NUM1,
//                                         particle_def.low,
//                                         BIGENERGY,
//                                         NUM1,
//                                         0,
//                                         1,
//                                         boost::bind(&CrossSection::FunctionToBuildDNdxInterpolant2D, this, _1, _2, i),
//                                         param_def.order_of_interpolation,
//                                         false,
//                                         false,
//                                         true,
//                                         param_def.order_of_interpolation,
//                                         false,
//                                         false,
//                                         false,
//                                         param_def.order_of_interpolation,
//                                         true,
//                                         false,
//                                         false);
//                     dndx_interpolant_1d_[i] =
//                         new Interpolant(NUM1,
//                                         particle_def.low,
//                                         BIGENERGY,
//                                         boost::bind(&CrossSection::FunctionToBuildDNdxInterpolant, this, _1, i),
//                                         param_def.order_of_interpolation,
//                                         false,
//                                         false,
//                                         true,
//                                         param_def.order_of_interpolation,
//                                         true,
//                                         false,
//                                         false);
//
//                     dndx_interpolant_2d_.at(i)->Save(output,raw);
//                     dndx_interpolant_1d_.at(i)->Save(output,raw);
//
//                 }
//             }
//             else
//             {
//                 storing_failed  =   true;
//                 log_warn("Can not open file %s for writing! Table will not be stored!",filename.str().c_str());
//             }
//             output.close();
//         }
//     }
//     if(pathname.empty() || storing_failed)
//     {
//         dndx_interpolant_1d_.resize( medium.GetNumComponents() );
//         dndx_interpolant_2d_.resize( medium.GetNumComponents() );
//         for(int i=0; i<(medium.GetNumComponents()); i++)
//         {
//             // component_ = i;
//             dndx_interpolant_2d_.at(i) =
//                 new Interpolant(NUM1,
//                                 particle_def.low,
//                                 BIGENERGY,
//                                 NUM1,
//                                 0,
//                                 1,
//                                 boost::bind(&CrossSection::FunctionToBuildDNdxInterpolant2D, this, _1, _2, i),
//                                 param_def.order_of_interpolation,
//                                 false,
//                                 false,
//                                 true,
//                                 param_def.order_of_interpolation,
//                                 false,
//                                 false,
//                                 false,
//                                 param_def.order_of_interpolation,
//                                 true,
//                                 false,
//                                 false);
//             dndx_interpolant_1d_.at(i) =
//                 new Interpolant(NUM1,
//                                 particle_def.low,
//                                 BIGENERGY,
//                                 boost::bind(&CrossSection::FunctionToBuildDNdxInterpolant, this, _1, i),
//                                 param_def.order_of_interpolation,
//                                 false,
//                                 false,
//                                 true,
//                                 param_def.order_of_interpolation,
//                                 true,
//                                 false,
//                                 false);
//         }
//     }
// }
