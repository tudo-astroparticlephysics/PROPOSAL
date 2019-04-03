
#include <cmath>

#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"


using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

WeakInteraction::WeakInteraction(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 double multiplier)
        : Parametrization(particle_def, medium, cuts, multiplier)
{
}

WeakInteraction::WeakInteraction(const WeakInteraction& param)
        : Parametrization(param)
{
}

WeakInteraction::~WeakInteraction() {}

bool WeakInteraction::compare(const Parametrization& parametrization) const
{
    return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

Parametrization::IntegralLimits WeakInteraction::GetIntegralLimits(double energy)
{
    IntegralLimits limits;

    double aux = (MP + MN) / 2; //for isoscalar targets
    aux = 2 * energy * aux + pow(aux, 2);

    limits.vMin = 1e6 / (aux ); //q^2_min = 1e6 MeV for experimental reasons
    limits.vUp = limits.vMin; //treat interaction as fully stochastic
    limits.vMax = 1;

    return limits;
}

// ------------------------------------------------------------------------- //
// Specific implementations
// ------------------------------------------------------------------------- //

WeakCooperSarkarMertsch::WeakCooperSarkarMertsch(const ParticleDef& particle_def,
                                                 const Medium& medium,
                                                 const EnergyCutSettings& cuts,
                                                 double multiplier)
        : WeakInteraction(particle_def, medium, cuts, multiplier)
        , interpolant_(this->medium_->GetNumComponents(), NULL)
{

    if(particle_def==EMinusDef::Get()||particle_def==MuMinusDef::Get()||particle_def==TauMinusDef::Get()){
        read_table(energies_, y_, dsigma_, false);
    }
    else if(particle_def==EPlusDef::Get()||particle_def==MuPlusDef::Get()||particle_def==TauPlusDef::Get()){
        read_table(energies_, y_, dsigma_, true);
    }
    else{
        log_fatal("Weak interaction table_read: Particle to propagate is not a SM charged lepton");
    }

    std::vector<Interpolant2DBuilder_array_as> builder2d(this->components_.size());
    Helper::InterpolantBuilderContainer builder_container2d(this->components_.size());

    for (unsigned int i = 0; i < this->components_.size(); ++i)
    {
        builder2d[i]
                .Setx1(energies_)
                .Setx2(y_)
                .Sety(dsigma_)
                .SetRomberg1(IROMB)
                .SetRational1(false)
                .SetRelative1(false)
                .SetRomberg2(IROMB)
                .SetRational2(false)
                .SetRelative2(false);

        builder_container2d[i].first  = &builder2d[i];
        builder_container2d[i].second = &interpolant_[i];
    }

    Helper::InitializeInterpolation("WeakInt", builder_container2d, std::vector<Parametrization*>(1, this), InterpolationDef());

}

WeakCooperSarkarMertsch::WeakCooperSarkarMertsch(const WeakCooperSarkarMertsch& param)
        : WeakInteraction(param)
        , interpolant_()
{
    interpolant_.resize(param.interpolant_.size());

    for (unsigned int i = 0; i < param.interpolant_.size(); ++i)
    {
        interpolant_[i] = new Interpolant(*param.interpolant_[i]);
    }
}

WeakCooperSarkarMertsch::~WeakCooperSarkarMertsch()
{
    for (std::vector<Interpolant*>::const_iterator iter = interpolant_.begin(); iter != interpolant_.end(); ++iter)
    {
        delete *iter;
    }

    interpolant_.clear();
}

bool WeakCooperSarkarMertsch::compare(const Parametrization& parametrization) const
{
    const WeakCooperSarkarMertsch* weak = static_cast<const WeakCooperSarkarMertsch*>(&parametrization);

    if (interpolant_.size() != weak->interpolant_.size())
        return false;

    for (unsigned int i = 0; i < interpolant_.size(); ++i)
    {
        if (*interpolant_[i] != *weak->interpolant_[i])
            return false;
    }

    return WeakInteraction::compare(parametrization);
}

double WeakCooperSarkarMertsch::DifferentialCrossSection(double energy, double v)
{
    return multiplier_ * medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule() * 1e-36 * std::max(interpolant_.at(this->component_index_)->InterpolateArray(std::log10(energy), v), 0.0);
    //factor 1e-36: conversion from pb to cm^2
}

void WeakCooperSarkarMertsch::read_table(std::vector<double> &energies, std::vector<std::vector<double> >  &y,
                                         std::vector<std::vector<double> > &dsigma, bool antiparticle) {
    int energylen = 111;
    int ylen = 100;
    double aux, aux2;

    //resize arrays
    energies.resize(energylen);
    y.resize(energylen);
    dsigma.resize(energylen);
    for(int i = 0; i < energylen; i++){
        dsigma[i].resize(ylen);
        y[i].resize(ylen);
    }

    //read data
    std::string pathname;
    if(antiparticle == false){
        pathname = "resources/dsigmady_nubar_CC_iso_NLO_HERAPDF15NLO_EIG.dat";
    }
    else{
        pathname = "resources/dsigmady_nu_CC_iso_NLO_HERAPDF15NLO_EIG.dat";
    }

    std::ifstream file(Helper::ResolvePath(pathname));

    if (!file.is_open())
    {
        log_fatal("Can not find file for weak interaction");
    }

    for(int i = 0; i < energylen; i++){
        energies[i] = 4. + i * 0.1;
        for(int k = 0; k<ylen; ++k){
            file >> aux;
            file >> aux2;
            y[i][k] = aux;
            dsigma[i][k] = aux2;
        }
    }

    file.close();
}

const std::string WeakCooperSarkarMertsch::name_ = "WeakCooperSarkarMertsch";

