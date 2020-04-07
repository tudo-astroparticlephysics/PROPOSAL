#include "PROPOSAL/crossection/CrossSectionBuilder.h"
//#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

CrossSectionBuilder::CrossSectionBuilder(std::string name, const ParticleDef& particle_def, std::shared_ptr<const Medium> medium) : CrossSection(), name(name){
    parametrization_ = new Parametrization_builder(particle_def, medium, 1.0);
    // Set zero return functions as default
    dEdx_function = [](double x)->double {
        (void) x; throw std::logic_error("dEdx_function must be set first");
    };
    dE2dx_function = [](double x)->double {
        (void)x; throw std::logic_error("dE2dx_function must be set first");
    };
    dNdx_function = [](double x)->double {
        (void)x; throw std::logic_error("dNdx_function must be set first");
    };
    dNdx_rnd_function = [](double x1, double x2)->double {
        (void)x1; (void)x2; throw std::logic_error("dNdx_rnd_function must be set first");
    };
    StochasticLoss_function = [](double x1, double x2, double x3)->double {
        (void)x1; (void)x2; (void)x3; throw std::logic_error("StochastcLoss_function must be set first");
    };
    CumulativeCrossSection_function = [](double x1, double x2, double x3)->double {
        (void)x1; (void)x2; (void)x3; throw std::logic_error("CumulativeCrossSection_function must be set first");
    };
}

bool CrossSectionBuilder::compare(const CrossSection& cross) const{
    const CrossSectionBuilder* cross_compare = static_cast<const CrossSectionBuilder*>(&cross);
    if(&dEdx_function != &cross_compare->dEdx_function)
        return false;
    if(&dE2dx_function != &cross_compare->dE2dx_function)
        return false;
    if(&dNdx_function != &cross_compare->dNdx_function)
        return false;
    if(&dNdx_rnd_function != &cross_compare->dNdx_rnd_function)
        return false;
    if(&StochasticLoss_function != &cross_compare->StochasticLoss_function)
        return false;
    if(&CumulativeCrossSection_function != &cross_compare->CumulativeCrossSection_function)
        return false;
}

double CrossSectionBuilder::CalculatedEdx(double energy){
    return dEdx_function(energy);
}

double CrossSectionBuilder::CalculatedE2dx(double energy){
    return dE2dx_function(energy);
}

double CrossSectionBuilder::CalculatedNdx(double energy){
    return dNdx_function(energy);
}

double CrossSectionBuilder::CalculatedNdx(double energy, double rnd){
    return dNdx_rnd_function(energy, rnd);
}

double CrossSectionBuilder::CalculateStochasticLoss(double energy, double rnd1, double rnd2){
    return StochasticLoss_function(energy, rnd1, rnd2);
}

double CrossSectionBuilder::CalculateCumulativeCrossSection(double energy, int component, double v){
    return CumulativeCrossSection_function(energy, component, v);
}

void CrossSectionBuilder::SetdEdx_function(Function func){
    dEdx_function = func;
}

void CrossSectionBuilder::SetdE2dx_function(Function func){
    dE2dx_function = func;
}

void CrossSectionBuilder::SetdNdx_function(Function func){
    dNdx_function = func;
}

void CrossSectionBuilder::SetdNdx_rnd_function(std::function<double(double, double)> func){
    dNdx_rnd_function = func;
}

void CrossSectionBuilder::SetStochasticLoss_function(std::function<double(double, double, double)> func){
    StochasticLoss_function = func;
}

void CrossSectionBuilder::SetCumulativeCrossSection_function(std::function<double(double, double, double)> func){
    CumulativeCrossSection_function = func;
}

double CrossSectionBuilder::CalculateStochasticLoss(double energy, double rnd1){
    (void) energy; (void) rnd1;
    return 0;
}

size_t CrossSectionBuilder::GetHash() const{
    std::size_t seed = 0;
    hash_combine(seed, name);

    return seed;
}

