
#include <functional>

#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

using std::bind;
using std::vector;
using std::placeholders::_1;

using namespace PROPOSAL;

CrossSectionIntegral::CrossSectionIntegral(unique_ptr<Parametrization>&& param,
    shared_ptr<const EnergyCutSettings> cut)
    : CrossSection(std::forward<unique_ptr<Parametrization>>(param), cut)
{
    std::cout << "Number of components: "
              << parametrization_->GetComponents().size() << std::endl;
    for (auto comp : parametrization_->GetComponents()) {
        /*     parametrization_->SetCurrentComponent(comp); */
        /*     dndx_integral_.emplace_back( */
        /*         bind(&CrossSectionIntegral::dndx_integral, this, _1, _2)); */
        dedx_integral_.emplace_back(
            bind(&CrossSectionIntegral::dedx_integral, this, _1));
        /*     de2dx_integral_.emplace_back( */
        /*         bind(&CrossSectionIntegral::de2dx_integral, this, _1)); */
    }
}

double CrossSectionIntegral::dndx_integral(double energy, double rnd)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_max = parametrization_->GetKinematicLimits(energy).vMax;

    integral_.IntegrateWithRandomRatio(v_cut, v_max,
        bind(&Parametrization::FunctionToDNdxIntegral, parametrization_.get(),
            energy, _1),
        4, rnd);

    return integral_.GetUpperLimit();
}

double CrossSectionIntegral::dedx_integral(double energy)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_min = parametrization_->GetKinematicLimits(energy).vMin;

    return integral_.Integrate(v_min, v_cut,
        bind(&Parametrization::FunctionToDEdxIntegral, parametrization_.get(),
            energy, _1),
        2);
}

double CrossSectionIntegral::de2dx_integral(double energy)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_min = parametrization_->GetKinematicLimits(energy).vMin;

    return integral_.Integrate(v_min, v_cut,
        bind(&Parametrization::FunctionToDE2dxIntegral, parametrization_.get(),
            energy, _1),
        2);
}

double CrossSectionIntegral::CalculatedNdx(double energy, double rnd)
{
    vector<double> rates;
    for (auto& dndx : dndx_integral_)
        rates.push_back(dndx(energy, 1.));

    auto total_rate = accumulate(rates.begin(), rates.end(), 0);

    size_t nth_component = 0;
    for (const auto& rate : rates) {
        total_rate -= rate / rnd;
        if (total_rate < 0)
            return dndx_integral_.at(nth_component)(energy, total_rate);
        nth_component += 1;
    }

    return total_rate;
}

double CrossSectionIntegral::CalculatedEdx(double energy)
{
    auto integrate_and_sum
        = [energy](double sum, function<double(double)> func) {
              return sum + func(energy);
          };

    std::cout << "size of dedx_integral_ : " << dedx_integral_.size()
              << std::endl;
    std::cout << dedx_integral_[0](energy) << std::endl;

    auto sum = accumulate(
        dedx_integral_.begin(), dedx_integral_.end(), 0, integrate_and_sum);

    return energy * sum;
}

double CrossSectionIntegral::CalculatedE2dx(double energy)
{
    auto integrate_and_sum
        = [energy](double sum, function<double(double)> func) {
              return sum + func(energy);
          };

    auto sum = accumulate(
        de2dx_integral_.begin(), de2dx_integral_.end(), 0, integrate_and_sum);

    return energy * energy * sum;
}
