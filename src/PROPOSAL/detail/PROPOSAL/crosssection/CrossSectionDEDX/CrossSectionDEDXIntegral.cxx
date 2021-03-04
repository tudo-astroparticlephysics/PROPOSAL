#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXIntegral.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/crosssection/parametrization/MupairProduction.h"
#include "PROPOSAL/crosssection/parametrization/Photonuclear.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/Constants.h"
#include <memory>

using namespace PROPOSAL;

namespace PROPOSAL {
namespace detail {

    template <typename Param, typename Target>
    dedx_integral_t _define_dedx_integral(Param const& param,
        ParticleDef const& p, Target const& t, EnergyCutSettings const& cut)
    {
        auto param_ptr = std::shared_ptr<Param>(param.clone());
        return [param_ptr, p, t, cut](double E) {
            auto i = Integral();
            auto lim = param_ptr->GetKinematicLimits(p, t, E);
            auto v_cut = cut.GetCut(lim, E);
            auto dEdx = [_param_ptr = param_ptr.get(), &p, &t, E](double v) {
                return _param_ptr->FunctionToDEdxIntegral(p, t, E, v);
            };
            return i.Integrate(lim.v_min, v_cut, dEdx, 2);
        };
    }

    template <typename Param, typename Target>
    dedx_integral_t _define_dedx_integral_log(Param const& param,
        ParticleDef const& p, Target const& t, EnergyCutSettings const& cut)
    {
        auto param_ptr = std::shared_ptr<crosssection::Parametrization<Target>>(param.clone());
        return [param_ptr, p, t, cut](double E) {
            auto i = Integral(IROMB, IMAXS, IPREC);
            auto lim = param_ptr->GetKinematicLimits(p, t, E);
            auto v_cut = cut.GetCut(lim, E);
            auto dEdx = [_param_ptr = param_ptr.get(), &p, &t, E](double v) {
                return _param_ptr->FunctionToDEdxIntegral(p, t, E, v);
            };
            return i.Integrate(lim.v_min, v_cut, dEdx, 4);
        };
    }

    dedx_integral_t define_dedx_integral(
        crosssection::Parametrization<Component> const& param,
        ParticleDef const& p, Component const& c, EnergyCutSettings const& cut)
    {
        return _define_dedx_integral(param, p, c, cut);
    }

    dedx_integral_t define_dedx_integral(
        crosssection::Parametrization<Medium> const& param,
        ParticleDef const& p, Medium const& m, EnergyCutSettings const& cut)
    {
        return _define_dedx_integral(param, p, m, cut);
    }

    dedx_integral_t define_dedx_integral(crosssection::IonizBetheBlochRossi const &param,
        ParticleDef const& p, Medium const& m, EnergyCutSettings const& cut)
    {
        auto param_ptr = std::make_shared<crosssection::IonizBetheBlochRossi>(param);
        return [param_ptr, p, m, cut](double E) {
            auto i = Integral(IROMB, IMAXS, IPREC);
            auto lim = param_ptr->GetKinematicLimits(p, m, E);
            auto v_cut = cut.GetCut(lim, E);
            auto dEdx = [_param_ptr = param_ptr.get(), &p, &m, E](double v) {
                return _param_ptr->FunctionToDEdxIntegral(p, m, E, v);
            };
            return i.Integrate(lim.v_min, v_cut, dEdx, 4) + param_ptr->IonizationLoss(p, m, E) / E;
        };
    }

    template <typename Param>
    dedx_integral_t _define_dedx_integral_ionization(Param const& param,
        ParticleDef const& p, Medium const& m) {
        using param_t = crosssection::Parametrization<Medium>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p, m](double E) {
            return param_ptr->FunctionToDEdxIntegral(p, m, E, 0.) / E;
        };
    }

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerBhabha const& param,
        ParticleDef const& p, Medium const& m, EnergyCutSettings const&)
    {
        return _define_dedx_integral_ionization(param, p, m);
    }

    dedx_integral_t define_dedx_integral(
        crosssection::IonizBergerSeltzerMoller const& param,
        ParticleDef const& p, Medium const& m, EnergyCutSettings const&)
    {
        return _define_dedx_integral_ionization(param, p, m);
    }

    dedx_integral_t define_dedx_integral(
        crosssection::ComptonKleinNishina const& param, ParticleDef const& p,
        Component const& c, EnergyCutSettings const& cut)
    {
        using param_t = crosssection::Parametrization<Component>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p, c, cut](double E) {
            // Integrate with the substitution t = ln(1-v) to avoid
            // numerical problems
            auto i = Integral();
            auto lim = param_ptr->GetKinematicLimits(p, c, E);
            auto v_cut = cut.GetCut(lim, E);
            double t_min = std::log(1. - v_cut);
            double t_max = std::log(1. - lim.v_min);
            auto dEdx = [ptr = param_ptr.get(), p, c, E](double t) {
                return std::exp(t)
                    * ptr->FunctionToDEdxIntegral(p, c, E, 1 - std::exp(t));
            };
            return i.Integrate(t_min, t_max, dEdx, 2);
        };
    }

    dedx_integral_t define_dedx_integral(
        crosssection::EpairProduction const& param, ParticleDef const& p,
        Component const& c, EnergyCutSettings const& cut)
    {
        using param_t = crosssection::Parametrization<Component>;
        auto param_ptr = std::shared_ptr<param_t>(param.clone());
        return [param_ptr, p, c, cut](double E) {
            auto lim = param_ptr->GetKinematicLimits(p, c, E);
            auto i = Integral();
            auto v_cut = cut.GetCut(lim, E);
            auto dEdx = [_param_ptr = param_ptr.get(), &p, &c, E](double v) {
                return _param_ptr->FunctionToDEdxIntegral(p, c, E, v);
            };
            auto dEdx_reverse= [_param_ptr = param_ptr.get(), &p, &c, E](double v) {
                return (1 - v) * _param_ptr->DifferentialCrossSection(p, c, E, 1 - v);
            };
            auto r1 = 0.8;
            auto rUp = v_cut * (1 - HALF_PRECISION);
            auto rflag = false;
            if (r1 < rUp)
                if (2 * param_ptr->FunctionToDEdxIntegral(p, c, E, r1) <
                     param_ptr->FunctionToDEdxIntegral(p, c, E, rUp))
                    rflag = true;
            if (rflag) {
                if (r1 > v_cut)
                    r1 = v_cut;
                if (r1 < lim.v_min)
                    r1 = lim.v_min;
                auto r2 = std::max(1 - v_cut, COMPUTER_PRECISION);
                if (r2 > 1 - r1)
                    r2 = 1 - r1;
                return (i.Integrate(lim.v_min, r1, dEdx, 4) + i.Integrate(1 - v_cut, r2, dEdx_reverse, 2) + i.Integrate(r2, 1 - r1, dEdx_reverse, 4));
            } else {
                return i.Integrate(lim.v_min, v_cut, dEdx, 4);
            }
        };
    }

    dedx_integral_t define_dedx_integral(crosssection::Photonuclear const& param,
        ParticleDef const& p, Component const& c, EnergyCutSettings const& cut)
    {
        return _define_dedx_integral_log(param, p, c, cut);
    }

    dedx_integral_t define_dedx_integral(crosssection::MupairProduction const& param,
        ParticleDef const& p, Component const& c, EnergyCutSettings const& cut)
    {
        return _define_dedx_integral_log(param, p, c, cut);
    }
}
}

double CrossSectionDEDXIntegral::Calculate(double E) const
{
    return dedx_integral(E) * E;
}
