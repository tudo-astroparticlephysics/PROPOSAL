#pragma once 

#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL { 
namespace stochastic_deflection { 
    class TsaiParametrizationBremsstrahlungGinneken : public Bremsstrahlung 
                    // ,public DefaultDeflection<TsaiParametrizationBremsstrahlungGinneken> 
                    { 
        static constexpr int n_rnd = 2;
        double mass;

        std::unique_ptr<Parametrization> clone() cont final 
        { 
            return std::unique_ptr<Parametrization>(
                std::make_unique<TsaiParametrizationBremsstrahlungGinneken>(*this));
        }

    public: 
        TsaiParametrizationBremsstrahlungGinneken(const ParticleDef& p_def, const Medium&) : mass(p_def.mass) {};

        size_t RequireRandomNumbers() const noexcept final { return n_rnd; }

        UnitSphericalVector CalculateStochasticDeflection(
            double e_i, double e_f, std::vector<double> const& rnd) const final;
    };
} // namespace stochastic_deflection 
} // namespace PROPOSAL

// A. Van Ginneken. “Energy loss and angular characteristics of high energy elec- tromagnetic processes”. In: Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment 251.1 (1986), pp. 21–39. issn: 0168-9002. doi: 10.1016/0168-9002(86)91146-0.