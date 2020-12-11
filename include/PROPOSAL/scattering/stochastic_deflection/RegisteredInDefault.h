#pragma once

#include "PROPOSAL/DefaultFactory.h"
#include "PROPOSAL/scattering/stochastic_deflection/parametrization/Parametrization.h"

namespace PROPOSAL {
namespace stochastic_deflection {

    //!
    //! Storage class to register interaction dependent default secondaries
    //! calculators. Storage is limited that only one secondary calculator type
    //! per interaction type can be stored. If different particle make use of
    //! different types of secondries calculator (not objects) for the same
    //! interaction type, the default register has to be modiefied.
    //!
    template <typename T> class RegisteredInDefault {
    protected:
        static bool s_registered;
        RegisteredInDefault() { (void)s_registered; }
    };

    //!
    //! Flag that a default secondary builder for these type of interaction is
    //! registered.
    //!
    template <typename T>
    bool RegisteredInDefault<T>::s_registered
        = DefaultFactory<Parametrization>::Register<T>(T::type);

} // namespace crosssection
} // namespace PROPOSAL
