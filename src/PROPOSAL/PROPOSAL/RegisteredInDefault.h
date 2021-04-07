#pragma once

#include "PROPOSAL/DefaultFactory.h"

#include <iostream>

namespace PROPOSAL {
    //!
    //! Storage class to register interaction dependent default secondaries
    //! calculators. Storage is limited that only one secondary calculator type
    //! per interaction type can be stored. If different particle make use of
    //! different types of secondries calculator (not objects) for the same
    //! interaction type, the default register has to be modiefied.
    //!
    template <typename T1, typename T2>
    class RegisteredInDefault {
    protected:
        static bool s_registered;
        RegisteredInDefault() { (void)s_registered; }
    };

    //!
    //! Flag that a default secondary builder for these type of interaction is
    //! registered.
    //!
    template <typename T1, typename T2>
    bool RegisteredInDefault<T1, T2>::s_registered
        = DefaultFactory<T1>::template Register<T2>(T2::type);

} // namespace PROPOSAL
