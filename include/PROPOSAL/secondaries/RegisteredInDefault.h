#pragma once

#include "PROPOSAL/secondaries/DefaultFactory.h"

namespace PROPOSAL {
namespace secondaries {

template <typename T> class RegisteredInDefault {
    protected:
        static bool s_registered;
        RegisteredInDefault() { (void)s_registered; }
};

template <typename T>
bool RegisteredInDefault<T>::s_registered = DefaultFactory::Register<T>(T::type);

} // namespace crosssection
} // namespace PROPOSAL
