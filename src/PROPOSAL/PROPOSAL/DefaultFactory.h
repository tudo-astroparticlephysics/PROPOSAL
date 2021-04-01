#pragma once

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include <iostream>

namespace PROPOSAL {

template <typename T> class DefaultFactory {
public:
    using TCreate = std::unique_ptr<T> (*)(ParticleDef const&, Medium const&);
    using m_t = std::map<InteractionType, TCreate>;

private:
    static std::unique_ptr<m_t> m;

public:
    DefaultFactory() = delete;

    template <typename T1> static bool Register(InteractionType t)
    {
        if (!m)
            m = std::make_unique<std::map<InteractionType, TCreate>>();

        auto it = m->find(t);
        if (it != m->end())
            return false;
        (*m)[t] = [](ParticleDef const& pa, Medium const& me) {
            return std::unique_ptr<T>(std::make_unique<T1>(pa, me));
        };
        return true;
    }

    template <typename... Args>
    static std::unique_ptr<T> Create(InteractionType t, Args... args)
    {
        if(!m)
            throw std::logic_error("No default registered at all.");

        auto it = m->find(t);
        if (it != m->end())
            return it->second(args...);
        std::ostringstream s;
        s << "No builder for this interaction type ("
          << Type_Interaction_Name_Map.find(t)->second << ") available.";
        throw std::logic_error(s.str());
    }
};

template <typename T>
std::unique_ptr<typename DefaultFactory<T>::m_t> DefaultFactory<T>::m = nullptr;

} // namespace PROPOSAL
