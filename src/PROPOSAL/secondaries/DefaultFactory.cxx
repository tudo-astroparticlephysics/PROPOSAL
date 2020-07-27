#include "PROPOSAL/secondaries/DefaultFactory.h"

using namespace PROPOSAL;

std::unordered_map<InteractionType, secondaries::TCreateMethod, InteractionType_hash>& secondaries::DefaultFactory::secondaries_map()
{
  static std::unordered_map<InteractionType, secondaries::TCreateMethod, InteractionType_hash> map;
  return map;
}

unique_ptr<secondaries::Parametrization> secondaries::DefaultFactory::Create(
    const InteractionType& type, ParticleDef p, Medium m)
{
    auto it = secondaries_map().find(type);
    if (it != secondaries_map().end())
        return it->second(p, m);
    std::ostringstream s;
    s << "No secondary builder for this interaction type ("
      << Type_Interaction_Name_Map.find(type)->second << ") available.";
    throw std::logic_error(s.str());
}
