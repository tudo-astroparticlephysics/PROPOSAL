#include "PROPOSAL/secondaries/DefaultFactory.h"

using namespace PROPOSAL;

unordered_map<InteractionType, secondaries::TCreateMethod>
    secondaries::DefaultFactory::secondaries_map{};

unique_ptr<secondaries::Parametrization> secondaries::DefaultFactory::Create(
    const InteractionType& type, ParticleDef p, Medium m)
{
    auto it = secondaries_map.find(type);
    if (it != secondaries_map.end())
        return it->second(p, m);
    std::ostringstream s;
    s << "No secondary builder for this interaction type ("
      << Type_Interaction_Name_Map.find(type)->second << ") available.";
    throw logic_error(s.str());
}
