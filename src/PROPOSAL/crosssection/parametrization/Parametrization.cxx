
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include <sstream>

using namespace PROPOSAL;
using std::make_tuple;
using std::string;

crosssection::Parametrization::Parametrization(
    InteractionType _type, const string& _name)
    : interaction_type(_type)
    , name(_name)
    , hash(0)
{
    hash_combine(hash, interaction_type, name);
}
