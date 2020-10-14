#include "PROPOSAL/secondaries/DefaultFactory.h"

using namespace PROPOSAL;

std::unordered_map<InteractionType, secondaries::TCreateMethod, InteractionType_hash> secondaries::DefaultFactory::secondaries_map = {};
