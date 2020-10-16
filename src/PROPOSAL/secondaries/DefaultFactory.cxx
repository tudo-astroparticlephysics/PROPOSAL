#include "PROPOSAL/secondaries/DefaultFactory.h"

using namespace PROPOSAL;

std::unique_ptr<std::map<InteractionType, secondaries::TCreateMethod>> secondaries::DefaultFactory::secondaries_map = nullptr;
