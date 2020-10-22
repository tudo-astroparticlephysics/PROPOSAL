#include "PROPOSAL/propagation_utility/Time.h"

using namespace PROPOSAL;

std::unique_ptr<Interpolant1DBuilder::Definition> Time::interpol_def = nullptr;
