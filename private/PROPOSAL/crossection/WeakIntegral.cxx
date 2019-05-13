
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/WeakIntegral.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

WeakIntegral::WeakIntegral(const WeakInteraction& param)
        : CrossSectionIntegral(DynamicData::WeakInt, param)
{
}

WeakIntegral::WeakIntegral(const WeakIntegral& weak)
        : CrossSectionIntegral(weak)
{
}

WeakIntegral::~WeakIntegral() {}
