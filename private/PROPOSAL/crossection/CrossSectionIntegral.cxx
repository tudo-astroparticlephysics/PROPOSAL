
#include <functional>

#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

using std::bind;
using std::vector;
using std::placeholders::_1;

namespace PROPOSAL {
double transform_relativ_loss(double v_cut, double v_max, double v)
{
    return v_cut * std::exp(v * std::log(v_max / v_cut));
}
} //namespace PROPOSAL;
