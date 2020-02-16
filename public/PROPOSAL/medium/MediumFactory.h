#include <memory>
#include <string>

#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

namespace PROPOSAL{
std::shared_ptr<const Medium> CreateMedium(std::string name, double density_correction= 1.0);
}
