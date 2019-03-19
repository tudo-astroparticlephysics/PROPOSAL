
#include <fstream>
#include <string>
#include <wordexp.h>
#include <PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h>

using namespace PROPOSAL;

int main()
{
    I3PropagatorServicePROPOSALPtr prop(new I3PropagatorServicePROPOSAL);
    std::ofstream output;
    const char* I3_BUILD = getenv("I3_BUILD");
    std::string s(I3_BUILD);
    output.open(s + "/PROPOSAL/resources/tables/.tables.auto_generated");
}