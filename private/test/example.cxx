#include "PROPOSAL/PROPOSAL.h"
#include <vector>
#include <iostream>

using namespace PROPOSAL;

int main(int argc, char* argv[])
{

    MuMinusDef muon;
    std::vector<std::unique_ptr<CrossSection>> crosssections;
    auto cut = std::make_shared<EnergyCutSettings>(500,0.05,true);
    auto medium = std::make_shared<const FrejusRock>();
    InterpolationDef inter_def;

    BremsKelnerKokoulinPetrukhin brems_param(muon, medium, 1.0, true);

    IonizBetheBlochRossi ioniz_param(muon, medium, cut, 1.0);

    /*BremsInterpolant brems_integral(brems_param, cut, inter_def); */
    IonizInterpolant brems_integral(ioniz_param, cut, inter_def);
    /* auto dedx = brems_integral.CalculatedEdx(1.e9); */
    /* std::cout << dedx << std::endl; */

    return 0;
}
