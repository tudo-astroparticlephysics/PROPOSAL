

#include <chrono>
#include <fstream>
#include <iostream>

#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

int main()
{

    int statistics = 100000;

    int exp_min = 3;
    int exp_max = 13;

    std::vector<std::vector<double> > run_times(exp_max - exp_min);

    for (int i = 0; i < exp_max - exp_min; ++i)
    {
        run_times[i].push_back(pow(10, i + exp_min));
    }

    // ----[ Interpolation ]--------------------------------- //

    InterpolationDef interpolation_def;
    // interpolation_def.path_to_tables =  "../src/resources/tables";
    interpolation_def.do_binary_tables = true;

    // ----[ Sector ]---------------------------------------- //

    std::vector<Sector::Definition> sector_defintions;
    Sector::Definition sec_def;

    sec_def.utility_def.epair_def.lpm_effect = true;
    sec_def.utility_def.brems_def.lpm_effect = true;

    sec_def.location = Sector::ParticleLocation::InsideDetector;

    sec_def.do_continuous_randomization = true;
    sec_def.do_exact_time_calculation   = true;

    sec_def.scattering_model = ScatteringFactory::NoScattering;

    sec_def.cut_settings.SetEcut(500);
    sec_def.cut_settings.SetVcut(0.05);

    sec_def.SetGeometry(Sphere(Vector3D(), 1e18, 0));
    sec_def.SetMedium(AntaresWater());

    sector_defintions.push_back(sec_def);

    // --------------------------------------------------------------------- //
    // Create Propagator
    // --------------------------------------------------------------------- //

    Propagator prop(MuMinusDef::Get(), sector_defintions, Sphere(Vector3D(), 1e18, 0), interpolation_def);

    Particle& particle = prop.GetParticle();
    particle.SetDirection(Vector3D(0, 0, -1));

    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;

    std::ofstream out_file;

    out_file.open("performance_test.txt");

    for (int i = 0; i < exp_max - exp_min; ++i)
    {
        out_file << run_times[i][0];
        out_file << '\t';
    }

    out_file << '\n';

    for (int j = 0; j < statistics; j++)
    {
        std::cout << "particle: " << j << '\n';

        for (int i = 0; i < exp_max - exp_min; ++i)
        {
            t1 = std::chrono::high_resolution_clock::now();

            particle.SetEnergy(pow(10, i + exp_min));
            particle.SetPropagatedDistance(0);
            particle.SetPosition(Vector3D(0, 0, 0));
            particle.SetDirection(Vector3D(0, 0, -1));

            prop.Propagate(10000);

            t2 = std::chrono::high_resolution_clock::now();

            out_file << (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count();
            out_file << '\t';
        }

        out_file << '\n';
    }

    out_file.close();
}
