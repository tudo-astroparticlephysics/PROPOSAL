import os
import proposal as pp
import numpy as np


particle_defs = [
    pp.particle.MuMinusDef(),
    # pp.particle.TauMinusDef(),
    # pp.particle.EMinusDef()
    ]

mediums = [
    pp.medium.Ice(),
    pp.medium.Hydrogen(),
    # pp.medium.Uranium()
    ]

cuts = [
    # pp.EnergyCutSettings(-1, -1),
    pp.EnergyCutSettings(500, 1),
    # pp.EnergyCutSettings(-1, 0.05),
    pp.EnergyCutSettings(500, 0.05)
    ]

statistics = 10
energies = np.logspace(4, 13, num=statistics)  # MeV
initial_energy = 1e12 # MeV
distance = 1000  # cm

position = pp.Cartesian3D(0, 0, 0)
direction = pp.Cartesian3D(1, 0, 0)

pp.InterpolationSettings.tables_path = "~/.local/share/PROPOSAL/tables"

def create_table_continous(dir_name):
    pp.RandomGenerator.get().set_seed(1234)
    with open(dir_name + "Sector_ContinousLoss.txt", "w") as file:
        for particle_def in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    args = {
                        "particle_def": particle_def,
                        "target": medium,
                        "interpolate": True,
                        "cuts": cut
                    }
                    cross = pp.crosssection.make_std_crosssection(**args)
                    int_calc = pp.make_interaction(cross, True)
                    dec_calc = pp.make_decay(cross, particle_def, True)

                    energy = initial_energy

                    buf = [""]
                    buf.append(str(particle_def.name))
                    buf.append(str(medium.name))
                    buf.append(str(cut.ecut))
                    buf.append(str(cut.vcut))
                    buf.append(str(energy))

                    for i in range(statistics):
                        rndd = pp.RandomGenerator.get().random_double()
                        rndi = pp.RandomGenerator.get().random_double()

                        energy = max(
                            dec_calc.energy_decay(energy, rndd, medium.mass_density),
                            int_calc.energy_interaction(energy, rndi)
                        )
                        buf.append(str(energy))

                    buf.append("\n")
                    file.write("\t".join(buf))

def create_table_stochastic(dir_name):
    pp.RandomGenerator.get().set_seed(1234)
    with open(dir_name + "Sector_StochasticLoss.txt", "w") as file:
        for particle_def in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    energy = initial_energy

                    args = {
                        "particle_def": particle_def,
                        "target": medium,
                        "interpolate": True,
                        "cuts": cut
                    }
                    cross = pp.crosssection.make_std_crosssection(**args)
                    int_calc = pp.make_interaction(cross, True)

                    buf = [""]
                    buf.append(str(particle_def.name))
                    buf.append(str(medium.name))
                    buf.append(str(cut.ecut))
                    buf.append(str(cut.vcut))
                    buf.append(str(initial_energy))

                    for i in range(statistics):
                        rnd = pp.RandomGenerator.get().random_double()
                        int_loss = int_calc.sample_loss(
                            energy, int_calc.rates(energy), rnd
                        )
                        energy = energy * (1 - int_loss.v_loss)
                        buf.append(str(energy))
                        buf.append(str(int_loss.type))
                        buf.append(str(rnd))

                    buf.append("\n")
                    file.write("\t".join(buf))

def create_table_energy_displacement(dir_name):
    pp.RandomGenerator.get().set_seed(1234)
    with open(dir_name + "Sector_Energy_Distance.txt", "w") as file:
        for particle_def in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    args = {
                        "particle_def": particle_def,
                        "target": medium,
                        "interpolate": True,
                        "cuts": cut
                    }
                    cross = pp.crosssection.make_std_crosssection(**args)
                    disp_calc = pp.make_displacement(cross, True)

                    for energy in energies:
                        buf = [""]
                        buf.append(str(particle_def.name))
                        buf.append(str(medium.name))
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(energy))

                        for disp in np.geomspace(1e0, 1e7, statistics):
                            try:
                                energy_disp = disp_calc.upper_limit_track_integral(
                                    energy, disp * medium.mass_density)
                            except:
                                energy_disp = disp_calc.lower_limit()
                            buf.append(str(disp))
                            buf.append(str(energy_disp))

                        buf.append("\n")
                        file.write("\t".join(buf))


def main(dir_name):
    print("Create continous testfiles.")
    create_table_continous(dir_name)
    print("Create stochastic testfiles.")
    create_table_stochastic(dir_name)
    print("Create energy displacement testfiles.")
    create_table_energy_displacement(dir_name)

if __name__ == "__main__":

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
