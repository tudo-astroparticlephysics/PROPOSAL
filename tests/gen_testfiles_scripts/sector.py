import proposal as pp
import numpy as np


particle_defs = [
    pp.particle.MuMinusDef(),
    pp.particle.TauMinusDef(),
    # pp.particle.EMinusDef()
    ]

mediums = [
    pp.medium.Ice(1.0),
    pp.medium.Hydrogen(1.0),
    # pp.medium.Uranium(1.0)
    ]

cuts = [
    # pp.EnergyCutSettings(-1, -1),
    pp.EnergyCutSettings(500, -1),
    # pp.EnergyCutSettings(-1, 0.05),
    pp.EnergyCutSettings(500, 0.05)
    ]

statistics = 10
energies = np.logspace(4, 13, num=statistics)  # MeV
initial_energy = 1e12 # MeV
distance = 1000  # cm

position = pp.Vector3D(0, 0, 0)
direction = pp.Vector3D(1, 0, 0)

# stopping_decay = True
# do_continuous_randomization = True
# do_exact_time_calculation = True
interpoldef = pp.InterpolationDef()
interpoldef.path_to_tables = "~/.local/share/PROPOSAL/tables"
interpoldef.path_to_tables_readonly = "~/.local/share/PROPOSAL/tables"


def create_table_continous(dir_name):
    pp.RandomGenerator.get().set_seed(1234)
    with open(dir_name + "Sector_ContinousLoss.txt", "w") as file:
        for particle_def in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    energy = initial_energy
                    sector_def = pp.SectorDefinition()
                    sector_def.cut_settings = cut
                    sector_def.medium = medium

                    sector = pp.Sector(particle_def, sector_def, interpoldef)

                    buf = [""]
                    buf.append(str(particle_def.name))
                    buf.append(str(medium.name))
                    buf.append(str(cut.ecut))
                    buf.append(str(cut.vcut))
                    buf.append(str(energy))

                    for i in range(statistics):
                        energy = max(
                            sector.energy_decay(energy,
                                pp.RandomGenerator.get().random_double()),
                            sector.energy_interaction(energy,
                                pp.RandomGenerator.get().random_double())
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
                    sector_def = pp.SectorDefinition()
                    sector_def.cut_settings = cut
                    sector_def.medium = medium

                    sector = pp.Sector(particle_def, sector_def, interpoldef)

                    buf = [""]
                    buf.append(str(particle_def.name))
                    buf.append(str(medium.name))
                    buf.append(str(cut.ecut))
                    buf.append(str(cut.vcut))
                    buf.append(str(initial_energy))

                    for i in range(statistics):
                        loss, interaction_type = sector.make_stochastic_loss(energy)
                        energy -= loss
                        buf.append(str(energy))
                        buf.append(str(interaction_type))
                        buf.append(str(pp.RandomGenerator.get().random_double()))

                    buf.append("\n")
                    file.write("\t".join(buf))

def create_table_energy_displacement(dir_name):
    pp.RandomGenerator.get().set_seed(1234)
    with open(dir_name + "Sector_Energy_Distance.txt", "w") as file:
        for particle_def in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    sector_def = pp.SectorDefinition()
                    sector_def.cut_settings = cut
                    sector_def.medium = medium

                    sector = pp.Sector(particle_def, sector_def, interpoldef)

                    for energy in energies:
                        buf = [""]
                        buf.append(str(particle_def.name))
                        buf.append(str(medium.name))
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(energy))

                        for disp in np.geomspace(1e0, 1e7, statistics):
                            energy_disp = sector.energy_distance(energy, disp)
                            buf.append(str(disp))
                            buf.append(str(energy_disp))

                        buf.append("\n")
                        file.write("\t".join(buf))

def create_table_propagate(dir_name):
    pp.RandomGenerator.get().set_seed(1234)
    with open(dir_name + "Sector_Propagate.txt", "w") as file:
        for particle_def in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    sector_def = pp.SectorDefinition()
                    sector_def.cut_settings = cut
                    sector_def.medium = medium

                    sector = pp.Sector(particle_def, sector_def, interpoldef)

                    for energy in energies:
                        buf = [""]
                        buf.append(str(particle_def.name))
                        buf.append(str(medium.name))
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(energy))

                        p_condition = pp.particle.DynamicData(0)
                        p_condition.position = pp.Vector3D(0,0,0)
                        p_condition.direction = pp.Vector3D(0,0,-1)
                        p_condition.propagated_distance = 0
                        p_condition.energy = energy
                        sec = sector.propagate(p_condition, 1000, 0)

                        buf.append(str(sec.number_of_particles))
                        if(sec.particles[-1].propagated_distance < 999.99):
                            buf.append(str(-sec.particles[-1].propagated_distance))
                        else:
                            buf.append(str(sec.particles[-1].energy))
                        buf.append("\n")
                        file.write("\t".join(buf))


def main(dir_name):
    print("Create continous testfiles.")
    create_table_continous(dir_name)
    print("Create stochastic testfiles.")
    create_table_stochastic(dir_name)
    print("Create energy displacement testfiles.")
    create_table_energy_displacement(dir_name)
    print("Create propagate testfiles.")
    create_table_propagate(dir_name)


if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
