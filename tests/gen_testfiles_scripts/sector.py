import pyPROPOSAL as pp
import numpy as np


particle_defs = [
    pp.MuMinusDef.get(),
    pp.TauMinusDef.get(),
    pp.EMinusDef.get()
]

mediums = [
    pp.Medium.Ice(1.0),
    pp.Medium.Hydrogen(1.0),
    pp.Medium.Uranium(1.0)
]

cuts = [
    pp.EnergyCutSettings(-1, -1),
    pp.EnergyCutSettings(500, -1),
    pp.EnergyCutSettings(-1, 0.05),
    pp.EnergyCutSettings(500, 0.05)
]

energies = np.logspace(4, 13, num=10)  # MeV
distance = 1000  # cm

position = pp.Vector3D(0, 0, 0)
direction = pp.Vector3D(1, 0, 0)

# stopping_decay = True
# do_continuous_randomization = True
# do_exact_time_calculation = True
interpoldef = pp.InterpolationDef()

pp.RandomGenerator.get().set_seed(1234)


def create_table_propagate(dir_name):

    with open(dir_name + "Sector_propagate.txt", "a") as f:

        for particle_def in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    sector_def = pp.SectorDefinition()
                    sector_def.cut_settings = cut
                    sector_def.medium = medium

                    particle = pp.Particle(particle_def)
                    sector = pp.Sector(particle, sector_def, interpoldef)

                    for energy in energies:
                        sector.particle.energy = energy
                        sector.particle.position = position
                        sector.particle.direction = direction
                        energy_final = sector.propagate(distance)

                        buf = [""]

                        buf.append(particle_def.name)
                        buf.append(medium.name)
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(energy))
                        buf.append(str(energy_final))
                        buf.append(str(distance))
                        buf.append("\n")

                        f.write("\t".join(buf))


def main(dir_name):
    create_table_propagate(dir_name)


if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    try:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))
    except OSError:
        print("Directory {} already exists".format(dir_name))

    main(dir_name)
