import os
import proposal as pp
import numpy as np

scat_names = [
    "Moliere",
    "HighlandIntegral",
    "Highland",
]

particle_defs = [
    pp.particle.MuMinusDef(),
    pp.particle.TauMinusDef(),
    pp.particle.EMinusDef()
]

mediums = [
    pp.medium.Ice(),
    pp.medium.Hydrogen(),
    pp.medium.Uranium()
]

cuts = [
    pp.EnergyCutSettings(np.inf, 1),
    pp.EnergyCutSettings(500, 1),
    pp.EnergyCutSettings(np.inf, 0.05),
    pp.EnergyCutSettings(500, 0.05)
]

energies = np.logspace(4, 13, num=10)  # MeV
energies_out = np.logspace(3, 12, num=10)
distance = 1000  # cm
position_init = pp.Cartesian3D(0, 0, 0)
direction_init = pp.Cartesian3D(1, 0, 0)

pp.RandomGenerator.get().set_seed(1234)

def create_table(dir_name):

    with open(dir_name + "Scattering_scatter.txt", "w") as file:

        for particle_def in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for scat_name in scat_names:

                        if scat_name == "HighlandIntegral" and particle_def.name == "MuMinus":
                            args = {
                                "particle_def": particle_def,
                                "target": medium,
                                "interpolate": False,
                                "cuts": cut
                            }
                            cross = pp.crosssection.make_std_crosssection(**args)
                            scattering = pp.make_multiple_scattering(
                                scat_name, particle_def, medium, cross, False)
                        elif scat_name != "HighlandIntegral":
                            scattering = pp.make_multiple_scattering(
                                scat_name, particle_def, medium)
                        else:
                            continue

                        buf = [""]
                        for jdx, energy in enumerate(energies):

                            rnd1 = pp.RandomGenerator.get().random_double()
                            rnd2 = pp.RandomGenerator.get().random_double()
                            rnd3 = pp.RandomGenerator.get().random_double()
                            rnd4 = pp.RandomGenerator.get().random_double()
                            coords = scattering.scatter(distance*medium.mass_density,
                                                            energy,
                                                            energies_out[jdx],
                                                            [rnd1, rnd2, rnd3, rnd4])

                            directions = pp.scattering.scatter_initial_direction(
                                direction_init, coords)

                            posi = position_init + distance * directions[0]
                            dire = directions[1]

                            buf.append(particle_def.name)
                            buf.append(medium.name)
                            buf.append(scat_name)
                            buf.append(str(cut.ecut))
                            buf.append(str(cut.vcut))
                            buf.append(str(energy))
                            buf.append(str(energies_out[jdx]))
                            buf.append(str(distance))
                            buf.append(str(rnd1))
                            buf.append(str(rnd2))
                            buf.append(str(rnd3))
                            buf.append(str(rnd4))
                            buf.append(str(posi.x))
                            buf.append(str(posi.y))
                            buf.append(str(posi.z))
                            buf.append(str(dire.x))
                            buf.append(str(dire.y))
                            buf.append(str(dire.z))
                            buf.append("\n")

                        file.write("\t".join(buf))


def main(dir_name):
    create_table(dir_name)

if __name__ == "__main__":

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
