import os
import proposal as pp
import numpy as np

particle_defs = [
    pp.particle.GammaDef(),
]

mediums = [
    pp.medium.Air(),
    pp.medium.Ice(),
    pp.medium.Uranium()
]

multiplier = 1.

params = [
    pp.parametrization.photomupair.BurkhardtKelnerKokoulin(),
]

secondary_params = [
    pp.secondaries.PhotoMuPairProductionBurkhardtKelnerKokoulin,
]

secondary_params_names = [
    "BurkhardtKelnerKokoulin"
]

energies = np.logspace(1, 13, num=13)


def create_tables(dir_name, **kwargs):

    pp.RandomGenerator.get().set_seed(1234)

    buf = {}

    for key in kwargs:
        if key == "dNdx" and kwargs[key] is True:
            f_dNdx = open(dir_name + "PhotoMuPair_dNdx.txt", "w")
            buf["dNdx"] = [f_dNdx, [""]]

    for particle in particle_defs:
        for medium in mediums:
            for param in params:
                args = {
                    "parametrization": param,
                    "particle_def": particle,
                    "target": medium,
                    "cuts": None,
                    "interpolate": False
                }

                xsection = pp.crosssection.make_crosssection(**args)

                for key in buf:
                    buf[key][1] = [""]

                    for energy in energies:
                        if key == "dNdx":
                            result = [str(xsection.calculate_dNdx(energy))]

                        buf[key][1].append(particle.name)
                        buf[key][1].append(medium.name)
                        buf[key][1].append(str(multiplier))
                        buf[key][1].append(str(energy))
                        buf[key][1].append(xsection.param_name)
                        buf[key][1].extend(result)
                        buf[key][1].append("\n")

                    buf[key][0].write("\t".join(buf[key][1]))


def create_tables_x(dir_name, **kwargs):

    pp.RandomGenerator.get().set_seed(1234)

    buf = {}

    f_x = open(dir_name + "PhotoMuPair_x.txt", "w")
    buf["x"] = [f_x, [""]]

    for particle in particle_defs:
        for medium in mediums:
            for param, param_name in zip(secondary_params, secondary_params_names):
                args = {
                    "particle_def": particle,
                    "medium": medium,
                }

                secondary_calc = param(**args)

                for key in buf:
                    buf[key][1] = [""]

                    for energy in energies:
                        if (energy <= 1000):
                            continue
                        rnd1 = pp.RandomGenerator.get().random_double()
                        rnd2 = pp.RandomGenerator.get().random_double()

                        components = medium.components
                        comp = components[int(rnd2*len(components))]

                        x = str(secondary_calc.calculate_x(energy, rnd1, comp))

                        buf[key][1].append(particle.name)
                        buf[key][1].append(medium.name)
                        buf[key][1].append(str(energy))
                        buf[key][1].append(param_name)
                        buf[key][1].append(x)
                        buf[key][1].append(str(rnd1))
                        buf[key][1].append(str(rnd2))
                        buf[key][1].append("\n")

                    buf[key][0].write("\t".join(buf[key][1]))

def main(dir_name):
    create_tables(dir_name, dNdx=True)
    create_tables_x(dir_name)

if __name__ == "__main__":

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
