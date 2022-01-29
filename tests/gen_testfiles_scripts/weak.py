import os
import proposal as pp
import numpy as np

particle_defs = [
    pp.particle.MuMinusDef(),
    pp.particle.TauMinusDef(),
    pp.particle.MuPlusDef()
]

mediums = [
    pp.medium.Ice(),
    pp.medium.Hydrogen(),
    pp.medium.Uranium()
]

multiplier = 1.

params = [
    pp.parametrization.weakinteraction.CooperSarkarMertsch(),
]

secondaries = [
    pp.secondaries.WeakCooperSarkarMertsch,
]

energies = np.logspace(4, 13, num=10)


def create_tables(dir_name, **kwargs):

    pp.RandomGenerator.get().set_seed(1234)

    buf = {}

    for key in kwargs:
        if key == "dNdx" and kwargs[key] is True:
            f_dNdx = open(dir_name + "Weak_dNdx.txt", "w")
            buf["dNdx"] = [f_dNdx, [""]]
        if key == "stoch" and kwargs[key] is True:
            f_stoch = open(dir_name + "Weak_e.txt", "w")
            buf["stoch"] = [f_stoch, [""]]

    for particle in particle_defs:
        for medium in mediums:
            for param, secondary in zip(params, secondaries):

                args = {
                    "parametrization": param,
                    "particle_def": particle,
                    "target": medium,
                    "cuts": None,
                    "interpolate": False
                }

                xsection = pp.crosssection.make_crosssection(**args)
                sec_calculator = secondary(particle, medium)

                for key in buf:
                    buf[key][1] = [""]

                    for energy in energies:
                        if key == "dNdx":
                            result = [str(xsection.calculate_dNdx(energy))]
                        if key == "stoch":
                            rnd1 = pp.RandomGenerator.get().random_double()
                            rnd2 = pp.RandomGenerator.get().random_double()

                            components = medium.components
                            comp = components[int(rnd2*len(components))]

                            result = sec_calculator.calculate_relative_loss(energy, rnd1, comp)
                            result = [str(rnd1), str(rnd2), str(result)]

                        buf[key][1].append(particle.name)
                        buf[key][1].append(medium.name)
                        buf[key][1].append(str(multiplier))
                        buf[key][1].append(str(energy))
                        buf[key][1].append(xsection.param_name)
                        buf[key][1].extend(result)
                        buf[key][1].append("\n")

                    buf[key][0].write("\t".join(buf[key][1]))


def main(dir_name):
    create_tables(dir_name, dNdx=True, stoch=True)

if __name__ == "__main__":

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
