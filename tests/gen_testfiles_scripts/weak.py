import pyPROPOSAL as pp
import numpy as np

particle_defs = [
    pp.particle.MuMinusDef.get(),
    pp.particle.TauMinusDef.get(),
    pp.particle.MuPlusDef.get()
]

mediums = [
    pp.medium.Ice(1.0),
    pp.medium.Hydrogen(1.0),
    pp.medium.Uranium(1.0)
]


multiplier = 1.

weak = [
    pp.parametrization.weakinteraction.CooperSarkarMertsch,
]


energies = np.logspace(4, 13, num=10)


interpoldef = pp.InterpolationDef()


def create_tables(dir_name, interpolate=False, **kwargs):

    pp.RandomGenerator.get().set_seed(1234)

    params = weak

    buf = {}

    for key in kwargs:
        if key == "dNdx" and kwargs[key] is True:
            f_dNdx = open(dir_name + "Weak_dNdx{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["dNdx"] = [f_dNdx, [""]]
        if key == "dNdx_rnd" and kwargs[key] is True:
            f_dNdx_rnd = open(dir_name + "Weak_dNdx_rnd{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["dNdx_rnd"] = [f_dNdx_rnd, [""]]
        if key == "stoch" and kwargs[key] is True:
            f_stoch = open(dir_name + "Weak_e{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["stoch"] = [f_stoch, [""]]

    # print(buf)

    for particle in particle_defs:
        for medium in mediums:
                for param in params:
                    if interpolate:
                        param_current = param(
                            particle,
                            medium,
                            multiplier)

                        xsection = pp.crosssection.WeakInterpolant(param_current, interpoldef)
                    else:
                        param_current = param(
                            particle,
                            medium,
                            multiplier
                            )

                        xsection = pp.crosssection.WeakIntegral(param_current)

                    for key in buf:
                        buf[key][1] = [""]

                        for energy in energies:
                            if key == "dNdx":
                                result = [str(xsection.calculate_dNdx(energy))]
                            if key == "dNdx_rnd":
                                rnd = pp.RandomGenerator.get().random_double()
                                result = xsection.calculate_dNdx_rnd(energy, rnd)
                                result = [str(rnd), str(result)]
                            if key == "stoch":
                                rnd1 = pp.RandomGenerator.get().random_double()
                                rnd2 = pp.RandomGenerator.get().random_double()
                                result = xsection.calculate_stochastic_loss(energy, rnd1, rnd2)
                                result = [str(rnd1), str(rnd2), str(result)]

                            buf[key][1].append(particle.name)
                            buf[key][1].append(medium.name)
                            buf[key][1].append(str(multiplier))
                            buf[key][1].append(str(energy))
                            buf[key][1].append(param_current.name)
                            buf[key][1].extend(result)
                            buf[key][1].append("\n")

                        buf[key][0].write("\t".join(buf[key][1]))




def main(dir_name):
    create_tables(dir_name, False, dNdx=True, dNdx_rnd=True, stoch=True)
    create_tables(dir_name, True, dNdx=True, dNdx_rnd=True, stoch=True)

if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
