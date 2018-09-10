import pyPROPOSAL as pp
import numpy as np

particle_defs = [
    pp.particle.MuMinusDef.get(),
    pp.particle.TauMinusDef.get(),
    pp.particle.EMinusDef.get()
]

mediums = [
    pp.medium.Ice(1.0),
    pp.medium.Hydrogen(1.0),
    pp.medium.Uranium(1.0)
]

cuts = [
    pp.EnergyCutSettings(-1, -1),
    pp.EnergyCutSettings(500, -1),
    pp.EnergyCutSettings(-1, 0.05),
    pp.EnergyCutSettings(500, 0.05)
]

multiplier = 1.

epair = [
    pp.parametrization.pairproduction.KelnerKokoulinPetrukhin,
    pp.parametrization.pairproduction.SandrockSoedingreksoRhode,
]

epair_interpol = [
    pp.parametrization.pairproduction.KelnerKokoulinPetrukhinInterpolant,
    pp.parametrization.pairproduction.SandrockSoedingreksoRhodeInterpolant,
]

lpms = [0, 1]

energies = np.logspace(4, 13, num=10)

interpoldef = pp.InterpolationDef()


def create_tables(dir_name, interpolate=False, **kwargs):

    pp.RandomGenerator.get().set_seed(1234)

    if interpolate:
        params = epair_interpol
    else:
        params = epair

    buf = {}

    for key in kwargs:
        if key == "dEdx" and kwargs[key] is True:
            f_dEdx = open(dir_name + "Epair_dEdx{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["dEdx"] = [f_dEdx, [""]]
        if key == "dNdx" and kwargs[key] is True:
            f_dNdx = open(dir_name + "Epair_dNdx{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["dNdx"] = [f_dNdx, [""]]
        if key == "dNdx_rnd" and kwargs[key] is True:
            f_dNdx_rnd = open(dir_name + "Epair_dNdx_rnd{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["dNdx_rnd"] = [f_dNdx_rnd, [""]]
        if key == "stoch" and kwargs[key] is True:
            f_stoch = open(dir_name + "Epair_e{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["stoch"] = [f_stoch, [""]]

    # print(buf)

    for particle in particle_defs:
        for medium in mediums:
            for cut in cuts:
                for lpm in lpms:
                    for param in params:

                        if interpolate:
                            param_current = param(
                                particle,
                                medium,
                                cut,
                                multiplier,
                                lpm,
                                interpoldef)

                            xsection = pp.crosssection.EpairInterpolant(param_current, interpoldef)
                        else:
                            param_current = param(
                                particle,
                                medium,
                                cut,
                                multiplier,
                                lpm)

                            xsection = pp.crosssection.EpairIntegral(param_current)

                        for key in buf:
                            buf[key][1] = [""]

                            for energy in energies:
                                if key == "dEdx":
                                    result = [str(xsection.calculate_dEdx(energy))]
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
                                buf[key][1].append(str(cut.ecut))
                                buf[key][1].append(str(cut.vcut))
                                buf[key][1].append(str(multiplier))
                                buf[key][1].append(str(lpm))
                                buf[key][1].append(str(energy))
                                buf[key][1].append(param_current.name)
                                buf[key][1].extend(result)
                                buf[key][1].append("\n")

                            buf[key][0].write("\t".join(buf[key][1]))


def create_table_dNdx(dir_name, interpolate=False):

    if interpolate:
        params = epair_interpol
    else:
        params = epair

    with open(dir_name + "Epair_dNdx{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for param in params:
                        for lpm in lpms:

                            if interpolate:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    lpm,
                                    interpoldef)

                                xsection = pp.crosssection.EpairInterpolant(param_current, interpoldef)
                            else:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    lpm)

                                xsection = pp.crosssection.EpairIntegral(param_current)

                            buf = [""]

                            for energy in energies:
                                dEdx = xsection.calculate_dNdx(energy)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(lpm))
                                buf.append(str(energy))
                                buf.append(str(dEdx))
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_dNdx_rnd(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(0)

    if interpolate:
        params = epair_interpol
    else:
        params = epair

    with open(dir_name + "Epair_dNdx_rnd{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for param in params:
                        for lpm in lpms:

                            if interpolate:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    lpm,
                                    interpoldef)

                                xsection = pp.crosssection.EpairInterpolant(param_current, interpoldef)
                            else:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    lpm)

                                xsection = pp.crosssection.EpairIntegral(param_current)

                            buf = [""]

                            for energy in energies:
                                rnd = pp.RandomGenerator.get().random_double()
                                dNdx = xsection.calculate_dNdx_rnd(energy, rnd)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(lpm))
                                buf.append(str(energy))
                                buf.append(str(rnd))
                                buf.append(str(dNdx))
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_stochastic_loss(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(0)

    if interpolate:
        params = epair_interpol
    else:
        params = epair

    with open(dir_name + "Epair_e{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for param in params:
                        for lpm in lpms:

                            if interpolate:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    lpm,
                                    interpoldef)

                                xsection = pp.crosssection.EpairInterpolant(param_current, interpoldef)
                            else:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    lpm)

                                xsection = pp.crosssection.EpairIntegral(param_current)

                            buf = [""]

                            for energy in energies:
                                rnd1 = pp.RandomGenerator.get().random_double()
                                rnd2 = pp.RandomGenerator.get().random_double()
                                stochastic_loss = xsection.calculate_dNdx_rnd(energy, rnd1, rnd2)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(lpm))
                                buf.append(str(energy))
                                buf.append(str(rnd1))
                                buf.append(str(rnd2))
                                buf.append(str(stochastic_loss))
                                buf.append("\n")

                            file.write("\t".join(buf))


def main(dir_name):
    create_tables(dir_name, False, dEdx=True, dNdx=True, dNdx_rnd=True, stoch=True)
    create_tables(dir_name, True, dEdx=True, dNdx=True, dNdx_rnd=True, stoch=True)

if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
