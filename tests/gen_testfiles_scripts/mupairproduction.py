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

mupair = [
    pp.parametrization.mupairproduction.KelnerKokoulinPetrukhin,
]

mupair_interpol = [
    pp.parametrization.mupairproduction.KelnerKokoulinPetrukhinInterpolant,
]


energies = np.logspace(4, 13, num=10)

vs = np.linspace(0.022, 0.98, 10) # choose v-range so that v_min < v < v_max is always fulflilled for all E in energies

interpoldef = pp.InterpolationDef()


def create_tables(dir_name, interpolate=False, **kwargs):

    pp.RandomGenerator.get().set_seed(1234)

    if interpolate:
        params = mupair_interpol
    else:
        params = mupair

    buf = {}

    for key in kwargs:
        if key == "dEdx" and kwargs[key] is True:
            f_dEdx = open(dir_name + "Mupair_dEdx{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["dEdx"] = [f_dEdx, [""]]
        if key == "dNdx" and kwargs[key] is True:
            f_dNdx = open(dir_name + "Mupair_dNdx{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["dNdx"] = [f_dNdx, [""]]
        if key == "dNdx_rnd" and kwargs[key] is True:
            f_dNdx_rnd = open(dir_name + "Mupair_dNdx_rnd{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["dNdx_rnd"] = [f_dNdx_rnd, [""]]
        if key == "stoch" and kwargs[key] is True:
            f_stoch = open(dir_name + "Mupair_e{}.txt".format("_interpol" if interpolate else ""), "a")
            buf["stoch"] = [f_stoch, [""]]

    # print(buf)

    for particle in particle_defs:
        for medium in mediums:
            for cut in cuts:
                    for param in params:
                        if interpolate:
                            param_current = param(
                                particle,
                                medium,
                                cut,
                                multiplier,
                                interpoldef)

                            xsection = pp.crosssection.MupairInterpolant(param_current, interpoldef)
                        else:
                            param_current = param(
                                particle,
                                medium,
                                cut,
                                multiplier
                                )

                            xsection = pp.crosssection.MupairIntegral(param_current)

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
                                buf[key][1].append(str(energy))
                                buf[key][1].append(param_current.name)
                                buf[key][1].extend(result)
                                buf[key][1].append("\n")

                            buf[key][0].write("\t".join(buf[key][1]))


def create_table_dNdx(dir_name, interpolate=False):

    if interpolate:
        params = mupair_interpol
    else:
        params = mupair

    with open(dir_name + "Mupair_dNdx{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for param in params:
                            if interpolate:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    interpoldef)

                                xsection = pp.crosssection.MupairInterpolant(param_current, interpoldef)
                            else:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier)

                                xsection = pp.crosssection.MupairIntegral(param_current)

                            buf = [""]

                            for energy in energies:
                                dEdx = xsection.calculate_dNdx(energy)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(dEdx))
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_dNdx_rnd(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(0)

    if interpolate:
        params = mupair_interpol
    else:
        params = mupair

    with open(dir_name + "Mupair_dNdx_rnd{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for param in params:
                            if interpolate:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    interpoldef)

                                xsection = pp.crosssection.MupairInterpolant(param_current, interpoldef)
                            else:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier)

                                xsection = pp.crosssection.MupairIntegral(param_current)

                            buf = [""]

                            for energy in energies:
                                rnd = pp.RandomGenerator.get().random_double()
                                dNdx = xsection.calculate_dNdx_rnd(energy, rnd)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(rnd))
                                buf.append(str(dNdx))
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_stochastic_loss(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(0)

    if interpolate:
        params = mupair_interpol
    else:
        params = mupair

    with open(dir_name + "Mupair_e{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for param in params:
                            if interpolate:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    interpoldef)

                                xsection = pp.crosssection.MupairInterpolant(param_current, interpoldef)
                            else:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier)

                                xsection = pp.crosssection.MupairIntegral(param_current)

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
                                buf.append(str(energy))
                                buf.append(str(rnd1))
                                buf.append(str(rnd2))
                                buf.append(str(stochastic_loss))
                                buf.append("\n")

                            file.write("\t".join(buf))

def create_table_rho(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(0)

    if interpolate:
        params = mupair_interpol
    else:
        params = mupair

    with open(dir_name + "Mupair_rho{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for v in vs:
                        for param in params:
                            if interpolate:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    interpoldef)

                                xsection = pp.crosssection.MupairInterpolant(param_current, interpoldef)
                            else:
                                param_current = param(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier)

                                xsection = pp.crosssection.MupairIntegral(param_current)

                            buf = [""]

                            for energy in energies:
                                rnd1 = pp.RandomGenerator.get().random_double()
                                rnd2 = pp.RandomGenerator.get().random_double()
                                particles = xsection.calculate_produced_particles(energy, v*energy, rnd1, rnd2)
                                E1 = particles[0].energy
                                E2 = particles[1].energy

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(v))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(param_current.name))
                                buf.append(str(rnd1))
                                buf.append(str(rnd2))
                                buf.append(str(E1))
                                buf.append(str(E2))
                                buf.append("\n")

                            file.write("\t".join(buf))


def main(dir_name):
    create_tables(dir_name, False, dEdx=True, dNdx=True, dNdx_rnd=True, stoch=True)
    create_tables(dir_name, True, dEdx=True, dNdx=True, dNdx_rnd=True, stoch=True)
    create_table_rho(dir_name, False)


if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
