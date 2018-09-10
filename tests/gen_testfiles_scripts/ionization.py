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

energies = np.logspace(4, 13, num=10)

interpoldef = pp.InterpolationDef()


def create_table_dEdx(dir_name, interpolate=False):

    with open(dir_name + "Ioniz_dEdx{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz_param = pp.parametrization.ionization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )
                    if interpolate:
                        Ioniz_Int = pp.crosssection.IonizInterpolant(Ioniz_param, interpoldef)
                    else:
                        Ioniz_Int = pp.crosssection.IonizIntegral(Ioniz_param)

                    buf = [""]

                    for energy in energies:
                        dEdx = Ioniz_Int.calculate_dEdx(energy)

                        buf.append(particle.name)
                        buf.append(medium.name)
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(multiplier))
                        buf.append(str(energy))
                        buf.append(str(dEdx))
                        buf.append("\n")

                    file.write("\t".join(buf))


def create_table_dNdx(dir_name, interpolate=False):

    with open(dir_name + "Ioniz_dNdx{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz_param = pp.parametrization.ionization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )
                    if interpolate:
                        Ioniz_Int = pp.crosssection.IonizInterpolant(Ioniz_param, interpoldef)
                    else:
                        Ioniz_Int = pp.crosssection.IonizIntegral(Ioniz_param)

                    buf = [""]

                    for energy in energies:
                        dNdx = Ioniz_Int.calculate_dNdx(energy)

                        buf.append(particle.name)
                        buf.append(medium.name)
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(multiplier))
                        buf.append(str(energy))
                        buf.append(str(dNdx))
                        buf.append("\n")

                    file.write("\t".join(buf))


def create_table_dNdx_rnd(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(1234)

    with open(dir_name + "Ioniz_dNdx_rnd{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz_param = pp.parametrization.ionization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )
                    if interpolate:
                        Ioniz_Int = pp.crosssection.IonizInterpolant(Ioniz_param, interpoldef)
                    else:
                        Ioniz_Int = pp.crosssection.IonizIntegral(Ioniz_param)

                    buf = [""]

                    for energy in energies:
                        rnd = pp.RandomGenerator.get().random_double()
                        dNdx_rnd = Ioniz_Int.calculate_dNdx_rnd(energy, rnd)

                        buf.append(particle.name)
                        buf.append(medium.name)
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(multiplier))
                        buf.append(str(energy))
                        buf.append(str(rnd))
                        buf.append(str(dNdx_rnd))
                        buf.append("\n")

                    file.write("\t".join(buf))


def create_table_stochastic_loss(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(1234)

    with open(dir_name + "Ioniz_e{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz_param = pp.parametrization.ionization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )
                    if interpolate:
                        Ioniz_Int = pp.crosssection.IonizInterpolant(Ioniz_param, interpoldef)
                    else:
                        Ioniz_Int = pp.crosssection.IonizIntegral(Ioniz_param)

                    buf = [""]

                    for energy in energies:
                        rnd1 = pp.RandomGenerator.get().random_double()
                        rnd2 = pp.RandomGenerator.get().random_double()
                        stochastic_loss = Ioniz_Int.calculate_stochastic_loss(energy, rnd1, rnd2)

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


def main(dir_name):
    create_table_dEdx(dir_name, interpolate=False)
    create_table_dNdx(dir_name, interpolate=False)
    create_table_dNdx_rnd(dir_name, interpolate=False)
    create_table_stochastic_loss(dir_name, interpolate=False)
    create_table_dEdx(dir_name, interpolate=True)
    create_table_dNdx(dir_name, interpolate=True)
    create_table_dNdx_rnd(dir_name, interpolate=True)
    create_table_stochastic_loss(dir_name, interpolate=True)


if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
