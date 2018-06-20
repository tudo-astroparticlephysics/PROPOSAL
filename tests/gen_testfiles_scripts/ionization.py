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

multiplier = 1.

energies = np.logspace(4, 13, num=10)

interpoldef = pp.InterpolationDef()


def create_table_dEdx(dir_name):

    with open(dir_name + "Ioniz_dEdx.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz = pp.Parametrization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )

                    Ioniz_Int = pp.CrossSection.IonizIntegral(Ioniz)

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

                    # print(buf)
                    f.write("\t".join(buf))


def create_table_dNdx(dir_name):

    with open(dir_name + "Ioniz_dNdx.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz = pp.Parametrization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )

                    Ioniz_Int = pp.CrossSection.IonizIntegral(Ioniz)

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

                    # print(buf)
                    f.write("\t".join(buf))


def create_table_dNdx_rnd(dirName):

    pp.RandomGenerator.get().set_seed(0)

    with open(dir_name + "Ioniz_dNdx_rnd.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz = pp.Parametrization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )

                    Ioniz_Int = pp.CrossSection.IonizIntegral(Ioniz)

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

                    # print(buf)
                    f.write("\t".join(buf))


def create_table_stochastic_loss(dir_name):

    pp.RandomGenerator.get().set_seed(5)

    with open(dir_name + "Ioniz_e.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz = pp.Parametrization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )

                    Ioniz_Int = pp.CrossSection.IonizIntegral(Ioniz)

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

                    # print(buf)
                    f.write("\t".join(buf))


def create_table_dEdx_interpol(dir_nameName):

    with open(dir_name + "Ioniz_dEdx_interpol.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz = pp.Parametrization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )

                    Ioniz_Interpol = pp.CrossSection.IonizInterpolant(Ioniz, interpoldef)

                    buf = [""]

                    for energy in energies:
                        dEdx = Ioniz_Interpol.calculate_dEdx(energy)

                        buf.append(particle.name)
                        buf.append(medium.name)
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(multiplier))
                        buf.append(str(energy))
                        buf.append(str(dEdx))
                        buf.append("\n")

                    # print(buf)
                    f.write("\t".join(buf))


def create_table_dNdx_interpol(dir_name):

    with open(dir_name + "Ioniz_dNdx_interpol.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz = pp.Parametrization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )

                    Ioniz_Interpol = pp.CrossSection.IonizInterpolant(Ioniz, interpoldef)

                    buf = [""]

                    for energy in energies:
                        dNdx = Ioniz_Interpol.calculate_dNdx(energy)

                        buf.append(particle.name)
                        buf.append(medium.name)
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(multiplier))
                        buf.append(str(energy))
                        buf.append(str(dNdx))
                        buf.append("\n")

                    # print(buf)
                    f.write("\t".join(buf))


def create_table_dNdx_rnd_interpol(dir_name):

    pp.RandomGenerator.get().set_seed(5)

    with open(dir_name + "Ioniz_dNdx_rnd_interpol.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz = pp.Parametrization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )

                    Ioniz_Interpol = pp.CrossSection.IonizInterpolant(Ioniz, interpoldef)

                    buf = [""]

                    for energy in energies:
                        rnd = pp.RandomGenerator.get().random_double()
                        dNdx_rnd = Ioniz_Interpol.calculate_dNdx_rnd(energy, rnd)

                        buf.append(particle.name)
                        buf.append(medium.name)
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(multiplier))
                        buf.append(str(energy))
                        buf.append(str(rnd))
                        buf.append(str(dNdx_rnd))
                        buf.append("\n")

                    # print(buf)
                    f.write("\t".join(buf))


def create_table_stochastic_loss_interpol(dir_name):

    pp.RandomGenerator.get().set_seed(5)

    with open(dir_name + "Ioniz_e_interpol.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    Ioniz = pp.Parametrization.Ionization(
                        particle,
                        medium,
                        cut,
                        multiplier
                    )

                    Ioniz_Interpol = pp.CrossSection.IonizInterpolant(Ioniz, interpoldef)

                    buf = [""]

                    for energy in energies:
                        rnd1 = pp.RandomGenerator.get().random_double()
                        rnd2 = pp.RandomGenerator.get().random_double()
                        stochastic_loss = Ioniz_Interpol.calculate_stochastic_loss(energy, rnd1, rnd2)

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

                    # print(buf)
                    f.write("\t".join(buf))


def main(dir_name):
    create_table_dEdx(dir_name)
    create_table_dNdx(dir_name)
    create_table_dNdx_rnd(dir_name)
    create_table_stochastic_loss(dir_name)
    create_table_dEdx_interpol(dir_name)
    create_table_dNdx_interpol(dir_name)
    create_table_dNdx_rnd_interpol(dir_name)
    create_table_stochastic_loss_interpol(dir_name)


if __name__ == "__main__":
    main()
