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

lpms = [0, 1]

energies = np.logspace(4, 13, num=10)

interpoldef = pp.InterpolationDef()


def create_table_dEdx():

    with open("Epair_dEdx.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for lpm  in lpms:

                        Epair = pp.Parametrization.EpairProduction.EpairProductionRhoIntegral(
                            particle,
                            medium,
                            cut,
                            multiplier,
                            lpm
                        )

                        Epair_Int = pp.CrossSection.EpairIntegral(Epair)

                        buf = [""]

                        for energy in energies:
                            dEdx = Epair_Int.calculate_dEdx(energy)

                            buf.append(particle.name)
                            buf.append(medium.name)
                            buf.append(str(cut.ecut))
                            buf.append(str(cut.vcut))
                            buf.append(str(multiplier))
                            buf.append(str(lpm))
                            buf.append(str(energy))
                            buf.append(str(dEdx))
                            buf.append("\n")

                        # print(buf)
                        f.write("\t".join(buf))


def create_table_dNdx():

    with open("Epair_dNdx.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for lpm  in lpms:

                        Epair = pp.Parametrization.EpairProduction.EpairProductionRhoIntegral(
                            particle,
                            medium,
                            cut,
                            multiplier,
                            lpm
                        )

                        Epair_Int = pp.CrossSection.EpairIntegral(Epair)

                        buf = [""]

                        for energy in energies:
                            dNdx = Epair_Int.calculate_dNdx(energy)

                            buf.append(particle.name)
                            buf.append(medium.name)
                            buf.append(str(cut.ecut))
                            buf.append(str(cut.vcut))
                            buf.append(str(multiplier))
                            buf.append(str(lpm))
                            buf.append(str(energy))
                            buf.append(str(dNdx))
                            buf.append("\n")

                        # print(buf)
                        f.write("\t".join(buf))


def create_table_dNdx_rnd():

    pp.RandomGenerator.get().set_seed(0)

    with open("Epair_dNdx_rnd.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for lpm  in lpms:

                        Epair = pp.Parametrization.EpairProduction.EpairProductionRhoIntegral(
                            particle,
                            medium,
                            cut,
                            multiplier,
                            lpm
                        )

                        Epair_Int = pp.CrossSection.EpairIntegral(Epair)

                        buf = [""]

                        for energy in energies:
                            rnd = pp.RandomGenerator.get().random_double()
                            dNdx_rnd = Epair_Int.calculate_dNdx_rnd(energy, rnd)

                            buf.append(particle.name)
                            buf.append(medium.name)
                            buf.append(str(cut.ecut))
                            buf.append(str(cut.vcut))
                            buf.append(str(multiplier))
                            buf.append(str(lpm))
                            buf.append(str(energy))
                            buf.append(str(rnd))
                            buf.append(str(dNdx_rnd))
                            buf.append("\n")

                        # print(buf)
                        f.write("\t".join(buf))


def create_table_stochastic_loss():

    pp.RandomGenerator.get().set_seed(5)

    with open("Epair_e.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for lpm  in lpms:

                        Epair = pp.Parametrization.EpairProduction.EpairProductionRhoIntegral(
                            particle,
                            medium,
                            cut,
                            multiplier,
                            lpm
                        )

                        Epair_Int = pp.CrossSection.EpairIntegral(Epair)

                        buf = [""]

                        for energy in energies:
                            rnd1 = pp.RandomGenerator.get().random_double()
                            rnd2 = pp.RandomGenerator.get().random_double()
                            stochastic_loss = Epair_Int.calculate_stochastic_loss(energy, rnd1, rnd2)

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

                        # print(buf)
                        f.write("\t".join(buf))


def create_table_dEdx_interpol():

    with open("Epair_dEdx_interpol.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for lpm  in lpms:

                        Epair = pp.Parametrization.EpairProduction.EpairProductionRhoInterpolant(
                            particle,
                            medium,
                            cut,
                            multiplier,
                            lpm,
                            interpoldef
                        )

                        Epair_Interpol = pp.CrossSection.EpairInterpolant(Epair, interpoldef)

                        buf = [""]

                        for energy in energies:
                            dEdx = Epair_Interpol.calculate_dEdx(energy)

                            buf.append(particle.name)
                            buf.append(medium.name)
                            buf.append(str(cut.ecut))
                            buf.append(str(cut.vcut))
                            buf.append(str(multiplier))
                            buf.append(str(lpm))
                            buf.append(str(energy))
                            buf.append(str(dEdx))
                            buf.append("\n")

                        # print(buf)
                        f.write("\t".join(buf))


def create_table_dNdx_interpol():

    with open("Epair_dNdx_interpol.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for lpm  in lpms:

                        Epair = pp.Parametrization.EpairProduction.EpairProductionRhoInterpolant(
                            particle,
                            medium,
                            cut,
                            multiplier,
                            lpm,
                            interpoldef
                        )

                        Epair_Interpol = pp.CrossSection.EpairInterpolant(Epair, interpoldef)

                        buf = [""]

                        for energy in energies:
                            dNdx = Epair_Interpol.calculate_dNdx(energy)

                            buf.append(particle.name)
                            buf.append(medium.name)
                            buf.append(str(cut.ecut))
                            buf.append(str(cut.vcut))
                            buf.append(str(multiplier))
                            buf.append(str(lpm))
                            buf.append(str(energy))
                            buf.append(str(dNdx))
                            buf.append("\n")

                        # print(buf)
                        f.write("\t".join(buf))


def create_table_dNdx_rnd_interpol():

    pp.RandomGenerator.get().set_seed(5)

    with open("Epair_dNdx_rnd_interpol.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for lpm  in lpms:

                        Epair = pp.Parametrization.EpairProduction.EpairProductionRhoInterpolant(
                            particle,
                            medium,
                            cut,
                            multiplier,
                            lpm,
                            interpoldef
                        )

                        Epair_Interpol = pp.CrossSection.EpairInterpolant(Epair, interpoldef)

                        buf = [""]

                        for energy in energies:
                            rnd = pp.RandomGenerator.get().random_double()
                            dNdx_rnd = Epair_Interpol.calculate_dNdx_rnd(energy, rnd)

                            buf.append(particle.name)
                            buf.append(medium.name)
                            buf.append(str(cut.ecut))
                            buf.append(str(cut.vcut))
                            buf.append(str(multiplier))
                            buf.append(str(lpm))
                            buf.append(str(energy))
                            buf.append(str(rnd))
                            buf.append(str(dNdx_rnd))
                            buf.append("\n")

                        # print(buf)
                        f.write("\t".join(buf))


def create_table_stochastic_loss_interpol():

    pp.RandomGenerator.get().set_seed(5)

    with open("Epair_e_interpol.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for lpm  in lpms:

                        Epair = pp.Parametrization.EpairProduction.EpairProductionRhoInterpolant(
                            particle,
                            medium,
                            cut,
                            multiplier,
                            lpm,
                            interpoldef
                        )

                        Epair_Interpol = pp.CrossSection.EpairInterpolant(Epair, interpoldef)

                        buf = [""]

                        for energy in energies:
                            rnd1 = pp.RandomGenerator.get().random_double()
                            rnd2 = pp.RandomGenerator.get().random_double()
                            stochastic_loss = Epair_Interpol.calculate_stochastic_loss(energy, rnd1, rnd2)

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

                        # print(buf)
                        f.write("\t".join(buf))


def main():
    create_table_dEdx()
    create_table_dNdx()
    create_table_dNdx_rnd()
    create_table_stochastic_loss()
    create_table_dEdx_interpol()
    create_table_dNdx_interpol()
    create_table_dNdx_rnd_interpol()
    create_table_stochastic_loss_interpol()


if __name__ == "__main__":
    main()
