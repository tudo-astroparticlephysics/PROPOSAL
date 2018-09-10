import pyPROPOSAL as pp
import numpy as np

photo_real = [
    pp.parametrization.photonuclear.Zeus,
    pp.parametrization.photonuclear.BezrukovBugaev,
    pp.parametrization.photonuclear.Rhode,
    pp.parametrization.photonuclear.Kokoulin
]

particle_defs = [
    pp.particle.MuMinusDef.get(),
    pp.particle.TauMinusDef.get()#,
    # pp.particle.EMinusDef.get()
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

hard_components = [0, 1]

photo_q2 = [
    pp.parametrization.photonuclear.AbramowiczLevinLevyMaor91,
    pp.parametrization.photonuclear.AbramowiczLevinLevyMaor97,
    pp.parametrization.photonuclear.ButkevichMikhailov,
    pp.parametrization.photonuclear.RenoSarcevicSu
]

photo_q2_interpol = [
    pp.parametrization.photonuclear.AbramowiczLevinLevyMaor91Interpolant,
    pp.parametrization.photonuclear.AbramowiczLevinLevyMaor97Interpolant,
    pp.parametrization.photonuclear.ButkevichMikhailovInterpolant,
    pp.parametrization.photonuclear.RenoSarcevicSuInterpolant
]

shadows = [
    pp.parametrization.photonuclear.ShadowDuttaRenoSarcevicSeckel(),
    pp.parametrization.photonuclear.ShadowButkevichMikhailov()
]

energies = np.logspace(4, 13, num=10)

interpoldef = pp.InterpolationDef()


def create_table_dEdx(dir_name, interpolate=False):

    with open(dir_name + "Photo_Real_dEdx{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for hard  in hard_components:
                        for parametrization in photo_real:

                            photo = parametrization(
                                particle,
                                medium,
                                cut,
                                multiplier,
                                hard)

                            if interpolate:
                                xsection = pp.crosssection.PhotoInterpolant(photo, interpoldef)
                            else:
                                xsection = pp.crosssection.PhotoIntegral(photo)

                            buf = [""]

                            for energy in energies:
                                dEdx = xsection.calculate_dEdx(energy)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(dEdx))
                                buf.append(photo.name)
                                buf.append(str(hard))
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_dNdx(dir_name, interpolate=False):

    with open(dir_name + "Photo_Real_dNdx{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for hard  in hard_components:
                        for parametrization in photo_real:

                            photo = parametrization(
                                particle,
                                medium,
                                cut,
                                multiplier,
                                hard)

                            if interpolate:
                                xsection = pp.crosssection.PhotoInterpolant(photo, interpoldef)
                            else:
                                xsection = pp.crosssection.PhotoIntegral(photo)

                            buf = [""]

                            for energy in energies:
                                dNdx = xsection.calculate_dNdx(energy)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(dNdx))
                                buf.append(photo.name)
                                buf.append(str(hard))
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_dNdx_rnd(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(1234)

    with open(dir_name + "Photo_Real_dNdx_rnd{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for hard  in hard_components:
                        rnd = pp.RandomGenerator.get().random_double()
                        for parametrization in photo_real:

                            photo = parametrization(
                                particle,
                                medium,
                                cut,
                                multiplier,
                                hard)

                            if interpolate:
                                xsection = pp.crosssection.PhotoInterpolant(photo, interpoldef)
                            else:
                                xsection = pp.crosssection.PhotoIntegral(photo)

                            buf = [""]

                            for energy in energies:
                                dNdx = xsection.calculate_dNdx_rnd(energy, rnd)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(rnd))
                                buf.append(str(dNdx))
                                buf.append(photo.name)
                                buf.append(str(hard))
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_stochastic_loss(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(1234)

    with open(dir_name + "Photo_Real_e{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for hard  in hard_components:
                        for parametrization in photo_real:

                            photo = parametrization(
                                particle,
                                medium,
                                cut,
                                multiplier,
                                hard)

                            if interpolate:
                                xsection = pp.crosssection.PhotoInterpolant(photo, interpoldef)
                            else:
                                xsection = pp.crosssection.PhotoIntegral(photo)

                            buf = [""]

                            for energy in energies:
                                rnd1 = pp.RandomGenerator.get().random_double()
                                rnd2 = pp.RandomGenerator.get().random_double()
                                stochastic_loss = xsection.calculate_stochastic_loss(energy, rnd1, rnd2)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(rnd1))
                                buf.append(str(rnd2))
                                buf.append(str(stochastic_loss))
                                buf.append(photo.name)
                                buf.append(str(hard))
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_dEdx_Q2(dir_name, interpolate=False):

    if interpolate:
        q2 = photo_q2_interpol
    else:
        q2 = photo_q2

    with open(dir_name + "Photo_Q2_dEdx{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for shadow in shadows:
                        for parametrization in q2:

                            if interpolate:
                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow,
                                    interpoldef)

                                xsection = pp.crosssection.PhotoInterpolant(photo, interpoldef)
                            else:
                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow)

                                xsection = pp.crosssection.PhotoIntegral(photo)

                            buf = [""]

                            for energy in energies:
                                dEdx = xsection.calculate_dEdx(energy)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(dEdx))
                                buf.append(photo.name)
                                buf.append(shadow.name)
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_dNdx_Q2(dir_name, interpolate=False):

    if interpolate:
        q2 = photo_q2_interpol
    else:
        q2 = photo_q2

    with open(dir_name + "Photo_Q2_dNdx{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for shadow in shadows:
                        for parametrization in q2:

                            if interpolate:
                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow,
                                    interpoldef)

                                xsection = pp.crosssection.PhotoInterpolant(photo, interpoldef)
                            else:
                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow)

                                xsection = pp.crosssection.PhotoIntegral(photo)

                            buf = [""]

                            for energy in energies:
                                dNdx = xsection.calculate_dNdx(energy)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(dNdx))
                                buf.append(photo.name)
                                buf.append(shadow.name)
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_dNdx_rnd_Q2(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(1234)

    if interpolate:
        q2 = photo_q2_interpol
    else:
        q2 = photo_q2

    with open(dir_name + "Photo_Q2_dNdx_rnd{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for shadow in shadows:
                        rnd = pp.RandomGenerator.get().random_double()
                        for parametrization in q2:

                            if interpolate:
                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow,
                                    interpoldef)

                                xsection = pp.crosssection.PhotoInterpolant(photo, interpoldef)
                            else:
                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow)

                                xsection = pp.crosssection.PhotoIntegral(photo)

                            buf = [""]

                            for energy in energies:
                                dNdx = xsection.calculate_dNdx_rnd(energy, rnd)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(rnd))
                                buf.append(str(dNdx))
                                buf.append(photo.name)
                                buf.append(shadow.name)
                                buf.append("\n")

                            file.write("\t".join(buf))


def create_table_stochastic_loss_Q2(dir_name, interpolate=False):

    pp.RandomGenerator.get().set_seed(1234)

    if interpolate:
        q2 = photo_q2_interpol
    else:
        q2 = photo_q2

    with open(dir_name + "Photo_Q2_e{}.txt".format("_interpol" if interpolate else ""), "a") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for shadow in shadows:
                        for parametrization in q2:

                            if interpolate:
                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow,
                                    interpoldef)

                                xsection = pp.crosssection.PhotoInterpolant(photo, interpoldef)
                            else:
                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow)

                                xsection = pp.crosssection.PhotoIntegral(photo)

                            buf = [""]

                            for energy in energies:
                                rnd1 = pp.RandomGenerator.get().random_double()
                                rnd2 = pp.RandomGenerator.get().random_double()
                                stochastic_loss = xsection.calculate_stochastic_loss(energy, rnd1, rnd2)

                                buf.append(particle.name)
                                buf.append(medium.name)
                                buf.append(str(cut.ecut))
                                buf.append(str(cut.vcut))
                                buf.append(str(multiplier))
                                buf.append(str(energy))
                                buf.append(str(rnd1))
                                buf.append(str(rnd2))
                                buf.append(str(stochastic_loss))
                                buf.append(photo.name)
                                buf.append(shadow.name)
                                buf.append("\n")

                            file.write("\t".join(buf))


def main(dir_name):
    # Integrate
    create_table_dEdx(dir_name)
    create_table_dNdx(dir_name)
    create_table_dNdx_rnd(dir_name)
    create_table_stochastic_loss(dir_name)
    create_table_dEdx_Q2(dir_name)
    create_table_dNdx_Q2(dir_name)
    create_table_dNdx_rnd_Q2(dir_name)
    create_table_stochastic_loss_Q2(dir_name)

    # Interpolate
    create_table_dEdx(dir_name, True)
    create_table_dNdx(dir_name, True)
    create_table_dNdx_rnd(dir_name, True)
    create_table_stochastic_loss(dir_name, True)
    create_table_dEdx_Q2(dir_name, True)
    create_table_dNdx_Q2(dir_name, True)
    create_table_dNdx_rnd_Q2(dir_name, True)
    create_table_stochastic_loss_Q2(dir_name, True)


if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
