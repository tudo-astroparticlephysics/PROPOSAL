import pyPROPOSAL as pp
import numpy as np

photo_real = [
    pp.Parametrization.Photonuclear.Zeus,
    pp.Parametrization.Photonuclear.BezrukovBugaev,
    pp.Parametrization.Photonuclear.Rhode,
    pp.Parametrization.Photonuclear.Kokoulin
]

particle_defs = [
    pp.MuMinusDef.get(),
    pp.TauMinusDef.get()#,
    # pp.EMinusDef.get()
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

hard_components = [0, 1]

photo_q2 = [
    pp.Parametrization.Photonuclear.AbramowiczLevinLevyMaor91,
    pp.Parametrization.Photonuclear.AbramowiczLevinLevyMaor97,
    pp.Parametrization.Photonuclear.ButkevichMikhailov,
    pp.Parametrization.Photonuclear.RenoSarcevicSu
]

shadows = [
    pp.Parametrization.Photonuclear.ShadowDuttaRenoSarcevicSeckel(),
    pp.Parametrization.Photonuclear.ShadowButkevichMikhailov()
]

energies = np.logspace(4, 13, num=10)

interpoldef = pp.InterpolationDef()

def create_table_dEdx():

    with open("TestFiles/Photo_Real_dEdx.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for hard  in hard_components:
                        for energy in energies:
                            for parametrization in photo_real:

                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    hard)

                                Photo_Int = pp.CrossSection.PhotoIntegral(photo)

                                buf = [""]

                                dEdx = Photo_Int.calculate_dEdx(energy)

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

                                # print(buf)
                                f.write("\t".join(buf))
#

def create_table_dNdx():

    with open("TestFiles/Photo_Real_dNdx.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for hard  in hard_components:
                        for energy in energies:
                            for parametrization in photo_real:

                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    hard)

                                Photo_Int = pp.CrossSection.PhotoIntegral(photo)

                                buf = [""]

                                dNdx = Photo_Int.calculate_dNdx(energy)

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

                                # print(buf)
                                f.write("\t".join(buf))
#

def create_table_dNdx_rnd():

    pp.RandomGenerator.get().set_seed(1234)

    with open("TestFiles/Photo_Real_dNdx_rnd.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for hard  in hard_components:
                        for energy in energies:
                            rnd = pp.RandomGenerator.get().random_double()
                            for parametrization in photo_real:

                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    hard)

                                Photo_Int = pp.CrossSection.PhotoIntegral(photo)

                                buf = [""]

                                dNdx = Photo_Int.calculate_dNdx_rnd(energy, rnd)

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

                                # print(buf)
                                f.write("\t".join(buf))
#

def create_table_stochastic_loss():

    pp.RandomGenerator.get().set_seed(1234)

    with open("TestFiles/Photo_Real_e.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for hard  in hard_components:
                        for energy in energies:
                            rnd1 = pp.RandomGenerator.get().random_double()
                            rnd2 = pp.RandomGenerator.get().random_double()
                            for parametrization in photo_real:

                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    hard)

                                Photo_Int = pp.CrossSection.PhotoIntegral(photo)

                                buf = [""]

                                stochastic_loss = Photo_Int.calculate_stochastic_loss(energy, rnd1, rnd2)

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

                                # print(buf)
                                f.write("\t".join(buf))
#

def create_table_dEdx_Q2():

    with open("TestFiles/Photo_Q2_dEdx.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for shadow in shadows:
                        for energy in energies:
                            for parametrization in photo_q2:

                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow)

                                Photo_Int = pp.CrossSection.PhotoIntegral(photo)

                                buf = [""]

                                dEdx = Photo_Int.calculate_dEdx(energy)

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

                                # print(buf)
                                f.write("\t".join(buf))
#

def create_table_dNdx_Q2():

    with open("TestFiles/Photo_Q2_dNdx.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for shadow in shadows:
                        for energy in energies:
                            for parametrization in photo_q2:

                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow)

                                Photo_Int = pp.CrossSection.PhotoIntegral(photo)

                                buf = [""]

                                dNdx = Photo_Int.calculate_dNdx(energy)

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

                                # print(buf)
                                f.write("\t".join(buf))
#

def create_table_dNdx_rnd_Q2():

    with open("TestFiles/Photo_Q2_dNdx_rnd.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for shadow in shadows:
                        for energy in energies:
                            rnd = pp.RandomGenerator.get().random_double()
                            for parametrization in photo_q2:

                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow)

                                Photo_Int = pp.CrossSection.PhotoIntegral(photo)

                                buf = [""]

                                dNdx = Photo_Int.calculate_dNdx_rnd(energy, rnd)

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

                                # print(buf)
                                f.write("\t".join(buf))
#

def create_table_stochastic_loss_Q2():

    with open("TestFiles/Photo_Q2_e.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:
                    for shadow in shadows:
                        for energy in energies:
                            rnd1 = pp.RandomGenerator.get().random_double()
                            rnd2 = pp.RandomGenerator.get().random_double()
                            for parametrization in photo_q2:

                                photo = parametrization(
                                    particle,
                                    medium,
                                    cut,
                                    multiplier,
                                    shadow)

                                Photo_Int = pp.CrossSection.PhotoIntegral(photo)

                                buf = [""]

                                stochastic_loss = Photo_Int.calculate_stochastic_loss(energy, rnd1, rnd2)

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

                                # print(buf)
                                f.write("\t".join(buf))
#



if __name__ == "__main__":
    # create_table_dEdx()
    # create_table_dNdx()
    # create_table_dNdx_rnd()
    # create_table_stochastic_loss()
    # create_table_dEdx_Q2()
    create_table_dNdx_Q2()
    create_table_dNdx_rnd_Q2()
    create_table_stochastic_loss_Q2()