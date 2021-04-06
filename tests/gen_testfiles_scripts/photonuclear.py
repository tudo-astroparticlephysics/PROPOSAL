import os
import proposal as pp
import numpy as np

particle_defs = [
    pp.particle.MuMinusDef(),
    pp.particle.TauMinusDef()#,
    # pp.particle.EMinusDef()
]

mediums = [
    pp.medium.Ice(),
    pp.medium.Hydrogen(),
    pp.medium.Uranium()
]

cuts = [
    pp.EnergyCutSettings(np.inf, 1),
    pp.EnergyCutSettings(500, 1),
    pp.EnergyCutSettings(np.inf, 0.05),
    pp.EnergyCutSettings(500, 0.05)
]

multiplier = 1.

hard_components = [0, 1]

shadows = [
    pp.parametrization.photonuclear.ShadowDuttaRenoSarcevicSeckel(),
    pp.parametrization.photonuclear.ShadowButkevichMikheyev()
]

names_shadow = [
    "ShadowDuttaRenoSarcevicSeckel",
    "ShadowButkevichMikheyev"
]

energies = np.logspace(4, 13, num=10)


def create_tables_real(dir_name, **kwargs):

    pp.RandomGenerator.get().set_seed(1234)

    buf = {}

    for key in kwargs:
        if key == "dEdx" and kwargs[key] is True:
            f_dEdx = open(dir_name + "Photo_Real_dEdx.txt", "w")
            buf["dEdx"] = [f_dEdx, [""]]
        if key == "dNdx" and kwargs[key] is True:
            f_dNdx = open(dir_name + "Photo_Real_dNdx.txt", "w")
            buf["dNdx"] = [f_dNdx, [""]]
        if key == "stoch" and kwargs[key] is True:
            f_dNdx_rnd = open(dir_name + "Photo_Real_e.txt", "w")
            buf["stoch"] = [f_dNdx_rnd, [""]]


    for particle in particle_defs:
        for medium in mediums:
            for cut in cuts:
                for hard  in hard_components:

                    parametrizations = [
                        pp.parametrization.photonuclear.Zeus(hard),
                        pp.parametrization.photonuclear.BezrukovBugaev(hard),
                        pp.parametrization.photonuclear.Rhode(hard),
                        pp.parametrization.photonuclear.Kokoulin(hard)
                    ]

                    for parametrization in parametrizations:
                        args = {
                            "parametrization": parametrization,
                            "particle_def": particle,
                            "target": medium,
                            "cuts": cut,
                            "interpolate": False
                        }

                        xsection = pp.crosssection.make_crosssection(**args)


                        for key in buf:
                            buf[key][1] = [""]

                            for energy in energies:
                                if key == "dEdx":
                                    result = [str(xsection.calculate_dEdx(energy) * medium.mass_density)]
                                if key == "dNdx":
                                    result = [str(xsection.calculate_dNdx(energy) * medium.mass_density)]
                                if key == "stoch":
                                    rnd1 = pp.RandomGenerator.get().random_double()
                                    rnd2 = pp.RandomGenerator.get().random_double()

                                    components = medium.components
                                    comp = components[int(rnd2*len(components))]
                                    dNdx_for_comp = xsection.calculate_dNdx(energy, comp.hash);

                                    if np.isfinite(cut.ecut) or cut.vcut < 1:
                                        result = xsection.calculate_stochastic_loss(
                                            comp.hash, energy, rnd1*dNdx_for_comp) * energy
                                    else:
                                        result = 0
                                    result = [str(rnd1), str(rnd2), str(result)]


                                buf[key][1].append(particle.name)
                                buf[key][1].append(medium.name)
                                buf[key][1].append(str(cut.ecut))
                                buf[key][1].append(str(cut.vcut))
                                buf[key][1].append(str(multiplier))
                                buf[key][1].append(str(energy))
                                buf[key][1].append(str(hard))
                                buf[key][1].append(xsection.param_name)
                                buf[key][1].extend(result)
                                buf[key][1].append("\n")

                            buf[key][0].write("\t".join(buf[key][1]))



def create_tables_Q2(dir_name, **kwargs):

    pp.RandomGenerator.get().set_seed(1234)

    buf = {}

    for key in kwargs:
        if key == "dEdx" and kwargs[key] is True:
            f_dEdx = open(dir_name + "Photo_Q2_dEdx.txt", "w")
            buf["dEdx"] = [f_dEdx, [""]]
        if key == "dNdx" and kwargs[key] is True:
            f_dNdx = open(dir_name + "Photo_Q2_dNdx.txt", "w")
            buf["dNdx"] = [f_dNdx, [""]]
        if key == "stoch" and kwargs[key] is True:
            f_dNdx_rnd = open(dir_name + "Photo_Q2_e.txt", "w")
            buf["stoch"] = [f_dNdx_rnd, [""]]


    for particle in particle_defs:
        for medium in mediums:
            for cut in cuts:
                for sidx, shadow in enumerate(shadows):

                    parametrizations = [
                        pp.parametrization.photonuclear.AbramowiczLevinLevyMaor91(shadow),
                        pp.parametrization.photonuclear.AbramowiczLevinLevyMaor97(shadow),
                        pp.parametrization.photonuclear.ButkevichMikheyev(shadow),
                        pp.parametrization.photonuclear.RenoSarcevicSu(shadow),
                        pp.parametrization.photonuclear.AbtFT(shadow),
                        pp.parametrization.photonuclear.BlockDurandHa(shadow)
                    ]
                    name_shadow = names_shadow[sidx]

                    for parametrization in parametrizations:
                        args = {
                            "parametrization": parametrization,
                            "particle_def": particle,
                            "target": medium,
                            "cuts": cut,
                            "interpolate": False
                        }

                        xsection = pp.crosssection.make_crosssection(**args)


                        for key in buf:
                            buf[key][1] = [""]

                            for energy in energies:
                                if key == "dEdx":
                                    result = [str(xsection.calculate_dEdx(energy) * medium.mass_density)]
                                if key == "dNdx":
                                    result = [str(xsection.calculate_dNdx(energy) * medium.mass_density)]
                                if key == "stoch":
                                    rnd1 = pp.RandomGenerator.get().random_double()
                                    rnd2 = pp.RandomGenerator.get().random_double()

                                    components = medium.components
                                    comp = components[int(rnd2*len(components))]
                                    dNdx_for_comp = xsection.calculate_dNdx(energy, comp.hash);

                                    if np.isfinite(cut.ecut) or cut.vcut < 1:
                                        result = xsection.calculate_stochastic_loss(
                                            comp.hash, energy, rnd1*dNdx_for_comp) * energy
                                    else:
                                        result = 0
                                    result = [str(rnd1), str(rnd2), str(result)]


                                buf[key][1].append(particle.name)
                                buf[key][1].append(medium.name)
                                buf[key][1].append(str(cut.ecut))
                                buf[key][1].append(str(cut.vcut))
                                buf[key][1].append(str(multiplier))
                                buf[key][1].append(str(energy))
                                buf[key][1].append(name_shadow)
                                buf[key][1].append(xsection.param_name)
                                buf[key][1].extend(result)
                                buf[key][1].append("\n")

                            buf[key][0].write("\t".join(buf[key][1]))



def main(dir_name):
    create_tables_real(dir_name, dEdx=True, dNdx=True, stoch=True)
    create_tables_Q2(dir_name, dEdx=True, dNdx=True, stoch=True)


if __name__ == "__main__":

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
