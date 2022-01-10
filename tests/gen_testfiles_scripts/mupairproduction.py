import os
import proposal as pp
import numpy as np

particle_defs = [
    pp.particle.MuMinusDef(),
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

params = [
    pp.parametrization.mupairproduction.KelnerKokoulinPetrukhin(),
]

energies = np.logspace(4, 13, num=10)

vs = np.linspace(0.022, 0.98, 10) # choose v-range so that v_min < v < v_max is always fulflilled for all E in energies


def create_tables(dir_name, **kwargs):

    pp.RandomGenerator.get().set_seed(1234)

    buf = {}

    for key in kwargs:
        if key == "dEdx" and kwargs[key] is True:
            f_dEdx = open(dir_name + "Mupair_dEdx.txt", "w")
            buf["dEdx"] = [f_dEdx, [""]]
        if key == "dNdx" and kwargs[key] is True:
            f_dNdx = open(dir_name + "Mupair_dNdx.txt", "w")
            buf["dNdx"] = [f_dNdx, [""]]
        if key == "stoch" and kwargs[key] is True:
            f_stoch = open(dir_name + "Mupair_e.txt", "w")
            buf["stoch"] = [f_stoch, [""]]


    for particle in particle_defs:
        for medium in mediums:
            for cut in cuts:
                for param in params:
                    args = {
                        "parametrization": param,
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
                                result = [str(xsection.calculate_dEdx(energy))]
                            if key == "dNdx":
                                result = [str(xsection.calculate_dNdx(energy))]
                            if key == "stoch":
                                rnd1 = pp.RandomGenerator.get().random_double()
                                rnd2 = pp.RandomGenerator.get().random_double()

                                components = medium.components
                                comp = components[int(rnd2*len(components))]
                                dNdx_for_comp = xsection.calculate_dNdx(energy, comp.hash);

                                if np.isfinite(cut.ecut) or cut.vcut < 1:
                                    result = xsection.calculate_stochastic_loss(
                                        comp.hash, energy, rnd1*dNdx_for_comp)
                                else:
                                    result = 0
                                result = [str(rnd1), str(rnd2), str(result)]


                            buf[key][1].append(particle.name)
                            buf[key][1].append(medium.name)
                            buf[key][1].append(str(cut.ecut))
                            buf[key][1].append(str(cut.vcut))
                            buf[key][1].append(str(multiplier))
                            buf[key][1].append(str(energy))
                            buf[key][1].append(xsection.param_name)
                            buf[key][1].extend(result)
                            buf[key][1].append("\n")

                        buf[key][0].write("\t".join(buf[key][1]))


def create_table_rho(dir_name):

    pp.RandomGenerator.get().set_seed(0)

    with open(dir_name + "Mupair_rho.txt", "w") as file:

        particle = pp.particle.MuMinusDef()
        for medium in mediums:
            for v in vs:
                # the constructors have not yet been pybinded
                # for param in params:
                param = pp.secondaries.make_secondary(
                    pp.particle.Interaction_Type.mupair, particle, medium)

                buf = [""]

                for energy in energies:
                    rnd1 = pp.RandomGenerator.get().random_double()
                    rnd2 = pp.RandomGenerator.get().random_double()

                    rho = param.calculate_rho(
                        energy, v, medium.components[0], rnd1, rnd2)
                    E1 = 0.5 * v * energy * (1 + rho)
                    E2 = 0.5 * v * energy * (1 - rho)

                    buf.append(particle.name)
                    buf.append(medium.name)
                    buf.append(str(v))
                    buf.append(str(energy))
                    buf.append(str(rnd1))
                    buf.append(str(rnd2))
                    buf.append(str(E1))
                    buf.append(str(E2))
                    buf.append("\n")

                file.write("\t".join(buf))


def main(dir_name):
    create_tables(dir_name, dEdx=True, dNdx=True, stoch=True)
    create_table_rho(dir_name)


if __name__ == "__main__":

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
