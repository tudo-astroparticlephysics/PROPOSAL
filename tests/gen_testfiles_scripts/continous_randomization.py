import os
import proposal as pp
import numpy as np


def create_table(dir_name):

    particle_defs = [
        pp.particle.MuMinusDef(),
        pp.particle.TauMinusDef(),
        pp.particle.EMinusDef()
    ]

    mediums = [
        pp.medium.Ice(),
        pp.medium.Hydrogen(),
        pp.medium.Uranium()
    ]

    cuts = [
        pp.EnergyCutSettings(np.inf, 1, True),
        pp.EnergyCutSettings(500, 1, True),
        pp.EnergyCutSettings(np.inf, 0.05, True),
        pp.EnergyCutSettings(500, 0.05, True)
    ]

    initial_energy = np.logspace(4, 13, num=10)  # MeV
    final_energy = initial_energy - 1000

    pp.RandomGenerator.get().set_seed(0)
    np.random.seed(123)

    with open(dir_name + "continous_randomization.txt", "w") as file:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    args = {
                        "particle_def": particle,
                        "target": medium,
                        "interpolate": False,
                        "cuts": cut
                    }

                    cross = pp.crosssection.make_std_crosssection(**args)
                    cont_rand = pp.make_contrand(cross, False)


                    buf = [""]

                    for i in range(len(initial_energy)):
                        rnd = np.random.random_sample()
                        rand_energy = cont_rand.randomize(
                            initial_energy[i],
                            final_energy[i],
                            rnd
                        )

                        buf.append(str(rnd))
                        buf.append(particle.name)
                        buf.append(medium.name)
                        buf.append(str(cut.ecut))
                        buf.append(str(cut.vcut))
                        buf.append(str(initial_energy[i]))
                        buf.append(str(final_energy[i]))
                        buf.append(str(rand_energy))
                        buf.append("\n")

                    file.write("\t".join(buf))


def main(dir_name):
    create_table(dir_name)

if __name__ == "__main__":

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
