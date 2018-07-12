import pyPROPOSAL as pp
import random


def create_table(dir_name):

    particle_defs = [
        pp.MuMinusDef.get(),
        pp.TauMinusDef.get(),
        pp.EMinusDef.get()
    ]

    cuts = [
        pp.EnergyCutSettings(-1, -1),
        pp.EnergyCutSettings(500, -1),
        pp.EnergyCutSettings(-1, 0.05),
        pp.EnergyCutSettings(500, 0.05)
    ]

    mediums = [
        pp.Medium.Ice(1.0),
        pp.Medium.Hydrogen(1.0),
        pp.Medium.Uranium(1.0)
    ]

    initial_energy = [10**x for x in range(4, 13)]
    final_energy = [energy - 1000 for energy in initial_energy]

    pp.RandomGenerator.get().set_seed(0)

    with open(dir_name + "continous_randomization.txt", "a") as f:

        for particle in particle_defs:
            for medium in mediums:
                for cut in cuts:

                    utility = pp.Utility(
                        particle,
                        medium,
                        cut,
                        pp.UtilityDefinition(),
                        pp.InterpolationDef()
                    )

                    cont_rand = pp.ContinuousRandomizer(
                        utility,
                        pp.InterpolationDef()
                    )

                    buf = [""]

                    for i in range(len(initial_energy)):
                        rnd = random.random()
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

                    print(buf)
                    f.write("\t".join(buf))


def main(dir_name):
    create_table(dir_name)


if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    try:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))
    except OSError:
        print("Directory {} already exists".format(dir_name))

    main(dir_name)
