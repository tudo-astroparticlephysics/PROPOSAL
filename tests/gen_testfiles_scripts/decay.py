import pyPROPOSAL as pp
import numpy as np 

def create_table(dir_name):
    """TODO: Docstring for create_table.
    Returns: TODO

    """

    statistics = int(1e6)
    NUM_bins = 50

    pp.RandomGenerator.get().set_seed(1234)

    mu = pp.particle.Particle(pp.particle.MuMinusDef.get())
    products = [pp.particle.EMinusDef.get(), pp.particle.NuMuDef.get(), pp.particle.NuEBarDef.get()]
    decay_channels = [pp.decay.LeptonicDecayChannelApprox(*products),
                      pp.decay.LeptonicDecayChannel(*products) 
                     ]
    decay_names = ["LeptonicDecayChannelApprox", "LeptonicDecayChannel"]

    histrogram_list = []


    for channel in decay_channels:
        prod_0_energies = []
        prod_1_energies = []
        prod_2_energies = []
        for i in range(statistics):
            mu.position = pp.Vector3D(0, 0, 0)
            mu.direction = pp.Vector3D(0, 0, -1)
            mu.energy = mu.particle_def.mass
            mu.propagated_distance = 0
            
            decay_products = channel.decay(mu)
            for p in decay_products:
                if p.particle_def==products[0]:
                    prod_0_energies.append(p.energy)
                elif p.particle_def==products[1]:
                    prod_1_energies.append(p.energy)
                elif p.particle_def==products[2]:
                    prod_2_energies.append(p.energy)
                else:
                    assert("This should never happen")
        
        histogram = []
        histogram.append(np.histogram(prod_0_energies, bins=NUM_bins, range=(0, mu.particle_def.mass/2))[0])
        histogram.append(np.histogram(prod_1_energies, bins=NUM_bins, range=(0, mu.particle_def.mass/2))[0])
        histogram.append(np.histogram(prod_2_energies, bins=NUM_bins, range=(0, mu.particle_def.mass/2))[0])

        histrogram_list.append(histogram)

    with open(dir_name + "Decay_Leptonic.txt", "a") as file:
        buf = [""]
        buf.append(str(statistics))
        buf.append(str(NUM_bins))
        buf.append("\n")
        file.write("\t".join(buf))

        for name, hist in zip(decay_names, histrogram_list):
            buf = [""]
            buf.append(name)
            buf.append("\n")
            for prod in hist:
                for bin_value in prod:
                    buf.append(str(bin_value))
                buf.append("\n")

            file.write("\t".join(buf))



def main(dir_name):
    create_table(dir_name)


if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
