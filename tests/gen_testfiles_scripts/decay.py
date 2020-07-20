import proposal as pp
import numpy as np

def matrix_element_evaluate(particle, products):
    G_F = 1.1663787*1e-2  # MeV
    # G_F = scipy.constants.value(u'Fermi coupling constant') * 1e1

    muon = particle
    electron = products[0]
    numu = products[1]
    nuebar = products[2]

    p1 = muon.energy * nuebar.energy - (muon.momentum * muon.direction) * (nuebar.momentum * nuebar.direction)
    p2 = electron.energy * numu.energy - (electron.momentum * electron.direction) * (numu.momentum * numu.direction)

    return 64 * G_F**2 * p1 * p2

def create_table(dir_name, particle_def, init_energy, decay_products, filename, statistics=int(1e6), NUM_bins=50):
    pp.RandomGenerator.get().set_seed(1234)

    init_particle = pp.particle.DynamicData(particle_def.particle_type)
    init_particle.energy = init_energy
    products = decay_products
    decay_channels = [pp.decay.LeptonicDecayChannelApprox(*products),
                      pp.decay.LeptonicDecayChannel(*products),
                      pp.decay.ManyBodyPhaseSpace(products, matrix_element_evaluate)
                     ]
    decay_names = ["LeptonicDecayChannelApprox", "LeptonicDecayChannel", "ManyBody"]

    histrogram_list = []

    v_max = ( particle_def.mass**2 + products[0].mass**2 ) / (2 * particle_def.mass)
    gamma = init_particle.energy / particle_def.mass
    betagamma = init_particle.momentum / particle_def.mass
    E_max = gamma * v_max + betagamma * np.sqrt(v_max**2 - products[0].mass**2)

    for channel in decay_channels:
        # print(particle_def.name, init_energy, channel)
        prod_0_energies = []
        prod_1_energies = []
        prod_2_energies = []
        for i in range(statistics):
            init_particle.position = pp.Vector3D(0, 0, 0)
            init_particle.direction = pp.Vector3D(0, 0, -1)
            init_particle.energy = init_energy
            init_particle.propagated_distance = 0

            decay_products = channel.decay(particle_def, init_particle)
            for p in decay_products.particles:
                if p.type==products[0].particle_type:
                    prod_0_energies.append(p.energy)
                elif p.type==products[1].particle_type:
                    prod_1_energies.append(p.energy)
                elif p.type==products[2].particle_type:
                    prod_2_energies.append(p.energy)
                else:
                    assert("This should never happen")

        histogram = []
        histogram.append(np.histogram(prod_0_energies, bins=NUM_bins, range=(0, E_max))[0])
        histogram.append(np.histogram(prod_1_energies, bins=NUM_bins, range=(0, E_max))[0])
        histogram.append(np.histogram(prod_2_energies, bins=NUM_bins, range=(0, E_max))[0])

        histrogram_list.append(histogram)

    with open(dir_name + filename, "w") as file:
        buf = [""]
        buf.append(str(statistics))
        buf.append(str(NUM_bins))
        buf.append(str(particle_def.name))
        buf.append(str(init_particle.energy))
        buf.append(str(products[0].name))
        buf.append(str(products[1].name))
        buf.append(str(products[2].name))
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
    create_table(dir_name, pp.particle.MuMinusDef(), pp.particle.MuMinusDef().mass, [pp.particle.EMinusDef(), pp.particle.NuMuDef(), pp.particle.NuEBarDef()], "Decay_MuMinus_rest.txt", int(1e6), 50)
    create_table(dir_name, pp.particle.MuMinusDef(), 1e5, [pp.particle.EMinusDef(), pp.particle.NuMuDef(), pp.particle.NuEBarDef()], "Decay_MuMinus_energy.txt", int(1e6), 50)
    create_table(dir_name, pp.particle.TauMinusDef(), pp.particle.TauMinusDef().mass, [pp.particle.EMinusDef(), pp.particle.NuTauDef(), pp.particle.NuEBarDef()], "Decay_TauMinus_rest.txt", int(1e6), 50)
    create_table(dir_name, pp.particle.TauMinusDef(), 1e5, [pp.particle.EMinusDef(), pp.particle.NuTauDef(), pp.particle.NuEBarDef()], "Decay_TauMinus_energy.txt", int(1e6), 50)

if __name__ == "__main__":

    import os

    dir_name = "TestFiles/"

    if os.path.isdir(dir_name):
        print("Directory {} already exists".format(dir_name))
    else:
        os.makedirs(dir_name)
        print("Directory {} created".format(dir_name))

    main(dir_name)
