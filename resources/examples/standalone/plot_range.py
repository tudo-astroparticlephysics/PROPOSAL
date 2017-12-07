import pyPROPOSAL

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("Matplotlib not installed!")


if __name__ == "__main__":

    # =========================================================
    # 	Propagate
    # =========================================================

    energy = 1e8  # MeV
    statistics = 1000

    particle_type = pyPROPOSAL.ParticleType.MuMinus
    medium_type = pyPROPOSAL.MediumType.Ice

    med = pyPROPOSAL.Medium(medium_type, 0.98)
    E = pyPROPOSAL.EnergyCutSettings()
    p = pyPROPOSAL.Propagator(
        medium=med,
        energy_cuts=E,
        particle_type=particle_type,
        path_to_tables="resources/tables",
        moliere=False,
        continuous_rand=True,
        exact_time=True,
        lpm=True,
        scattering_model=2,
    )

    mu_length = list()
    n_secondarys = list()

    for i in range(statistics):
        p.reset_particle()
        p.particle.energy = energy
        d = p.propagate()

        mu_length.append(p.particle.propagated_distance / 100)
        n_secondarys.append(len(d))

    # =========================================================
    # 	Plot lenghts
    # =========================================================

    fig_length = plt.figure()
    ax = fig_length.add_subplot(111)

    ax.hist(mu_length, histtype="step", log=True, bins=50)

    ax.set_title("{} muons with energy {} TeV in {}".format(
        statistics,
        energy / 1e6,
        p.collections[0].medium.name.lower()
    ))
    ax.set_xlabel(r'range / m')
    ax.set_ylabel(r'count')

    fig_length.tight_layout()
    fig_length.savefig("muon_lenghts.pdf")

    # =========================================================
    # 	Plot secondarys
    # =========================================================

    fig_secondarys = plt.figure()
    ax = fig_secondarys.add_subplot(111)

    ax.hist(n_secondarys, histtype="step", log=True, bins=1000)

    ax.set_title("{} muons with energy {} TeV in {}".format(
        statistics,
        energy / 1e6,
        p.collections[0].medium.name.lower()
    ))
    ax.set_xlabel(r'number of interactions')
    ax.set_ylabel(r'count')

    fig_secondarys.tight_layout()
    fig_secondarys.savefig("muon_secondarys.pdf")
