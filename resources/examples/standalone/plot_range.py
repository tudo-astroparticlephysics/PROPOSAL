import pyPROPOSAL as pp

try:
    import matplotlib as mpl
    # mpl.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("Matplotlib not installed!")


if __name__ == "__main__":

    # =========================================================
    # 	Propagate
    # =========================================================

    energy = 1e8  # MeV
    statistics = 1000

    geo_def = pp.GeometryDefinition()
    geo_def.shape = pp.Shape.Sphere
    geo_def.outer_radius = 1e20
    geo_def.inner_radius = 0

    med_def = pp.MediumDefinition()
    med_def.type = pp.MediumType.Ice

    sec_def = pp.SectorDefinition()
    sec_def.medium_def = med_def
    sec_def.geometry_def = geo_def
    sec_def.particle_location = pp.ParticleLocation.inside_detector

    sec_def.scattering_model = pp.ScatteringModel.moliere

    sec_def.e_cut = 500
    sec_def.v_cut = 0.05

    interpolation_def = pp.InterpolationDef()

    prop = pp.Propagator(
        particle_def=pp.MuMinusDef.get(),
        config_file="resources/config_ice.json"
        # sector_defs=[sec_def],
        # detector=pp.Sphere(pp.Vector3D(), 0, 1e20),
        # interpolation_def=interpolation_def
    )

    mu = prop.particle

    mu_length = []
    n_secondarys = []

    for i in range(statistics):

        mu.position = pp.Vector3D(0, 0, 0)
        mu.direction = pp.Vector3D(0, 0, -1)
        mu.energy = energy
        mu.propagated_distance = 0

        d = prop.propagate()

        mu_length.append(mu.propagated_distance / 100)
        n_secondarys.append(len(d))

    # =========================================================
    # 	Plot lenghts
    # =========================================================

    fig_length = plt.figure()
    ax = fig_length.add_subplot(111)

    ax.hist(mu_length, histtype="step", log=True, bins=50)

    ax.set_title("{} muons with energy {} TeV".format(
        statistics,
        energy / 1e6
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

    ax.hist(n_secondarys, histtype="step", log=True, bins=50)

    ax.set_title("{} muons with energy {} TeV".format(
        statistics,
        energy / 1e6
    ))
    ax.set_xlabel(r'number of interactions')
    ax.set_ylabel(r'count')

    fig_secondarys.tight_layout()
    fig_secondarys.savefig("muon_secondarys.pdf")
