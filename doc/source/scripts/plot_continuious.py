import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import proposal as pp

def muons(energy, statistics, vcut, do_continuous_randomization, dist):

    sec_def = pp.SectorDefinition()
    sec_def.medium = pp.medium.StandardRock(1.0)
    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), 1e20, 0)
    sec_def.particle_location = pp.ParticleLocation.inside_detector

    sec_def.scattering_model = pp.scattering.ScatteringModel.Moliere
    sec_def.do_continuous_randomization = do_continuous_randomization

    sec_def.cut_settings.ecut = 0
    sec_def.cut_settings.vcut = vcut

    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = ""

    prop = pp.Propagator(
            particle_def=pp.particle.MuMinusDef.get(),
            sector_defs=[sec_def],
            detector=pp.geometry.Sphere(pp.Vector3D(), 1e20, 0),
            interpolation_def=interpolation_def
    )

    mu = prop.particle

    mu_energies = []

    for i in range(statistics):

        mu.position = pp.Vector3D(0, 0, 0)
        mu.direction = pp.Vector3D(0, 0, -1)
        mu.energy = energy
        mu.propagated_distance = 0

        d = prop.propagate(dist * 100)

        mu_energies.append(mu.energy)

    return mu_energies


if __name__ == "__main__":

    # =========================================================
    # 	Save energies
    # =========================================================

    energy = 1e8
    statistics = int(1e4)
    dist = 300
    binning = 100

    energies_1 = muons(energy, statistics, 0.05, False, dist)
    energies_2 = muons(energy, statistics, 0.0001, False, dist)
    energies_3 = muons(energy, statistics, 0.05, True, dist)

    font_size = 10

    # =========================================================
    # 	Plot energies
    # =========================================================

    binning = 100

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.hist(
        energies_1,
        histtype="step",
        log=True,
        bins=binning,
        label=r"$v_{cut} = 0.05$"
    )

    ax.hist(
        energies_2,
        histtype="step",
        log=True,
        bins=binning,
        label=r"$v_{cut} = 10^{-4}$"
    )

    ax.hist(
        energies_3,
        histtype="step",
        log=True,
        bins=binning,
        label=r"$v_{cut} = 0.05$ with cont."
    )

    ax.set_xlabel(r'finale energie / MeV')
    ax.set_ylabel(r'counts')

    ax.legend(loc='upper left')

    fig.savefig("../pyBindings/figures/mu_continuous.png")

