import proposal as pp
import matplotlib.pyplot as plt
from tqdm import tqdm


if __name__ == "__main__":

    # =========================================================
    # 	Propagate
    # =========================================================

    energy = 1e7  # MeV
    statistics = 10000

    sec_def = pp.SectorDefinition()
    sec_def.medium = pp.medium.Ice(1.0)
    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), 1e20, 0)
    sec_def.particle_location = pp.ParticleLocation.inside_detector

    sec_def.scattering_model = pp.scattering.ScatteringModel.Moliere
    sec_def.crosssection_defs.brems_def.lpm_effect = False
    sec_def.crosssection_defs.epair_def.lpm_effect = False

    sec_def.cut_settings.ecut = 500
    sec_def.cut_settings.vcut = 0.05

    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = "~/.local/share/PROPOSAL/tables"
    interpolation_def.path_to_tables_readonly = "~/.local/share/PROPOSAL/tables"

    prop = pp.Propagator(
            particle_def=pp.particle.MuMinusDef(),
            sector_defs=[sec_def],
            detector=pp.geometry.Sphere(pp.Vector3D(), 1e20, 0),
            interpolation_def=interpolation_def
    )

    mu = pp.particle.DynamicData(pp.particle.MuMinusDef().particle_type)
    mu.position = pp.Vector3D(0, 0, 0)
    mu.direction = pp.Vector3D(0, 0, -1)
    mu.energy = energy
    mu.propagated_distance = 0

    mu_length = []
    n_secondarys = []

    for i in tqdm(range(statistics)):

        d = prop.propagate(mu).particles

        mu_length.append(d[-1].position.magnitude() / 100)
        n_secondarys.append(len(d))

    # =========================================================
    # 	Plot lenghts
    # =========================================================

    fig_length = plt.figure()
    ax = fig_length.add_subplot(111)

    ax.hist(mu_length, histtype="step", log=True, bins=50)

    ax.set_title("{} {}'s with energy {} TeV".format(
        statistics,
        prop.particle_def.name,
        energy / 1e6
    ))
    ax.set_xlabel(r'range / m')
    ax.set_ylabel(r'count')

    fig_length.tight_layout()
    fig_length.savefig(
        "{}_lenghts.pdf".format(prop.particle_def.name)
    )

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
    fig_secondarys.savefig(
        "{}_secondaries.pdf".format(prop.particle_def.name)
    )
