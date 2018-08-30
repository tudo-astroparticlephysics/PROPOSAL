import pyPROPOSAL as pp

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("Matplotlib not installed!")


if __name__ == "__main__":

    # =========================================================
    # 	Propagate
    # =========================================================

    energy = 1e8  # MeV
    statistics = 100

    sec_def = pp.SectorDefinition()
    sec_def.medium = pp.medium.Ice(1.0)
    sec_def.geometry = pp.geometry.Sphere(pp.Vector3D(), 1e20, 0)
    sec_def.particle_location = pp.ParticleLocation.inside_detector

    sec_def.do_continuous_energy_loss_output = True

    sec_def.scattering_model = pp.scattering.ScatteringModel.HighlandIntegral
    sec_def.crosssection_defs.brems_def.lpm_effect = False
    sec_def.crosssection_defs.epair_def.lpm_effect = False

    sec_def.cut_settings.ecut = 500
    sec_def.cut_settings.vcut = 0.05

    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = "~/.local/share/PROPOSAL/tables"

    prop = pp.Propagator(
            particle_def=pp.particle.MuMinusDef.get(),
            sector_defs=[sec_def],
            detector=pp.geometry.Sphere(pp.Vector3D(), 1e20, 0),
            interpolation_def=interpolation_def
    )

    mu = prop.particle

    pp.RandomGenerator.get().set_seed(1234)

    for i in range(statistics):
        mu.position = pp.Vector3D(0, 0, 0)
        mu.direction = pp.Vector3D(0, 0, -1)
        mu.energy = energy
        mu.propagated_distance = 0

        secondaries = prop.propagate()

        for idx, sec in enumerate(secondaries):
            if idx < 1:
                continue
            if idx > len(secondaries) - 2:
                break
            if sec.id == pp.particle.Data.ContinuousEnergyLoss:
                if secondaries[idx-1].id == pp.particle.Data.ContinuousEnergyLoss or secondaries[idx+1].id == pp.particle.Data.ContinuousEnergyLoss:
                    print("2 Continuous Losses in a row")
                    continue
                energy_diff = secondaries[idx-1].parent_particle_energy - secondaries[idx-1].energy - secondaries[idx+1].parent_particle_energy
                if abs(energy_diff - sec.energy) > 1e-3:
                    print("energy loss differs", energy_diff, sec.energy)
                time_diff = secondaries[idx+1].time - secondaries[idx-1].time
                if abs(time_diff - sec.time) > 1e-3:
                    print("time differs", energy_diff, sec.energy)
