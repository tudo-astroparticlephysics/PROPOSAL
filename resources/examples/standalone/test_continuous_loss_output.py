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
    sec_def.medium = pp.Medium.Ice(1.0)
    sec_def.geometry_def = pp.Sphere(pp.Vector3D(), 1e20, 0)
    sec_def.particle_location = pp.ParticleLocation.inside_detector

    sec_def.do_continuous_energy_loss_output = True;

    sec_def.scattering_model = pp.ScatteringModel.HighlandIntegral
    sec_def.crosssection_defs.brems_def.lpm_effect = False
    sec_def.crosssection_defs.epair_def.lpm_effect = False

    sec_def.cut_settings.ecut = 500
    sec_def.cut_settings.vcut = 0.05

    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = "~/.local/share/PROPOSAL/tables"

    prop = pp.Propagator(
            particle_def=pp.MuMinusDef.get(),
            sector_defs=[sec_def],
            detector=pp.Sphere(pp.Vector3D(), 1e20, 0),
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
            # print(sec.id)
            if sec.id == pp.Data.ContinuousEnergyLoss:
                if secondaries[idx-1].id == pp.Data.ContinuousEnergyLoss or secondaries[idx+1].id == pp.Data.ContinuousEnergyLoss:
                    print("2 Continuous Losses in a row")
                    continue
                if secondaries[idx+1].id == pp.Data.Particle:
                    # print("a decay")
                    continue

                energy_diff = secondaries[idx-1].parent_particle_energy - secondaries[idx-1].energy - secondaries[idx+1].parent_particle_energy
                if abs(energy_diff - sec.energy) > 1e-3:
                    print("energy loss differs", energy_diff, sec.energy)
                time_diff = secondaries[idx+1].time - secondaries[idx-1].time
                if abs(time_diff - sec.time) > 1e-3:
                    print("time differs", energy_diff, sec.energy)
