import pyPROPOSAL as pp

import matplotlib.pyplot as plt

from tqdm import tqdm

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
    interpolation_def.path_to_tables_readonly = "~/.local/share/PROPOSAL/tables"

    mu_def = pp.particle.MuMinusDef()
    prop = pp.Propagator(
            particle_def=mu_def,
            sector_defs=[sec_def],
            detector=pp.geometry.Sphere(pp.Vector3D(), 1e20, 0),
            interpolation_def=interpolation_def
    )

    mu = pp.particle.DynamicData(mu_def.particle_type)

    mu.position = pp.Vector3D(0, 0, 0)
    mu.direction = pp.Vector3D(0, 0, -1)
    mu.energy = energy
    mu.propagated_distance = 0

    pp.RandomGenerator.get().set_seed(1234)
    for i in tqdm(range(statistics)):

        secondaries = prop.propagate(mu).particles

        for idx, sec in enumerate(secondaries):
            if idx < 1:
                continue
            if idx > len(secondaries) - 2:
                break
            if sec.type == int(pp.particle.Interaction_Type.ContinuousEnergyLoss):
                if secondaries[idx-1].type == int(pp.particle.Interaction_Type.ContinuousEnergyLoss) or secondaries[idx+1].type == int(pp.particle.Interaction_Type.ContinuousEnergyLoss):
                    print("2 Continuous Losses in a row")
                    continue
                energy_diff = secondaries[idx-1].energy - secondaries[idx+1].parent_particle_energy
                continuou_energy_lost = sec.parent_particle_energy - sec.energy
                if abs(energy_diff - continuou_energy_lost) > 1e-3:
                    print("energy loss differs", energy_diff, continuou_energy_lost)
                time_diff = secondaries[idx+1].time - secondaries[idx-1].time
                if abs(time_diff - sec.time) > 1e-3:
                    print("time differs", time_diff, sec.time)
