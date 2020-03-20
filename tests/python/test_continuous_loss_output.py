import proposal as pp
import os
from pytest import approx


ContinuousEnergyLoss = pp.particle.Interaction_Type.ContinuousEnergyLoss
table_path = os.path.expanduser("~/.local/share/PROPOSAL/tables")


def test_proposal():

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
    interpolation_def.path_to_tables = table_path
    interpolation_def.path_to_tables_readonly = table_path

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
    for i in range(statistics):

        secondaries = prop.propagate(mu).particles

        for idx, sec in enumerate(secondaries[1:-1]):
            if sec.type == int(ContinuousEnergyLoss):
                if (
                    secondaries[idx - 1].type == int(ContinuousEnergyLoss)
                    or secondaries[idx + 1].type == int(ContinuousEnergyLoss)
                ):
                    print("2 Continuous Losses in a row")
                    continue

                energy_diff = secondaries[idx - 1].energy - secondaries[idx + 1].parent_particle_energy
                continuou_energy_lost = sec.parent_particle_energy - sec.energy
                assert energy_diff == approx(continuou_energy_lost, abs=1e-3), "Energy differs"

                time_diff = secondaries[idx + 1].time - secondaries[idx - 1].time
                assert time_diff == approx(sec.time, abs=1e-3), "Time differs"


if __name__ == '__main__':
    test_proposal()
