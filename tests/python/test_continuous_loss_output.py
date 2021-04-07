import proposal as pp
import os
from pytest import approx


#ContinuousEnergyLoss = pp.particle.Interaction_Type.ContinuousEnergyLoss
#table_path = os.path.expanduser("~/.local/share/PROPOSAL/tables")


def test_proposal():

    # =========================================================
    # 	Propagate
    # =========================================================

    # initialize propagator

    args = {
        "particle_def": pp.particle.MuMinusDef(),
        "target": pp.medium.Ice(),
        "interpolate": True,
        "cuts": pp.EnergyCutSettings(500, 0.05)
    }

    cross = pp.crosssection.make_std_crosssection(**args)

    collection = pp.PropagationUtilityCollection()
    collection.displacement = pp.make_displacement(cross, True)
    collection.interaction = pp.make_interaction(cross, True)
    collection.time = pp.make_time(cross, args["particle_def"], True)
    utility = pp.PropagationUtility(collection = collection)

    detector = pp.geometry.Sphere(pp.Cartesian3D(0,0,0), 1e20)
    density_distr = pp.density_distribution.density_homogeneous(args["target"].mass_density)

    prop = pp.Propagator(args["particle_def"], [(detector, utility, density_distr)])

    # intialize initial state

    statistics = 100

    init_state = pp.particle.ParticleState()
    init_state.energy = 1e8 # initial energy in MeV
    init_state.position = pp.Cartesian3D(0, 0, 0)
    init_state.direction = pp.Cartesian3D(0, 0, -1)
    init_state.time = 0

    pp.RandomGenerator.get().set_seed(1234)
    for i in range(statistics):

        output = prop.propagate(init_state)
        energies = output.track_energies()
        times = output.track_times()

        E_old = init_state.energy
        t_old = init_state.time

        for E, t in zip(energies, times):
        	energy_diff = E_old - E
        	time_diff = t - t_old
        	E_old = E
        	t_old = t

        	assert (energy_diff >= 0)
        	assert (time_diff >= 0)

if __name__ == '__main__':
    test_proposal()
