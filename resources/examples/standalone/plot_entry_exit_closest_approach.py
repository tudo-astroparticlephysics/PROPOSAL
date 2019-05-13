import numpy as np
from matplotlib import pyplot as plt
import pyPROPOSAL as pp

prop = pp.Propagator(
    particle_def=pp.particle.MuMinusDef.get(),
    config_file="resources/config.json"
)

mu = prop.particle
start_position = pp.Vector3D(1e5, 0, 1e2)
mu.position = start_position
direction = pp.Vector3D(-1, 0, 0)
direction.spherical_from_cartesian()
mu.direction = direction
start_energy = 1e9
mu.energy = start_energy
mu.propagated_distance = 0
mu.time = 0

secondarys = prop.propagate()

nsecs = len(secondarys) -3
secs_dists = np.empty(nsecs)
secs_enrgy = np.empty(nsecs)
mu_energies = np.empty(nsecs)
secs_ids = np.empty(nsecs)

for idx in range(nsecs):
    secs_dists[idx] = (start_position - secondarys[idx].position).magnitude()
    secs_enrgy[idx] = secondarys[idx].energy
    mu_energies[idx] = secondarys[idx].parent_particle_energy
    if secondarys[idx].id == pp.particle.Data.Epair:
        secs_ids[idx] = 0
    elif secondarys[idx].id == pp.particle.Data.Brems:
        secs_ids[idx] = 1
    elif secondarys[idx].id == pp.particle.Data.DeltaE:
        secs_ids[idx] = 2
    elif secondarys[idx].id == pp.particle.Data.NuclInt:
        secs_ids[idx] = 3

print(mu)
dist_start_detector = (prop.detector.position - start_position).magnitude()

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(secs_dists[secs_ids == 0], secs_enrgy[secs_ids == 0]/1e3, 'b.', label='EPair')
ax.plot(secs_dists[secs_ids == 1], secs_enrgy[secs_ids == 1]/1e3, 'r.', label='Brems')
ax.plot(secs_dists[secs_ids == 2], secs_enrgy[secs_ids == 2]/1e3, 'm.', label='Ioniz')
ax.plot(secs_dists[secs_ids == 3], secs_enrgy[secs_ids == 3]/1e3, 'g.', label='NuclI')
ax.plot(secs_dists, mu_energies/1e3)
ax.axhline(0.5, color='r', label='ecut')

ax.axvline((start_position - mu.entry_point).magnitude(), color='k', ls='-', label='detector')
ax.axvline((start_position - mu.closet_approach_point).magnitude(), color='k', ls='--', label='closest approach')
ax.axvline((start_position - mu.exit_point).magnitude(), color='k', ls='-')
ax.set_yscale('log')
ax.set_ylim([0.1, start_energy/1e3])
ax.set_ylabel('Energy / GeV')
ax.set_xlabel('Distance / cm')
ax.legend()
plt.show()
