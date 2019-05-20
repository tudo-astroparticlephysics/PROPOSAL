import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import pyPROPOSAL as pp

def propagate_particle(propagator,
                       position=[-1e5, 0, 1e4],
                       direction=[1, 0, 0],
                       energy=1e9):
    propagator.particle.position = pp.Vector3D(position[0],
                                         position[1],
                                         position[2])
    tmp_dir = pp.Vector3D(direction[0],
                          direction[1],
                          direction[2])
    tmp_dir.spherical_from_cartesian()
    propagator.particle.direction = tmp_dir

    propagator.particle.energy = energy
    propagator.particle.propagated_distance = 0
    propagator.particle.time = 0

    return propagator.propagate()


def main():
    prop = pp.Propagator(
        particle_def=pp.particle.MuMinusDef.get(),
        config_file="resources/config.json"
    )

    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(3, 1)
    ax1 = fig.add_subplot(gs[:-1])
    ax = fig.add_subplot(gs[-1], sharex=ax1)
    # ax1 = fig.add_subplot(111)

    ax1.plot(np.array([-prop.detector.radius,
                       prop.detector.radius,
                       prop.detector.radius,
                       -prop.detector.radius,
                       -prop.detector.radius]),
             np.array([-prop.detector.height,
                       -prop.detector.height,
                       prop.detector.height,
                       prop.detector.height,
                       -prop.detector.height])/2,
             color='k', label='detector')

    ax1.set_xlabel('x coord. / cm')
    ax1.set_ylabel('z coord. / cm')
    ax1.set_xlim([-1e5, 1e5])

    labels = ['EPair', 'Brems', 'Ioniz', 'NuclInt']

    start_positions = [[-1e5,0,1e4],
                       [-1e5,0,2e4],
                       [-3e4,0,3e4],
                       [1e4,0,4e4],]
    start_energies = [1e9, 3e5, 1e5, 1e5]

    for jdx in range(4):
        secondarys = propagate_particle(prop,
                           position=start_positions[jdx],
                           direction=[1, 0, 0],
                           energy=start_energies[jdx])

        print(prop.particle)
        nsecs = len(secondarys) - 2 # to get rid of the decay
        positions = np.empty((nsecs, 3))
        secs_enrgy = np.empty(nsecs)
        mu_energies = np.empty(nsecs)
        secs_ids = np.empty(nsecs)

        for idx in range(nsecs):
            positions[idx] = np.array([secondarys[idx].position.x,
                                       secondarys[idx].position.y,
                                       secondarys[idx].position.z])
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

        for idx in range(4):
            ax.plot(positions[:,0][secs_ids == idx],
                    secs_enrgy[secs_ids == idx]/1e3,
                    ls='None',
                    marker='.',
                    label=labels[idx])
        ax.plot(positions[:,0], mu_energies/1e3, label=r'$E_{\mu}$')
        ax.axhline(0.5, color='r', label='ecut')
        ax.axvline(prop.particle.entry_point.x, color='g', ls='-', label='entry/exit')
        ax.axvline(prop.particle.exit_point.x, color='g', ls='-')
        ax.axvline(prop.particle.closet_approach_point.x, color='b', ls='dotted', label='closest approach')
        ax.set_yscale('log')
        ax.set_ylabel('Energy / GeV')
        ax.set_xlabel('x coord. / cm')
        # ax.legend()

        plt.subplots_adjust(hspace=.0)
        plt.setp(ax1.get_xticklabels(), visible=False)

        ax1.plot(positions[:,0], positions[:,2], label='muon {}'.format(jdx))
        ax1.plot([prop.particle.entry_point.x, prop.particle.exit_point.x],
                 [prop.particle.entry_point.z, prop.particle.exit_point.z],
                 ls='None', marker='x', label='entry/exit {}'.format(jdx))
        ax1.plot(prop.particle.closet_approach_point.x,
                 prop.particle.closet_approach_point.z,
                 ls='None', marker='x', label='closet approach {}'.format(jdx))
        ax1.plot([prop.particle.entry_point.x, prop.particle.closet_approach_point.x, prop.particle.exit_point.x],
                 [prop.particle.entry_point.z, prop.particle.closet_approach_point.z, prop.particle.exit_point.z],
                 ls='dotted', label='approx line {}'.format(jdx))

    ax1.legend()

    fig.savefig('entry_exit_points.png')

    plt.show()

if __name__ == '__main__':
    main()