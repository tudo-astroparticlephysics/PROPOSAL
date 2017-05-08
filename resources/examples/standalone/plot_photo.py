
import pyPROPOSAL

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    from cycler import cycler
except ImportError:
    raise ImportError("Matplotlib not installed!")

try:
    import numpy as np
except ImportError:
    raise ImportError(
        "Numpy not installed! Needed to calculate the detector cylinder"
    )

import math


if __name__ == "__main__":

    ptype = pyPROPOSAL.ParticleType.MuMinus
    mu = pyPROPOSAL.Particle(ptype)
    med = pyPROPOSAL.Medium()
    E_cut = pyPROPOSAL.EnergyCutSettings(-1, -1)

    photo = list()
    dEdx_photo = list()
    energy = list()
    energy_range = np.arange(0, 12, 0.2)

    for log_E in energy_range:
        E = mu.mass + math.pow(10, log_E)
        energy.append(E)

    params = [
        pyPROPOSAL.ParametrizationType.PhotoKokoulinShadowBezrukovHard,
        pyPROPOSAL.ParametrizationType.PhotoRhodeShadowBezrukovHard,
        pyPROPOSAL.ParametrizationType.PhotoBezrukovBugaevShadowBezrukovHard,
        pyPROPOSAL.ParametrizationType.PhotoZeusShadowBezrukovHard,
    ]

    for param in params:
        p = pyPROPOSAL.Photonuclear(mu, med, E_cut)
        p.parametrization = param

        p.enable_dEdx_interpolation()
        photo.append(p)

        dEdx = list()
        for E in energy:
            mu.energy = E
            dEdx.append(p.calculate_dEdx() / E)

        dEdx_photo.append(dEdx)

    # =========================================================
    # 	Plot
    # =========================================================

    inch_to_cm = 2.54
    golden_ratio = 1.61803
    width = 29.7  # cm

    fig = plt.figure(
        figsize=(width / inch_to_cm, width / inch_to_cm / golden_ratio)
    )

    fig.suptitle(
        "energyloss of {} with mass {} MeV in {}".format(
            mu.name,
            mu.particle.mass,
            med.name.lower()
        )
    )

    ax = fig.add_subplot(111)
    ax.set_prop_cycle(cycler('color', ['c', 'm', 'y', 'k']))
    ax.grid(which='both')

    for dEdx, param in zip(dEdx_photo, params):
        ax.loglog(
            energy,
            dEdx,
            linestyle='-',
            label=param
        )

    ax.set_xlabel(r'$E$ / MeV')
    ax.set_ylabel(r'energyloss per energy / $\rm{g}^{-1} \rm{cm}^2$')

    ax.legend(loc='best')

    fig.savefig('photo.pdf')
    plt.show()
