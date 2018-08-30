
import pyPROPOSAL as pp
import pyPROPOSAL.parametrization as parametrization

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
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

    mu = pp.particle.MuMinusDef.get()
    medium = pp.medium.Ice(1.0)  # With densitiy correction
    cuts = pp.EnergyCutSettings(-1, -1)  # ecut, vcut

    photo = []
    dEdx_photo = []
    energy = [mu.mass + 10**x for x in np.arange(0, 12, 0.2)]

    interpolation_def = pp.InterpolationDef()

    # =========================================================
    # 	Constructor args for parametrizations
    #
    #   - particle
    #   - medium
    #   - cut
    #   - multiplier
    #   - shadowing parametrization
    #   - interpolation definition
    # =========================================================

    param_defs = [
            mu,
            medium,
            cuts,
            1.0,
            parametrization.photonuclear.ShadowButkevichMikhailov(),
            interpolation_def
        ]

    params = [
        parametrization.photonuclear.AbramowiczLevinLevyMaor97Interpolant(
            *param_defs
        ),
        parametrization.photonuclear.ButkevichMikhailovInterpolant(
            *param_defs
        )
    ]

    # =========================================================
    # 	Create x sections out of their parametrizations
    # =========================================================

    crosssections = []

    for param in params:
        crosssections.append(pp.crosssection.PhotoInterpolant(
            param,
            interpolation_def
        ))

    print(crosssections[0])

    # =========================================================
    # 	Calculate DE/dx at the given energies
    # =========================================================

    for cross in crosssections:
        dEdx = []
        for E in energy:
            dEdx.append(cross.calculate_dEdx(E) / E)

        dEdx_photo.append(dEdx)

    # =========================================================
    # 	Plot
    # =========================================================

    fig = plt.figure()

    fig.suptitle(
        "energyloss of {} with mass {} MeV in {}".format(
            mu.name,
            mu.mass,
            medium.name.lower()
        )
    )

    ax = fig.add_subplot(111)
    ax.grid(which='both')

    for dEdx, param in zip(dEdx_photo, params):
        ax.loglog(
            energy,
            dEdx,
            linestyle='-',
            label=param.name
        )

    ax.set_xlabel(r'$E$ / MeV')
    ax.set_ylabel(r'energyloss per energy / $\rm{g}^{-1} \rm{cm}^2$')

    ax.legend(loc='best')

    fig.savefig('photo.pdf')
    plt.show()
