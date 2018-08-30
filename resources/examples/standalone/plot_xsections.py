
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

    dEdx = []
    # energy = [mu.mass + 10**x for x in np.arange(0, 12, 0.2)]
    energy = np.logspace(2, 12, 100)

    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = "~/.local/share/PROPOSAL/tables"

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

    photo_param = parametrization.photonuclear.AbramowiczLevinLevyMaor97Interpolant(
        mu,
        medium,
        cuts,
        1.0,
        parametrization.photonuclear.ShadowButkevichMikhailov(),
        interpolation_def
    )

    # =========================================================
    # 	Create x sections out of their parametrizations
    # =========================================================

    crosssections = []

    brems = pp.crosssection.BremsInterpolant(
        parametrization.bremsstrahlung.KelnerKokoulinPetrukhin(
            mu,
            medium,
            cuts,
            1.0,
            True
        ),
        interpolation_def
        )

    photo = pp.crosssection.PhotoInterpolant(
        parametrization.photonuclear.AbramowiczLevinLevyMaor97Interpolant(
            mu,
            medium,
            cuts,
            1.0,
            parametrization.photonuclear.ShadowButkevichMikhailov(),
            interpolation_def
        ),
        interpolation_def
        )

    epair = pp.crosssection.EpairInterpolant(
        parametrization.pairproduction.KelnerKokoulinPetrukhinInterpolant(
            mu,
            medium,
            cuts,
            1.0,
            True,
            interpolation_def
        ),
        interpolation_def
        )

    ioniz = pp.crosssection.IonizInterpolant(
        parametrization.ionization.Ionization(
            mu,
            medium,
            cuts,
            1.0
        ),
        interpolation_def
        )

    # =========================================================
    # 	Calculate DE/dx at the given energies
    # =========================================================

    dEdx_brems = []
    dEdx_photo = []
    dEdx_epair = []
    dEdx_ioniz = []
    dEdx_sum = []

    for E in energy:
        b = brems.calculate_dEdx(E)
        dEdx_brems.append(b)
        p = photo.calculate_dEdx(E)
        dEdx_photo.append(p)
        e = epair.calculate_dEdx(E)
        dEdx_epair.append(e)
        i = ioniz.calculate_dEdx(E)
        dEdx_ioniz.append(i)

        dEdx_sum.append(b + p + e + i)

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

    ax.loglog(
        energy,
        dEdx_brems,
        linestyle='-',
        label="bremsstrahlung"
    )
    ax.loglog(
        energy,
        dEdx_photo,
        linestyle='-',
        label="photonuclear"
    )
    ax.loglog(
        energy,
        dEdx_epair,
        linestyle='-',
        label="epair production"
    )
    ax.loglog(
        energy,
        dEdx_ioniz,
        linestyle='-',
        label="ionization"
    )

    ax.loglog(
        energy,
        dEdx_sum,
        linestyle='-',
        label="sum"
    )

    ax.set_xlabel(r'$E$ / MeV')
    ax.set_ylabel(r'energyloss $dE/dx$ / $\rm{g}^{-1} \rm{cm}^2$')

    ax.legend(loc='best')

    fig.savefig('xsections.pdf')
    plt.show()
