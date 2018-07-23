
import pyPROPOSAL as pp
import pyPROPOSAL.Parametrization as Parametrization

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
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

    mu = pp.MuMinusDef.get()
    medium = pp.Medium.StandardRock(1.0)  # With densitiy correction
    cuts = pp.EnergyCutSettings(-1, -1)  # ecut, vcut

    dEdx_photo = []
    # energy = [mu.mass + 10**x for x in np.arange(0, 12, 0.2)]
    # energy = [mu.mass + 10**x for x in np.arange(0, 3.5, 0.01)]
    # energy = np.linspace(4e2, 3e3, 100)
    energy = np.logspace(2, 9, 100)
    # energy =  np.linspace(0, 1, 100)

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

    # param_defs = [
    #         mu,
    #         medium,
    #         cuts,
    #         1.0,
    #         False,
    #         interpolation_def
    #     ]
    #
    # params = [
    #     Parametrization.EpairProduction.KelnerInterpolant(
    #         *param_defs
    #     ),
    #     Parametrization.EpairProduction.RhodeSandrockSoedingreksoInterpolant(
    #         *param_defs
    #     )
    # ]

    param_defs = [
            mu,
            medium,
            cuts,
            1.0,
            False,
        ]

    params = [
        Parametrization.EpairProduction.Kelner(
            *param_defs
        ),
        Parametrization.EpairProduction.RhodeSandrockSoedingrekso(
            *param_defs
        )
    ]

    # =========================================================
    # 	Create x sections out of their parametrizations
    # =========================================================

    crosssections = []

    # for param in params:
    #     print(param.name)
    #     crosssections.append(pp.CrossSection.EpairInterpolant(
    #         param,
    #         interpolation_def
    #     ))

    for param in params:
        print(param.name)
        crosssections.append(pp.CrossSection.EpairIntegral(
            param,
        ))

    print(crosssections[0])

    # =========================================================
    # 	Calculate DE/dx at the given energies
    # =========================================================

    for cross in crosssections:
        dEdx = []
        for E in energy:
            dEdx.append(cross.calculate_dEdx(E))

        print(dEdx)
        dEdx_photo.append(dEdx)

    # for param in params:
    #     dEdx = []
    #     for E in energy:
    #         dEdx.append(param.differential_crosssection(1e5, E))
    #
    #     print(dEdx)
    #
    #     dEdx_photo.append(dEdx)

    # =========================================================
    # 	Plot
    # =========================================================

    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[2, 1])

    ax = fig.add_subplot(gs[0])

    # ax.grid(which='both')

    for dEdx, param in zip(dEdx_photo, params):
        ax.loglog(
            energy,
            dEdx,
            linestyle='-',
            label=param.name
        )

    ax.set_ylabel(r'dEdx / $\rm{MeV}\rm{g}^{-1} \rm{cm}^2$')

    ax.legend(loc='best')

    # ====[ ratio ]============================================

    ax = fig.add_subplot(gs[1], sharex=ax)

    start = 0
    ax.semilogx(
        energy[start:],
        np.array(dEdx_photo)[1][start:] / np.array(dEdx_photo[0][start:]),
        linestyle='-',
        label="RSS / KKP"
    )

    ax.xaxis.grid(which='major', ls=":")
    ax.yaxis.grid(which='minor', ls=":")
    ax.set_ylim(top=1.03, bottom=0.97)
    ax.set_xlim(left=1e2, right=1e9)

    ax.yaxis.set_minor_formatter(ticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))

    ax.set_xlabel(r'$E$ / MeV')
    ax.set_ylabel(r'RSS / KPP')

    ax.axhline(1.0, color='k', lw=0.5, ls='-.')

    ax.legend(loc='best')

    fig.savefig('epair.png')
    plt.show()
