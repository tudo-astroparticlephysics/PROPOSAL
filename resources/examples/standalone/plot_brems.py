
import pyPROPOSAL as pp
import pyPROPOSAL.parametrization as parametrization

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

    mu = pp.particle.MuMinusDef.get()
    medium = pp.medium.StandardRock(1.0)  # With densitiy correction
    cuts = pp.EnergyCutSettings(-1, -1)  # ecut, vcut

    dEdx_photo = []
    energy = np.logspace(2, 9, 100)

    interpolation_def = pp.InterpolationDef()

    # =========================================================
    # 	Constructor args for parametrizations
    #
    #   - particle
    #   - medium
    #   - cut
    #   - multiplier
    #   - lpm effect
    #   - interpolation definition
    # =========================================================

    param_defs = [
            mu,
            medium,
            cuts,
            1.0,
            False,
        ]

    params = [
        parametrization.bremsstrahlung.KelnerKokoulinPetrukhin(
            *param_defs
        ),
        parametrization.bremsstrahlung.SandrockSoedingreksoRhode(
            *param_defs
        )
    ]

    # =========================================================
    # 	Create x sections out of their parametrizations
    # =========================================================

    crosssections = []

    for param in params:
        # print(param.name)
        crosssections.append(pp.crosssection.BremsIntegral(
            param,
        ))

    # print(crosssections[0])

    # =========================================================
    # 	Calculate DE/dx at the given energies
    # =========================================================

    for cross in crosssections:
        dEdx = []
        for E in energy:
            dEdx.append(cross.calculate_dEdx(E)/E)

        dEdx_photo.append(dEdx)
        # print(dEdx)

    # =========================================================
    # 	Plot
    # =========================================================

    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[1, 1])

    ax = fig.add_subplot(gs[0])

    for dEdx, param in zip(dEdx_photo, params):
        ax.loglog(
            energy,
            dEdx,
            linestyle='-',
            label="".join([c for c in param.name[1:] if c.isupper()])
        )

    ax.set_ylabel(r'-1/E dEdx / $\rm{g}^{-1} \rm{cm}^2$')

    ax.legend(loc='best')

    # ====[ ratio ]============================================

    ax = fig.add_subplot(gs[1], sharex=ax)

    start = 0
    ax.semilogx(
        energy[start:],
        np.array(dEdx_photo)[1][start:] / np.array(dEdx_photo[0][start:]),
        linestyle='-',
        label="{} / {}".format(
            "".join([c for c in params[1].name[1:] if c.isupper()]),
            "".join([c for c in params[0].name[1:] if c.isupper()])
        )
    )

    ax.xaxis.grid(which='major', ls=":")
    ax.yaxis.grid(which='minor', ls=":")
    ax.set_ylim(top=1.1, bottom=0.7)
    ax.set_xlim(left=1e2, right=1e9)

    ax.set_xlabel(r'$E$ / MeV')
    ax.set_ylabel(r'ratio')

    ax.axhline(1.0, color='k', lw=0.5, ls='-.')

    ax.legend(loc='best')

    fig.savefig('brems.png')
    plt.show()
