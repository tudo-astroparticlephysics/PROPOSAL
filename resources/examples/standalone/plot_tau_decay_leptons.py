import pyPROPOSAL as pp
import time

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("Matplotlib not installed!")

import numpy as np

leptons = [
    pp.particle.NuTauDef.get(),
    pp.particle.NuTauBarDef.get(),
    pp.particle.NuMuDef.get(),
    pp.particle.NuMuBarDef.get(),
    pp.particle.NuEDef.get(),
    pp.particle.NuEBarDef.get(),
    pp.particle.TauMinusDef.get(),
    pp.particle.TauPlusDef.get(),
    pp.particle.MuMinusDef.get(),
    pp.particle.MuPlusDef.get(),
    pp.particle.EMinusDef.get(),
    pp.particle.EPlusDef.get(),
]


def filter_hadr(secondarys):
    prods = [p for p in secondarys if p.id == pp.particle.Data.Particle]
    E = [p.energy for p in prods if p.particle_def not in leptons]
    return sum(E)


def filter_particle(secondarys, particle):
    prods = [p for p in secondarys if p.id == pp.particle.Data.Particle]
    E = [p.energy for p in prods if p.particle_def == particle]
    return sum(E)


def filter_lep(secondarys):
    prods = [p for p in secondarys if p.id == pp.particle.Data.Particle]
    E = [p.energy for p in prods if p.particle_def in leptons]
    return sum(E)


def evaluate(particle, products):
    G_F = 1.1663787*1e-2  # MeV

    mu = particle
    e = products[0]
    numu = products[1]
    nue = products[2]

    p1 = mu.energy * nue.energy - (mu.momentum * mu.direction) * (nue.momentum * nue.direction)
    p2 = e.energy * numu.energy - (e.momentum * e.direction) * (numu.momentum * numu.direction)

    return 64 * G_F**2 * p1 * p2


if __name__ == "__main__":

    # =========================================================
    # 	Save energies
    # =========================================================

    statistics = int(1e5)
    binning = 50

    mu = pp.particle.Particle(pp.particle.MuMinusDef.get())
    mu.direction = pp.Vector3D(0, 0, -1)

    products = [pp.particle.EMinusDef.get(), pp.particle.NuMuDef.get(), pp.particle.NuEBarDef.get()]

    lep_ME = pp.decay.ManyBodyPhaseSpace(products, evaluate)
    lep = pp.decay.LeptonicDecayChannelApprox(*products)

    E_lep_ME = []
    E_lep = []

    passed_time = 0.0

    for i in range(statistics):
        mu.position = pp.Vector3D(0, 0, 0)
        mu.direction = pp.Vector3D(0, 0, -1)
        mu.energy = mu.particle_def.mass
        mu.propagated_distance = 0

        t = time.time()

        d = lep_ME.decay(mu)
        E_lep_ME.append(filter_particle(d, pp.particle.EMinusDef.get()))

        mu.position = pp.Vector3D(0, 0, 0)
        mu.direction = pp.Vector3D(0, 0, -1)
        mu.energy = mu.particle_def.mass
        mu.propagated_distance = 0

        d = lep.decay(mu)
        E_lep.append(filter_particle(d, pp.particle.EMinusDef.get()))

        passed_time += time.time() - t

    print("needed time = {}s".format(passed_time))

    # =========================================================
    # 	Plot
    # =========================================================

    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[2, 1])

    ax = fig.add_subplot(gs[0])

    hist_ME = ax.hist(
        np.array(E_lep) / mu.energy,
        histtype="step",
        # log=True,
        bins=binning,
        color='k',
        ls='-',
        label="decay width"
    )

    hist = ax.hist(
        np.array(E_lep_ME) / mu.energy,
        histtype="step",
        # log=True,
        bins=binning,
        color='k',
        ls='-.',
        label="phase space"
    )

    ax.set_ylabel(r'count')

    ax.legend(loc='upper left')

    # ====[ ratio ]=============================================

    ax = fig.add_subplot(gs[1], sharex=ax)

    diff = hist_ME[0] / hist[0]
    diff[np.isnan(diff)] = 1.0

    bins = hist[1]
    bincenters = 0.5*(bins[1:]+bins[:-1])

    ax.plot(
        bincenters,
        np.nan_to_num(diff),
        drawstyle="steps-mid",
        color='k',
        lw=1.0,
        label="phase space / decay width"
    )

    ax.set_xlabel(r'final energy / MeV')
    ax.set_ylabel(r'ratio')

    ax.legend(loc='lower right')

    ax.axhline(1.0, color='k', lw=0.5, ls='-.')

    fig.tight_layout()

    fig.savefig("leptonic.pdf")

    plt.show()
