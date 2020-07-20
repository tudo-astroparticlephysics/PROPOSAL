import proposal as pp
import time
from tqdm import tqdm

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt

import numpy as np

leptons = [
    pp.particle.NuTauDef(),
    pp.particle.NuTauBarDef(),
    pp.particle.NuMuDef(),
    pp.particle.NuMuBarDef(),
    pp.particle.NuEDef(),
    pp.particle.NuEBarDef(),
    pp.particle.TauMinusDef(),
    pp.particle.TauPlusDef(),
    pp.particle.MuMinusDef(),
    pp.particle.MuPlusDef(),
    pp.particle.EMinusDef(),
    pp.particle.EPlusDef(),
]


def filter_hadr(secondarys):
    prods = [p for p in secondarys if p.type == pp.particle.Data.Particle]
    E = [p.energy for p in prods if p.particle_def not in leptons]
    return sum(E)


def filter_particle(secondarys, particle_type):
    prods = [p for p in secondarys if p.type == int(particle_type)]
    E = [p.energy for p in prods]
    return sum(E)


def filter_lep(secondarys):
    prods = [p for p in secondarys if p.type == pp.particle.Data.Particle]
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
    binning = np.linspace(0, 0.5, 50)

    pdef = pp.particle.MuMinusDef()
    mu = pp.particle.DynamicData(pdef.particle_type)
    mu.direction = pp.Vector3D(0, 0, -1)

    products = [pp.particle.EMinusDef(), pp.particle.NuMuDef(), pp.particle.NuEBarDef()]

    lep_ME = pp.decay.ManyBodyPhaseSpace(products, evaluate)
    lep = pp.decay.LeptonicDecayChannelApprox(*products)

    E_lep_ME = []
    E_lep = []

    passed_time = 0.0

    for i in tqdm(range(statistics)):
        mu.position = pp.Vector3D(0, 0, 0)
        mu.direction = pp.Vector3D(0, 0, -1)
        mu.energy = pdef.mass
        mu.propagated_distance = 0

        t = time.time()

        d = lep_ME.decay(pdef, mu).particles
        E_lep_ME.append(filter_particle(d, pp.particle.Particle_Type.EMinus))

        mu.position = pp.Vector3D(0, 0, 0)
        mu.direction = pp.Vector3D(0, 0, -1)
        mu.energy = pdef.mass
        mu.propagated_distance = 0

        d = lep.decay(pdef, mu).particles
        E_lep.append(filter_particle(d, pp.particle.Particle_Type.EMinus))

        passed_time += time.time() - t

    print("needed time = {}s".format(passed_time))

    # =========================================================
    # 	Plot
    # =========================================================

    fig = plt.figure()
    gs = GridSpec(2, 1, height_ratios=[2, 1])

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
