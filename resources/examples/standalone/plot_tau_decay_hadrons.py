import pyPROPOSAL as pp
import time

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import colors
except ImportError:
    raise ImportError("Matplotlib not installed!")

import numpy as np

data_formfactor = np.array([
    [0.6, 0.7, 0.75, 1.0, 1.6, 2.45],
    [0.493, 0.471, 0.407, 0.351, 0.243, 0.167]
])


form_factor = np.poly1d(np.polyfit(1e6 * data_formfactor[0], data_formfactor[1], 3))

# def form_factor(x):
#     term = 1. + 0.0306 * x + (0.0194 * x**3) / (1. + x)
#     return np.exp(-1.171 * x**(0.536)) * term

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


def scalar_product(a, b):
    return a[0] * b[0] - a[1] * b[1]


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def evaluate(particle, products):
    G_F = 1.1663787*1e-2  # MeV
    V_ud = 0.97427

    tau = particle
    pi1 = products[0]
    pi2 = products[1]
    nu = products[2]

    P = (tau.energy, tau.momentum * tau.direction)
    p1 = (pi1.energy, pi1.momentum * pi1.direction)
    p2 = (pi2.energy, pi2.momentum * pi2.direction)
    k = (nu.energy, nu.momentum * nu.direction)

    p12 = add(p1, (-1.0 * p2[0], -1.0 * p2[1]))
    q = add(p1, p2)
    q2 = scalar_product(q, q)
    term1 = scalar_product(k, p12)
    term2 = scalar_product(P, p12)
    term3 = scalar_product(P, k) * scalar_product(p12, p12)

    Y = 2. * term1 * term2 - term3

    return 4. * (G_F * V_ud * form_factor(q2))**2 * Y

if __name__ == "__main__":

    # =========================================================
    # 	Save energies
    # =========================================================

    statistics = int(1e5)
    binning = 50

    tau = pp.particle.Particle(pp.particle.TauMinusDef.get())
    tau.direction = pp.Vector3D(0, 0, -1)

    products = [
        pp.particle.Pi0Def.get(),
        pp.particle.PiMinusDef.get(),
        pp.particle.NuTauDef.get()
    ]

    products_particles = [pp.particle.Particle(p) for p in products]
    for p in products_particles:
        p.direction = pp.Vector3D(0, 0, -1)
        p.energy = 1e2

    print(evaluate(tau, products_particles))

    ME = pp.decay.ManyBodyPhaseSpace(products, evaluate)

    E_lep_pi0 = []
    E_lep_pim = []
    E_lep_nu = []

    s1 = []
    s2 = []

    passed_time = 0.0

    for i in range(statistics):
        tau.position = pp.Vector3D(0, 0, 0)
        tau.direction = pp.Vector3D(0, 0, -1)
        tau.energy = tau.particle_def.mass
        tau.propagated_distance = 0

        t = time.time()

        d = ME.decay(tau)

        E_lep_pi0.append(filter_particle(d, pp.particle.Pi0Def.get()))
        E_lep_pim.append(filter_particle(d, pp.particle.PiMinusDef.get()))
        E_lep_nu.append(filter_particle(d, pp.particle.NuTauDef.get()))

        p1 = add((d[0].energy, d[0].momentum * d[0].direction), (d[1].energy, d[1].momentum * d[1].direction))
        p2 = add((d[1].energy, d[1].momentum * d[1].direction), (d[2].energy, d[2].momentum * d[2].direction))

        s1.append(
            scalar_product(p1, p1)
        )

        s2.append(
            scalar_product(p2, p2)
        )

        passed_time += time.time() - t

    print("needed time = {}s".format(passed_time))

    # =========================================================
    # 	Plot
    # =========================================================

    fig = plt.figure()
    ax = fig.add_subplot(111)

    hist = ax.hist(
        np.array(E_lep_pi0) / tau.energy,
        histtype="step",
        bins=binning,
        color='k',
        ls='-',
        label=r"$\pi^0$"
    )

    hist = ax.hist(
        np.array(E_lep_pim) / tau.energy,
        histtype="step",
        bins=binning,
        color='k',
        ls='-.',
        label=r"$\pi^-$"
    )

    hist = ax.hist(
        np.array(E_lep_nu) / tau.energy,
        histtype="step",
        bins=binning,
        color='k',
        ls=':',
        label=r"$\nu$"
    )

    ax.set_ylabel(r'count')
    ax.legend(loc='upper left')

    fig.tight_layout()

    fig.savefig("hadronic.pdf")

    # =========================================================
    # 	Dalitz Plot
    # =========================================================

    fig = plt.figure()
    ax = fig.add_subplot(111)

    hist = ax.hist2d(
        s1,
        s2,
        bins=500,
        norm=colors.LogNorm()
    )

    plt.colorbar(hist[3], ax=ax, pad=0.01)

    ax.set_xlabel(r'$s_1 = (p_1 + p_2)^2$ / $\rm{GeV}^2$')
    ax.set_ylabel(r'$s_2 = (p_2 + p_3)^2$ / $\rm{GeV}^2$')

    fig.tight_layout()

    fig.savefig("dalitz.pdf")

    plt.show()
