import pyPROPOSAL
import math

try:
    import matplotlib
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
except ImportError:
    print("pylab not installed.  no plots for you.")

try:
    import numpy as np
except ImportError:
    print("numpy not installed.  no plots for you.")

try:
    import multiprocessing as mp
except ImportError:
    print("No Multiprocessing will be done.")

# =========================================================
# 	POPOSAL
# =========================================================


class EnergyProbability(object):

    def __init__(self, ptype, med, cuts, table_dir, stats, processes):

        self.statistics = stats
        self.processes = processes
        self.stats_per_job = self.statistics / self.processes
        self.prop = pyPROPOSAL.Propagator(med, cuts, ptype, table_dir)

        # self.prop = []
        # for p in range(self.cores):
        #     self.prop.append(
        #         pyPROPOSAL.Propagator(med, cuts, ptype, table_dir)
        #     )

        self.epair_primary_energy = mp.Manager().list()
        self.epair_secondary_energy = mp.Manager().list()
        self.brems_primary_energy = mp.Manager().list()
        self.brems_secondary_energy = mp.Manager().list()
        self.ioniz_primary_energy = mp.Manager().list()
        self.ioniz_secondary_energy = mp.Manager().list()
        self.photo_primary_energy = mp.Manager().list()
        self.photo_secondary_energy = mp.Manager().list()

    def propagate(self, stats, process_nr):
        # self.basis.append([i])
        E_max_log = 14

        for i in range(stats):
            # self.prop[process_nr].reset_particle()
            # self.prop[process_nr].particle.energy = math.pow(10, E_max_log)
            # secondarys = self.prop[process_nr].propagate()
            self.prop.reset_particle()
            self.prop.particle.energy = math.pow(10, E_max_log)
            secondarys = self.prop.propagate()

            for sec in secondarys:
                log_sec_energy = math.log10(sec.energy)
                log_energy = math.log10(sec.parent_particle_energy)

                if sec.type is pyPROPOSAL.ParticleType.EPair:
                    self.epair_primary_energy.append(log_energy)
                    self.epair_secondary_energy.append(log_sec_energy)
                if sec.type is pyPROPOSAL.ParticleType.Brems:
                    self.brems_primary_energy.append(log_energy)
                    self.brems_secondary_energy.append(log_sec_energy)
                if sec.type is pyPROPOSAL.ParticleType.DeltaE:
                    self.ioniz_primary_energy.append(log_energy)
                    self.ioniz_secondary_energy.append(log_sec_energy)
                if sec.type is pyPROPOSAL.ParticleType.NuclInt:
                    self.photo_primary_energy.append(log_energy)
                    self.photo_secondary_energy.append(log_sec_energy)

    def run(self):

        jobs = []
        for i in xrange(self.cores):
            print("Job: ", i)
            job = mp.Process(
                target=self.propagate,
                args=(self.stats_per_job, i, )
            )
            jobs.append(job)
            job.start()
        for j in jobs:
            j.join()


def Propagate(
    # statistics,
    prop,
    epair_primary_energy,
    epair_secondary_energy,
    brems_primary_energy,
    brems_secondary_energy,
    ioniz_primary_energy,
    ioniz_secondary_energy,
    photo_primary_energy,
    photo_secondary_energy
):
    """TODO: Docstring for Propagate.
    Returns: TODO

    """

    # statistics = 10
    E_max_log = 14

    # epair_primary_energy = []
    # epair_secondary_energy = []
    #
    # brems_primary_energy = []
    # brems_secondary_energy = []
    #
    # ioniz_primary_energy = []
    # ioniz_secondary_energy = []
    #
    # photo_primary_energy = []
    # photo_secondary_energy = []

    prop.reset_particle()
    prop.particle.energy = math.pow(10, E_max_log)
    secondarys = prop.propagate()

    for sec in secondarys:
        log_sec_energy = math.log10(sec.energy)
        log_energy = math.log10(sec.parent_particle_energy)

        if sec.type is pyPROPOSAL.ParticleType.EPair:
            epair_primary_energy.append(log_energy)
            epair_secondary_energy.append(log_sec_energy)
        if sec.type is pyPROPOSAL.ParticleType.Brems:
            brems_primary_energy.append(log_energy)
            brems_secondary_energy.append(log_sec_energy)
        if sec.type is pyPROPOSAL.ParticleType.DeltaE:
            ioniz_primary_energy.append(log_energy)
            ioniz_secondary_energy.append(log_sec_energy)
        if sec.type is pyPROPOSAL.ParticleType.NuclInt:
            photo_primary_energy.append(log_energy)
            photo_secondary_energy.append(log_sec_energy)

    # for i in statistics:
    #     prop.reset_particle()
    #     prop.particle.energy = math.pow(10, E_max_log)
    #     secondarys = prop.propagate()
    #
    #     for sec in secondarys:
    #         if sec.energy == 0:
    #             continue
    #         if sec.parent_particle_energy == 0:
    #             continue
    #
    #         log_sec_energy = math.log10(sec.energy)
    #         log_energy = math.log10(sec.parent_particle_energy)
    #
    #         if sec.type is pyPROPOSAL.ParticleType.EPair:
    #             epair_primary_energy.append(log_energy)
    #             epair_secondary_energy.append(log_sec_energy)
    #         if sec.type is pyPROPOSAL.ParticleType.Brems:
    #             brems_primary_energy.append(log_energy)
    #             brems_secondary_energy.append(log_sec_energy)
    #         if sec.type is pyPROPOSAL.ParticleType.DeltaE:
    #             ioniz_primary_energy.append(log_energy)
    #             ioniz_secondary_energy.append(log_sec_energy)
    #         if sec.type is pyPROPOSAL.ParticleType.NuclInt:
    #             photo_primary_energy.append(log_energy)
    #             photo_secondary_energy.append(log_sec_energy)

    return (
        epair_primary_energy,
        epair_secondary_energy,
        brems_primary_energy,
        brems_primary_energy,
        ioniz_primary_energy,
        ioniz_secondary_energy,
        photo_primary_energy,
        photo_secondary_energy
    )


def plot_hist(ax, prim, sec):

    hist = ax.hist2d(
        prim,
        sec,
        bins=1000,
        norm=LogNorm(),
    )

    ax.set_xlim(xmin=2, xmax=14)
    ax.set_ylim(ymin=-2, ymax=14)
    ax.grid(ls=":", lw=0.2)

    ax.text(0.05, 0.95, 'count = {:g}'.format(np.sum(hist[0])),
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes, fontsize=10)

    plt.colorbar(hist[3], ax=ax, pad=0.01)

    return ax


if __name__ == "__main__":

    ptype = pyPROPOSAL.ParticleType.EMinus
    mu = pyPROPOSAL.Particle(ptype)
    med = pyPROPOSAL.Medium("standard_rock")
    cuts = pyPROPOSAL.EnergyCutSettings()

    # epair_primary_energy = []
    # epair_secondary_energy = []
    #
    # brems_primary_energy = []
    # brems_secondary_energy = []
    #
    # ioniz_primary_energy = []
    # ioniz_secondary_energy = []
    #
    # photo_primary_energy = []
    # photo_secondary_energy = []

    P = EnergyProbability(ptype, med, cuts, "../../resources/tables", 1000, 3)
    P.run()

    epair_primary_energy = P.epair_primary_energy
    epair_secondary_energy = P.epair_secondary_energy

    brems_primary_energy = P.brems_primary_energy
    brems_secondary_energy = P.brems_secondary_energy

    ioniz_primary_energy = P.ioniz_primary_energy
    ioniz_secondary_energy = P.ioniz_secondary_energy

    photo_primary_energy = P.photo_primary_energy
    photo_secondary_energy = P.photo_secondary_energy

    # statistics = range(10000)
    #
    # pool = multiprocessing.Pool()
    # res = pool.map(
    #     Propagate(
    #         # statistics,
    #         prop,
    #         epair_primary_energy,
    #         epair_secondary_energy,
    #         brems_primary_energy,
    #         brems_secondary_energy,
    #         ioniz_primary_energy,
    #         ioniz_secondary_energy,
    #         photo_primary_energy,
    #         photo_secondary_energy
    #     ),
    #     statistics
    # )

    # pool = multiprocessing.Pool(3)
    # results = pool.map(Propagate(prop))

    # epair_primary_energy = res[0]
    # epair_secondary_energy = res[1]
    # brems_primary_energy = res[2]
    # brems_primary_energy = res[3]
    # ioniz_primary_energy = res[4]
    # ioniz_secondary_energy = res[5]
    # photo_primary_energy = res[6]
    # photo_secondary_energy = res[7]
    #
    # print(res)

    # =========================================================
    # 	Plot
    # =========================================================

    tex_preamble = [
        r"\usepackage{amsmath}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
    ]

    font_size = 10

    params = {
        'backend': 'pdf',
        'font.family': 'serif',
        'font.size': 12,
        'text.usetex': True,
        'text.latex.preamble': tex_preamble,
        'axes.labelsize': font_size,
        'legend.numpoints': 1,
        'legend.shadow': False,
        'legend.fontsize': font_size,
        'xtick.labelsize': font_size,
        'ytick.labelsize': font_size,
        'axes.unicode_minus': True
    }

    plt.rcParams.update(params)

    inch_to_cm = 2.54
    golden_ratio = 1.61803
    width = 29.7  # cm

    fig = plt.figure(
        figsize=(width / inch_to_cm, width / inch_to_cm / golden_ratio)
    )
    fig.subplots_adjust(wspace=0.1, hspace=0.3)

    # =========================================================
    # 	Ionization
    # =========================================================

    ax1 = fig.add_subplot(221)
    ax1 = plot_hist(ax1, ioniz_primary_energy, ioniz_secondary_energy)
    ax1.set_title("ionization")
    ax1.set_ylabel(r'secondary particle energy log($E$/MeV)')

    # =========================================================
    # 	Brems
    # =========================================================

    ax2 = fig.add_subplot(222)
    ax2 = plot_hist(ax2, brems_primary_energy, brems_secondary_energy)
    ax2.set_title("bremsstrahlung")

    # =========================================================
    # 	Photonuclear
    # =========================================================

    ax3 = fig.add_subplot(223)
    ax3 = plot_hist(ax3, photo_primary_energy, photo_secondary_energy)
    ax3.set_title("photonuclear")
    ax3.set_xlabel(r'primary particle energy log($E$/MeV)')
    ax3.set_ylabel(r'secondary particle energy log($E$/MeV)')

    # =========================================================
    # 	Epair
    # =========================================================

    ax4 = fig.add_subplot(224)
    ax4 = plot_hist(ax4, epair_primary_energy, epair_secondary_energy)
    ax4.set_title("pair production")
    ax4.set_xlabel(r'primary particle energy log($E$/MeV)')

    fig.savefig("energy_probability.pdf")
