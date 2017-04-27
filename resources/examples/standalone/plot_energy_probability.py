import sys
import pyPROPOSAL
import math
import time
import datetime

try:
    import matplotlib
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from mpl_toolkits.axes_grid1 import make_axes_locatable

except ImportError:
    print("pylab not installed.  no plots for you.")


class ProgressBar(object):

    """Docstring for ProgressBar. """

    def __init__(self, loops, bar_lenght=50, start=0):

        self._bar_lenght = bar_lenght
        self._loops = loops
        self._start = float(start)
        self._current_loop = start

        self._started_process = False
        self._start_time = None

        self._status = ""
        self._text = "\rPercent: [{0}] {1}% Time: {2} {3}"
        self._bar_full = "="
        self._bar_empty = " "

    def reset(self):
        self._current_loop = self._start
        self._status = ""
        self._started_process = False

    def start(self):
        self._started_process = True
        self._start_time = time.time()

    def update(self):
        if self._started_process is False:
            print("Pleas start ProgressBar before updating it!")
            return

        self._current_loop += 1.0
        progress = self._current_loop / self._loops

        if progress >= 1.0:
            self._status = "Done..."

        block = int(round(self._bar_lenght * progress))

        text = self._text.format(
            self._bar_full*block + self._bar_empty*(self._bar_lenght - block),
            progress*100,
            str(datetime.timedelta(seconds=(time.time() - self._start_time))),
            self._status
        )

        sys.stdout.write(text)
        sys.stdout.flush()


# =========================================================
#   POPOSAL
# =========================================================

start_time = time.time()

ptype = pyPROPOSAL.ParticleType.MuMinus
mu = pyPROPOSAL.Particle(ptype)
# ptype = pyPROPOSAL.ParticleType.STauMinus
# mu = pyPROPOSAL.Particle(ptype)
mu.mass = 100000
# mu.mass = 10000
# med = pyPROPOSAL.Medium("standard_rock")
# cuts = pyPROPOSAL.EnergyCutSettings()

# prop = pyPROPOSAL.Propagator(med, cuts, ptype, "../../resources/tables")
prop = pyPROPOSAL.Propagator("resources/configuration_IceOnly", mu)

statistics = 1000
E_max_log = 14

epair_primary_energy = []
epair_secondary_energy = []

brems_primary_energy = []
brems_secondary_energy = []

ioniz_primary_energy = []
ioniz_secondary_energy = []

photo_primary_energy = []
photo_secondary_energy = []

length = []
n_daughters = []

progress = ProgressBar(statistics)
progress.start()

for i in range(statistics):
    # update_progress(i / statistics)
    progress.update()
    prop.reset_particle()
    prop.particle.energy = math.pow(10, E_max_log)
    secondarys = prop.propagate()

    length.append(prop.particle.propagated_distance)
    n_daughters.append(len(secondarys))

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

# =========================================================
#   Plot
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


def plot_hist(ax, prim, sec):

    if len(prim) != 0 and len(sec) != 0:
        hist = ax.hist2d(
            prim,
            sec,
            bins=1000,
            norm=LogNorm(),
        )
    else:
        hist = ax.hist2d(
            prim,
            sec,
            bins=1000,
        )

    ax.set_xlim(xmin=2, xmax=14)
    ax.set_ylim(ymin=-2, ymax=14)
    ax.grid(ls=":", lw=0.2)

    textstr = "count = {:g}\ntotal energy loss = {:g} MeV".format(
        sum([sum(x) for x in hist[0]]),
        sum(sec)
    )
    props = dict(facecolor='white', alpha=0.8)
    ax.text(0.03, 0.95, textstr,
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes, fontsize=10, bbox=props)

    # ax.text(0.05, 0.95, 'count = {:g}'.format(sum([
    #     sum(x) for x in hist[0]
    # ])),
    #         verticalalignment='top', horizontalalignment='left',
    #         transform=ax.transAxes, fontsize=10)
    # ax.text(0.05, 0.90, 'total energy loss = {:g}'.format(sum(sec)),
    #         verticalalignment='top', horizontalalignment='left',
    #         transform=ax.transAxes, fontsize=10)

    plt.colorbar(hist[3], ax=ax, pad=0.01)

    return ax

inch_to_cm = 2.54
golden_ratio = 1.61803
width = 29.7  # cm

fig = plt.figure(
    figsize=(width / inch_to_cm, width / inch_to_cm / golden_ratio)
)
fig.subplots_adjust(wspace=0.1, hspace=0.3)
fig.suptitle("{} with mass {} MeV in {}".format(
    prop.particle.name,
    prop.particle.mass,
    prop.collections[0].medium.name.lower()
))

# =========================================================
#   Ionization
# =========================================================

ax1 = fig.add_subplot(221)
ax1 = plot_hist(ax1, ioniz_primary_energy, ioniz_secondary_energy)
ax1.set_title("ionization")
ax1.set_ylabel(r'secondary particle energy log($E$/MeV)')

# =========================================================
#   Brems
# =========================================================

ax2 = fig.add_subplot(222)
ax2 = plot_hist(ax2, brems_primary_energy, brems_secondary_energy)
ax2.set_title("bremsstrahlung")

# =========================================================
#   Photonuclear
# =========================================================

ax3 = fig.add_subplot(223)
ax3 = plot_hist(ax3, photo_primary_energy, photo_secondary_energy)
ax3.set_title("photonuclear")
ax3.set_xlabel(r'primary particle energy log($E$/MeV)')
ax3.set_ylabel(r'secondary particle energy log($E$/MeV)')

# =========================================================
#   Epair
# =========================================================

ax4 = fig.add_subplot(224)
ax4 = plot_hist(ax4, epair_primary_energy, epair_secondary_energy)
ax4.set_title("pair production")
ax4.set_xlabel(r'primary particle energy log($E$/MeV)')

fig.savefig("energy_probability_{}_{}_{}_stats_{}.pdf".format(
    prop.particle.name,
    prop.particle.mass,
    prop.collections[0].medium.name.lower(),
    statistics
))

# =========================================================
#   Lenght
# =========================================================

fig_length = plt.figure(
    figsize=(width / inch_to_cm, width / inch_to_cm / golden_ratio)
)
fig_length.suptitle("propagation lenght of {} with mass {} MeV in {}".format(
    prop.particle.name,
    prop.particle.mass,
    prop.collections[0].medium.name.lower()
))

ax_length = fig_length.add_subplot(111)
ax_length.hist(length, histtype="step", log=True, bins=100)
ax_length.set_xlabel(r'$l_{}(\rm{{m}})$'.format(prop.particle.name))

fig_length.savefig("lenght_{}_{}_{}_stats_{}.pdf".format(
    prop.particle.name,
    prop.particle.mass,
    prop.collections[0].medium.name.lower(),
    statistics
))

# write execution time
end_time = time.time()

with open("time_{}_stats_{}.txt".format(mu.name, statistics), "w") as f:
    f.write("execution time: {}".format(end_time - start_time))
