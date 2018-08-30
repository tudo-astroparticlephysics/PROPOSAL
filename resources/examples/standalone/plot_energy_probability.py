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
    print("Matplotlib not installed!")

import numpy as np
# np.set_printoptions(threshold='nan')


class ProgressBar(object):

    def __init__(self, loops, bar_lenght=50, start=0, **keywords):

        self._bar_lenght = bar_lenght
        self._bar = []
        self._loops = loops
        self._start = float(start)
        self._current_loop = start

        self._started_process = False
        self._start_time = None

        self._pacman = False

        self._status = ""
        self._text = "\rPercent: [{0}] {1}% Time: {2} Iteration: {3}/{4} {5}"

        self._bar_full = "="
        self._bar_empty = " "

        for key, value in keywords.iteritems():
            if key is "pacman":
                assert type(value) is bool
                self._pacman = value

        if self._pacman:
            self._bar_full = "-"
            self._bar_empty = "o"

            current = self._bar_empty
            for i in range(self._bar_lenght):
                if current is self._bar_empty:
                    current = " "
                    self._bar.append(current)
                else:
                    current = self._bar_empty
                    self._bar.append(current)
        else:
            for i in range(self._bar_lenght):
                self._bar.append(self._bar_empty)

        self._current_pac_state = "C"
        self._current_pac_block = 0

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
            self._status = "Done...\n"

        if self._pacman:
            block = int((self._bar_lenght - 1) * progress)

            if self._current_pac_block < block:
                self._current_pac_block = block
                if self._current_pac_state is "c":
                    self._current_pac_state = "C"
                else:
                    self._current_pac_state = "c"
            else:
                pass

            self._bar[block] = '\033[1m' + "\033[93m" + \
                               self._current_pac_state + '\033[0m'
            self._bar[:block] = block * [self._bar_full]
        else:
            block = int(self._bar_lenght * progress)
            self._bar[:block] = block * [self._bar_full]

        text = self._text.format(
            "".join(self._bar),
            progress*100,
            str(datetime.timedelta(seconds=(time.time() - self._start_time))),
            int(self._current_loop),
            self._loops,
            self._status
        )

        sys.stdout.write(text)
        sys.stdout.flush()


def plot_hist(ax, prim, sec):

    if len(prim) != 0 and len(sec) != 0:
        hist = ax.hist2d(
            prim,
            sec,
            bins=100,
            norm=LogNorm(),
        )
    else:
        hist = ax.hist2d(
            prim,
            sec,
            bins=100,
        )

    ax.set_xlim(xmin=2, xmax=14)
    ax.set_ylim(ymin=-2, ymax=14)
    ax.grid(ls=":", lw=0.2)

    textstr = "count = {:g}\ntotal energy loss = {:.1e} MeV".format(
        sum([sum(x) for x in hist[0]]),
        sum([math.pow(10, x) for x in sec])
    )
    props = dict(facecolor='white', alpha=0.8)
    ax.text(0.03, 0.95, textstr,
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes, fontsize=10, bbox=props)

    plt.colorbar(hist[3], ax=ax, pad=0.01)

    return ax

if __name__ == "__main__":

    import sys

    # =========================================================
    # 	Commandline args
    # =========================================================

    statistics = 100
    config_file = "resources/config_ice.json"

    if len(sys.argv) == 2:
        statistics = int(sys.argv[1])
    elif len(sys.argv) == 3:
        statistics = int(sys.argv[1])
        config_file = sys.argv[2]

    # =========================================================
    #   POPOSAL
    # =========================================================

    prop = pyPROPOSAL.Propagator(
        particle_def=pyPROPOSAL.particle.MuMinusDef.get(),
        config_file=config_file
    )

    mu = prop.particle

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
    n_secondarys = []

    progress = ProgressBar(statistics, pacman=True)
    progress.start()

    for i in range(statistics):
        progress.update()

        mu.position = pyPROPOSAL.Vector3D(0, 0, 0)
        mu.direction = pyPROPOSAL.Vector3D(0, 0, -1)
        mu.energy = math.pow(10, E_max_log)
        mu.propagated_distance = 0

        secondarys = prop.propagate()

        length.append(mu.propagated_distance / 100)
        n_secondarys.append(len(secondarys))

        for sec in secondarys:
            log_sec_energy = math.log10(sec.energy)
            log_energy = math.log10(sec.parent_particle_energy)

            if sec.id == pyPROPOSAL.particle.Data.Epair:
                epair_primary_energy.append(log_energy)
                epair_secondary_energy.append(log_sec_energy)
            if sec.id == pyPROPOSAL.particle.Data.Brems:
                brems_primary_energy.append(log_energy)
                brems_secondary_energy.append(log_sec_energy)
            if sec.id == pyPROPOSAL.particle.Data.DeltaE:
                ioniz_primary_energy.append(log_energy)
                ioniz_secondary_energy.append(log_sec_energy)
            if sec.id == pyPROPOSAL.particle.Data.NuclInt:
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

    inch_to_cm = 2.54
    golden_ratio = 1.61803
    width = 29.7  # cm

    fig = plt.figure(
        figsize=(width / inch_to_cm, width / inch_to_cm / golden_ratio)
    )
    fig.subplots_adjust(wspace=0.1, hspace=0.3)
    fig.suptitle(
        r"{:g} {} with mass {:.1f} MeV and energy $10^{{{}}}$ MeV in {}"
        .format(
            statistics,
            mu.particle_def.name,
            mu.particle_def.mass,
            E_max_log,
            "ice"
        )
    )

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

    fig.text(
        0.95,
        0.7,
        "total energy loss = {:.4e} MeV".format(
            sum([
                sum([math.pow(10, x) for x in ioniz_secondary_energy]),
                sum([math.pow(10, x) for x in brems_secondary_energy]),
                sum([math.pow(10, x) for x in photo_secondary_energy]),
                sum([math.pow(10, x) for x in epair_secondary_energy])
            ])
        ),
        verticalalignment='top',
        horizontalalignment='right',
        fontsize=12,
        rotation='vertical'
    )

    fig.savefig("energy_probability_{}_{}_{}_stats_{}.pdf".format(
        mu.particle_def.name,
        mu.particle_def.mass,
        "ice",
        statistics
    ))

    # =========================================================
    #   Length
    # =========================================================

    fig_length = plt.figure(
        figsize=(width / inch_to_cm, width / inch_to_cm / golden_ratio)
    )
    fig_length.suptitle(
        "propagation lenght of {} with mass {} MeV in {}".format(
            mu.particle_def.name,
            mu.particle_def.mass,
            "ice"
        )
    )

    ax_length = fig_length.add_subplot(111)
    ax_length.hist(length, histtype="step", log=True, bins=100)
    ax_length.set_xlabel(r'range / m')
    ax_length.set_ylabel(r'count')

    fig_length.savefig("lenght_{}_{}_{}_stats_{}.pdf".format(
        mu.particle_def.name,
        mu.particle_def.mass,
        "ice",
        statistics
    ))

    # =========================================================
    #   Secondarys
    # =========================================================

    fig_secondarys = plt.figure(
        figsize=(width / inch_to_cm, width / inch_to_cm / golden_ratio)
    )
    fig_secondarys.suptitle(
        "propagation lenght of {} with mass {} MeV in {}".format(
            mu.particle_def.name,
            mu.particle_def.mass,
            "ice"
        )
    )

    ax_secondarys = fig_secondarys.add_subplot(111)
    ax_secondarys.hist(n_secondarys, histtype="step", log=True, bins=100)
    ax_secondarys.set_xlabel(r'number of interactions')
    ax_secondarys.set_ylabel(r'count')

    fig_secondarys.savefig("secondarys_{}_{}_{}_stats_{}.pdf".format(
        mu.particle_def.name,
        mu.particle_def.mass,
        "ice",
        statistics
    ))

    # write execution time
    # end_time = time.time()

    # with open("time_{}_stats_{}.txt".format(mu.name, statistics), "w") as f:
    #     f.write("execution time: {}".format(end_time - start_time))
