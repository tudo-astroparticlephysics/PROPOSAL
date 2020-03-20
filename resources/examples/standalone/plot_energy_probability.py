import sys
import proposal as pp
import math
import time

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import numpy as np
from tqdm import tqdm


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

    ax.set_xlim(left=2, right=14)
    ax.set_ylim(bottom=-2, top=14)
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

    mu_def = pp.particle.MuMinusDef()
    prop = pp.Propagator(
        particle_def=mu_def,
        config_file=config_file
    )

    E_max_log = 14

    mu = pp.particle.DynamicData(mu_def.particle_type)

    mu.position = pp.Vector3D(0, 0, 0)
    mu.direction = pp.Vector3D(0, 0, -1)
    mu.energy = math.pow(10, E_max_log)
    mu.propagated_distance = 0
    mu.time = 0

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

    for i in tqdm(range(statistics)):

        secondarys = prop.propagate(mu).particles

        length.append(secondarys[-1].propagated_distance / 100)
        n_secondarys.append(len(secondarys))

        for sec in secondarys:
            log_sec_energy = math.log10(sec.parent_particle_energy-sec.energy)
            log_energy = math.log10(sec.parent_particle_energy)

            if sec.type == int(pp.particle.Interaction_Type.Epair):
                epair_primary_energy.append(log_energy)
                epair_secondary_energy.append(log_sec_energy)
            elif sec.type == int(pp.particle.Interaction_Type.Brems):
                brems_primary_energy.append(log_energy)
                brems_secondary_energy.append(log_sec_energy)
            elif sec.type == int(pp.particle.Interaction_Type.DeltaE):
                ioniz_primary_energy.append(log_energy)
                ioniz_secondary_energy.append(log_sec_energy)
            elif sec.type == int(pp.particle.Interaction_Type.NuclInt):
                photo_primary_energy.append(log_energy)
                photo_secondary_energy.append(log_sec_energy)
            # else:
                # print("No valid interaction id {}".format(sec.id))

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
            mu_def.name,
            mu_def.mass,
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
        mu_def.name,
        mu_def.mass,
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
            mu_def.name,
            mu_def.mass,
            "ice"
        )
    )

    ax_length = fig_length.add_subplot(111)
    ax_length.hist(length, histtype="step", log=True, bins=100)
    ax_length.set_xlabel(r'range / m')
    ax_length.set_ylabel(r'count')

    fig_length.savefig("lenght_{}_{}_{}_stats_{}.pdf".format(
        mu_def.name,
        mu_def.mass,
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
            mu_def.name,
            mu_def.mass,
            "ice"
        )
    )

    ax_secondarys = fig_secondarys.add_subplot(111)
    ax_secondarys.hist(n_secondarys, histtype="step", log=True, bins=100)
    ax_secondarys.set_xlabel(r'number of interactions')
    ax_secondarys.set_ylabel(r'count')

    fig_secondarys.savefig("secondarys_{}_{}_{}_stats_{}.pdf".format(
        mu_def.name,
        mu_def.mass,
        "ice",
        statistics
    ))

    # write execution time
    # end_time = time.time()

    # with open("time_{}_stats_{}.txt".format(mu.name, statistics), "w") as f:
    #     f.write("execution time: {}".format(end_time - start_time))
