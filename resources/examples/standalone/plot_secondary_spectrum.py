import sys
import os
import pyPROPOSAL as pp
import math
import time
import datetime

try:
    import matplotlib
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import AutoMinorLocator
    from mpl_toolkits.axes_grid1 import make_axes_locatable

except ImportError:
    print("Matplotlib not installed!")

try:
    import numpy as np
except ImportError:
    print("Numpy not installed!")

try:
    from sklearn.utils import check_random_state
except ImportError:
    print("SkLearn not installed!")


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
#

def power_law_sampler(gamma, xlow, xhig, n, random_state=None):
    r"""
    Sample n events from a power law with index gamma between xlow and xhig
    by using the analytic inversion method. The power law pdf is given by

    .. math::
       \mathrm{pdf}(\gamma) = x^{-\gamma} / \mathrm{norm}

    where norm ensures an area under curve of one. Positive spectral index
    gamma means a falling spectrum.

    Note: When :math:`\gamma=1` the integral is

    .. math::
       \int 1/x \mathrm{d}x = ln(x) + c

    This case is also handled.

    Sampling of power laws over multiple order of magnitude with the rejection
    method is VERY inefficient.

    Parameters
    ----------
    gamma : float
        Power law index.
    xlow, xhig : float
        Border of the pdf, needed for proper normalization.
    n : int
        Number of events to be sampled.
    random_state : seed, optional
        Turn seed into a np.random.RandomState instance. See
        `sklearn.utils.check_random_state`. (default: None)

    Returns
    -------
    sample : float array
        Array with the requested n numbers drawn distributed as a power law
        with the given parameters.
    """
    rndgen = check_random_state(random_state)
    # Get uniform random number which is put in the inverse function
    u = rndgen.uniform(size=int(n))

    if gamma == 1:
        return np.exp(u * np.log(xhig / xlow)) * xlow
    else:
        radicant = (u * (xhig**(1. - gamma) - xlow**(1. - gamma)) +
                    xlow**(1. - gamma))
        return radicant**(1. / (1. - gamma))
#


def propagate_muons():

    # start_time = time.time()

    mu_def = pp.particle.MuMinusDef.get()
    geometry = pp.geometry.Sphere(pp.Vector3D(), 1.e20, 0.0)
    ecut = 500
    vcut = 5e-2

    sector_def = pp.SectorDefinition()
    sector_def.cut_settings = pp.EnergyCutSettings(ecut, vcut)
    sector_def.medium = pp.medium.StandardRock(1.0)
    sector_def.geometry = geometry
    sector_def.scattering_model = pp.scattering.ScatteringModel.NoScattering
    sector_def.crosssection_defs.brems_def.lpm_effect = False
    sector_def.crosssection_defs.epair_def.lpm_effect = False
    # sector_def.crosssection_defs.photo_def.parametrization = pp.parametrization.photonuclear.PhotoParametrization.BezrukovBugaev
    # sector_def.do_stochastic_loss_weighting = True
    # sector_def.stochastic_loss_weighting = -0.1

    detector = geometry

    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = "~/.local/share/PROPOSAL/tables"

    prop = pp.Propagator(mu_def, [sector_def], detector, interpolation_def)

    statistics_log = 4
    statistics = int(10**statistics_log)
    propagation_length = 1e4 # cm
    E_min_log = 10.0
    E_max_log = 10.0
    spectral_index = 1.0
    pp.RandomGenerator.get().set_seed(1234)

    # muon_energies = np.logspace(E_min_log, E_max_log, statistics)
    # muon_energies = power_law_sampler(spectral_index, 10**E_min_log, 10**E_max_log, statistics)
    muon_energies = np.ones(statistics)*10**10

    epair_secondary_energy = []
    brems_secondary_energy = []
    ioniz_secondary_energy = []
    photo_secondary_energy = []

    progress = ProgressBar(statistics, pacman=True)
    progress.start()

    for mu_energy in muon_energies:
        progress.update()

        prop.particle.position = pp.Vector3D(0, 0, 0)
        prop.particle.direction = pp.Vector3D(0, 0, -1)
        prop.particle.propagated_distance = 0
        prop.particle.energy = mu_energy

        secondarys = prop.propagate(propagation_length)

        for sec in secondarys:
            log_sec_energy = math.log10(sec.energy)

            if sec.id == pp.particle.Data.Epair:
                epair_secondary_energy.append(log_sec_energy)
            if sec.id == pp.particle.Data.Brems:
                brems_secondary_energy.append(log_sec_energy)
            if sec.id == pp.particle.Data.DeltaE:
                ioniz_secondary_energy.append(log_sec_energy)
            if sec.id == pp.particle.Data.NuclInt:
                photo_secondary_energy.append(log_sec_energy)

    # =========================================================
    #   Write
    # =========================================================

    dir_prefix = ""
    np.savez(
        os.path.join(
            dir_prefix,
            'data_sec_dist_{}_{}_Emin_{}_Emax_{}'.format(
                prop.particle.particle_def.name,
                sector_def.medium.name.lower(),
                E_min_log,
                E_max_log,
                ecut,
                vcut)
        ),
        brems=brems_secondary_energy,
        epair=epair_secondary_energy,
        photo=photo_secondary_energy,
        ioniz=ioniz_secondary_energy,
        statistics=[statistics],
        E_min=[E_min_log],
        E_max=[E_max_log],
        spectral_index=[spectral_index],
        distance=[prop.particle.propagated_distance / 100],
        medium_name=[sector_def.medium.name.lower()],
        particle_name=[prop.particle.particle_def.name],
        ecut=[ecut],
        vcut=[vcut]
    )
#


def plot_secondary_spectrum():

    # =========================================================
    #   Plot
    # =========================================================

    # tex_preamble = [
    #     r"\usepackage{amsmath}",
    #     r"\usepackage[utf8]{inputenc}",
    #     r"\usepackage[T1]{fontenc}",
    # ]

    # font_size = 10

    # params = {
    #     'backend': 'pdf',
    #     'font.family': 'serif',
    #     'font.size': 12,
    #     'text.usetex': True,
    #     'text.latex.preamble': tex_preamble,
    #     'axes.labelsize': font_size,
    #     'legend.numpoints': 1,
    #     'legend.shadow': False,
    #     'legend.fontsize': font_size,
    #     'xtick.labelsize': font_size,
    #     'ytick.labelsize': font_size,
    #     'axes.unicode_minus': True
    # }

    # plt.rcParams.update(params)

    inch_to_cm = 2.54
    golden_ratio = 1.61803
    width = 29.7  # cm

    # # =========================================================
    # #   All hists together $10^{{{:.0g}}}$
    # # =========================================================

    npzfile = np.load("data_sec_dist_MuMinus_standardrock_Emin_10.0_Emax_10.0.npz")

    ioniz_secondary_energy = npzfile['ioniz']
    brems_secondary_energy = npzfile['brems']
    photo_secondary_energy = npzfile['photo']
    epair_secondary_energy = npzfile['epair']

    statistics = npzfile['statistics'][0]
    E_min_log = npzfile['E_min'][0]
    E_max_log = npzfile['E_max'][0]
    spectral_index = npzfile['spectral_index'][0]
    distance = npzfile['distance'][0]
    medium_name = npzfile['medium_name'][0]
    particle_name = npzfile['particle_name'][0]
    ecut = npzfile['ecut'][0]
    vcut = npzfile['vcut'][0]


    fig_all = plt.figure(
        figsize=(width / inch_to_cm, width / inch_to_cm / golden_ratio)
    )
    fig_all.suptitle(r"{:g} {} from $10^{{{:.2g}}}$ to $10^{{{:.2g}}}$ MeV from $E^{{-{:.2g}}}$ spectrum propagated {} m in {}".format(
        statistics,
        particle_name,
        E_min_log,
        E_max_log,
        spectral_index,
        distance,
        medium_name,
    ))

    ax_all = fig_all.add_subplot(111)
    ax_all.hist(
        [
            ioniz_secondary_energy,
            photo_secondary_energy,
            brems_secondary_energy,
            epair_secondary_energy,
            np.concatenate((
                ioniz_secondary_energy,
                brems_secondary_energy,
                photo_secondary_energy,
                epair_secondary_energy)
            )
        ],
        histtype='step',
        log=True,
        bins=100,
        label=['Ionization', 'Photonuclear', 'Bremsstrahlung', 'Pair Production', 'Sum']
    )
    # ax_all.set_ylim(ymin=0)
    minor_locator = AutoMinorLocator()
    ax_all.xaxis.set_minor_locator(minor_locator)
    ax_all.legend()
    ax_all.set_xlabel(r'energy loss / log($E$/MeV)')
    ax_all.set_ylabel(r'$N$')

    fig_all.savefig("all_{}_stats_{}_Emin_{}_Emax_{}_index_{}.pdf".format(
        medium_name,
        statistics,
        E_min_log,
        E_max_log,
        spectral_index
    ))
#


def plot_theory_curve():

    npzfile = np.load("data_sec_dist_MuMinus_standardrock_Emin_10.0_Emax_10.0.npz")

    ioniz_secondary_energy = npzfile['ioniz']
    brems_secondary_energy = npzfile['brems']
    photo_secondary_energy = npzfile['photo']
    epair_secondary_energy = npzfile['epair']

    all_secondary_energy = np.concatenate((
        ioniz_secondary_energy,
        brems_secondary_energy,
        photo_secondary_energy,
        epair_secondary_energy)
    )

    sum_hist = np.sum(all_secondary_energy)
    print(sum_hist)

    list_secondary_energies = [
        ioniz_secondary_energy,
        brems_secondary_energy,
        photo_secondary_energy,
        epair_secondary_energy,
        all_secondary_energy
    ]

    list_secondary_energies_label = [
        'Ionization',
        'Photonuclear',
        'Bremsstrahlung',
        'Pair Production',
        'Sum'
    ]

    statistics = npzfile['statistics'][0]
    E_min_log = npzfile['E_min'][0]
    E_max_log = npzfile['E_max'][0]
    spectral_index = npzfile['spectral_index'][0]
    distance = npzfile['distance'][0]
    medium_name = npzfile['medium_name'][0]
    particle_name = npzfile['particle_name'][0]
    ecut = npzfile['ecut'][0]
    vcut = npzfile['vcut'][0]

    particle_def = pp.particle.MuMinusDef.get()
    medium = pp.medium.StandardRock(1.0)
    energy_cuts = pp.EnergyCutSettings(500, -1)
    multiplier = 1.
    lpm = False
    shadow_effect = pp.parametrization.photonuclear.ShadowButkevichMikhailov()
    add_pertubative = True
    interpolation_def = pp.InterpolationDef()
    interpolation_def.path_to_tables = "~/.local/share/PROPOSAL/tables"

    ioniz = pp.parametrization.Ionization(
        particle_def=particle_def,
        medium=medium,
        energy_cuts=energy_cuts,
        multiplier=multiplier)

    epair = pp.parametrization.pairproduction.EpairProductionRhoInterpolant(
        particle_def=particle_def,
        medium=medium,
        energy_cuts=energy_cuts,
        multiplier=multiplier,
        lpm_effect=lpm,
        interpolation_def=interpolation_def)

    brems = pp.parametrization.bremsstrahlung.KelnerKokoulinPetrukhin(
        particle_def=particle_def,
        medium=medium,
        energy_cuts=energy_cuts,
        multiplier=multiplier,
        lpm_effect=lpm)

    photo = pp.parametrization.photonuclear.AbramowiczLevinLevyMaor97Interpolant(
        particle_def=particle_def,
        medium=medium,
        energy_cuts=energy_cuts,
        multiplier=multiplier,
        shadow_effect=shadow_effect,
        interpolation_def=interpolation_def)

    photo2 = pp.parametrization.photonuclear.BezrukovBugaev(
        particle_def=particle_def,
        medium=medium,
        energy_cuts=energy_cuts,
        multiplier=multiplier,
        add_pertubative=add_pertubative)

    losses_params = [ioniz, epair, brems, photo2]

    muon_energy = 10**10  # MeV

    inch_to_cm = 2.54
    golden_ratio = 1.61803
    width = 29.7  # cm

    num_bins = 100
    v_bins = np.linspace(np.log10(500), np.log10(muon_energy), num_bins)
    v_bins_log = np.logspace(np.log10(500./muon_energy), np.log10(1.), num_bins)

    fig = plt.figure(figsize=(width / inch_to_cm, width / inch_to_cm / golden_ratio))
    ax = fig.add_subplot(111)

    for idx, secondary_list in enumerate(list_secondary_energies):
        ax.hist(
            secondary_list,
            weights=np.ones(len(secondary_list))/sum_hist,
            histtype='step',
            log=True,
            bins=v_bins,
            label=list_secondary_energies_label[idx]
        )

    all_cross_sections = np.empty((len(losses_params), num_bins))

    for idx, param in enumerate(losses_params):
        all_cross_sections[idx] = np.array([param.differential_crosssection(muon_energy, v)*v for v in v_bins_log])

    sum_cross_sections = np.sum(all_cross_sections, axis=0)
    print(sum(all_cross_sections[np.isfinite(all_cross_sections)]))
    print(sum(sum_cross_sections[np.isfinite(sum_cross_sections)]))

    for cross_section in np.append(all_cross_sections, [sum_cross_sections], axis=0):
        ax.plot(
            v_bins,
            cross_section/sum(sum_cross_sections[np.isfinite(sum_cross_sections)]),
            drawstyle="steps-pre"
        )

    # ax.set_xscale("log")
    minor_locator = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minor_locator)
    ax.legend()
    ax.set_xlabel(r'energy loss / log($E$/MeV)')
    ax.set_ylabel(r'$N$')

    ax.set_ylim(ymin=1e-8)
    ax.set_yscale("log")
    ax.legend()
    fig.savefig("theory_curve.pdf")

if __name__ == "__main__":
    propagate_muons()
    plot_secondary_spectrum()
    # plot_theory_curve()
#
