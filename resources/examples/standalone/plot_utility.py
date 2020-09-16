import numpy as np
import proposal as pp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def __transitions(energies):
    return np.geomspace(0.4, 1, len(energies))


def _plot_dsigmadv(ax, param, particle, target, energies, color="C1"):
    ax.set_yscale("log")
    if isinstance(target, (tuple, np.ndarray, list)):
        for i, t_i in enumerate(target):
            _plot_dsigmadv(ax, param, particle, t_i, energies, f"C{i}")
        return

    v = np.linspace(0, 1, 100)
    for energy, alpha in zip(energies, __transitions(energies)):
        sigma = []
        for v_i in v:
            sigma.append(
                param.differential_crosssection(particle, target, energy, v_i)
            )
        (line,) = ax.plot(v, sigma, alpha=alpha, color=color)
    line.set_label(target.name)
    ax.legend(loc="best")
    ax.set_xlabel("v")
    ax.set_ylabel(r"d$\sigma$/dv")


def _plot_dEdx(ax, param, particle, target):

    energies = np.geomspace(1e3, 1e12, 40)
    for cut in [
        pp.EnergyCutSettings(np.Infinity, 1, False),
        pp.EnergyCutSettings(500, 0.05, False),
        pp.EnergyCutSettings(np.Infinity, 0.05, False),
    ]:
        cross = pp.crosssection.make_crosssection(
            param, particle, target, cut, False
        )
        ax.plot(
            energies,
            cross.calculate_dEdx(energies),
            label=f"cut = ({cut.ecut}, {cut.vcut}) ",
        )
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlabel("Energy / MeV")
        ax.set_ylabel(r"dE/dx / [MeV cm^2/g]")
    ax.legend(loc="best")

def _plot_dNdx(ax, param, particle, target):
    for cut in [
            pp.EnergyCutSettings(500, 0.05, False),
            pp.EnergyCutSettings(np.Infinity, 0.05, False),
            ]:
        cross = pp.crosssection.make_crosssection(param, particle, target, cut, False)
        energies = np.geomspace(1e3, 1e12, 40)
        dNdx = []
        for e in energies:
            rates = cross.calculate_dNdx(e)
            dNdx.append(sum(rates.values()))
        ax.plot(
            energies,
            dNdx,
            label=f"cut = ({cut.ecut}, {cut.vcut})",
        )
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("Energy / MeV")
    ax.legend(loc="best")




def plot_crosssection(fig, param, particle, medium, energies=(1e3, 1e12)):

    if isinstance(param, (tuple, np.ndarray, list)):
        print(param, "is a list")
        for p_i in param:
            plot_crosssection(
                ax, param, p_i, particle, medium.components, energies
            )

    ax1 = fig.add_subplot(131)
    _plot_dsigmadv(ax1, param, particle, medium.components, energies)
    ax2 = fig.add_subplot(132)
    _plot_dEdx(ax2, param, particle, medium)
    ax3 = fig.add_subplot(133)
    _plot_dNdx(ax3, param, particle, medium)
