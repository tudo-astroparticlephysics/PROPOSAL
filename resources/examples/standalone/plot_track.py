
import math

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    raise ImportError("Matplotlib not installed!")

try:
    import numpy as np
except ImportError:
    raise ImportError(
        "Numpy not installed! Needed to calculate the detector cylinder"
    )

import pyPROPOSAL


def plot_3D_cylinder(
    ax,
    radius,
    height,
    elevation=0,
    resolution=100,
    color='r',
    alpha=0.1,
    x_center=0,
    y_center=0
):
    x = np.linspace(x_center-radius, x_center+radius, resolution)
    z = np.linspace(elevation, elevation+height, resolution)
    X, Z = np.meshgrid(x, z)

    Y = np.sqrt(radius**2 - (X - x_center)**2) + y_center

    ax.plot_wireframe(
        X,
        Y,
        Z,
        linewidth=1.0,
        color=color,
        alpha=alpha,
        rcount=1,
        ccount=10
    )
    ax.plot_wireframe(
        X,
        (2*y_center-Y),
        Z,
        linewidth=1.0,
        color=color,
        alpha=alpha,
        rcount=1,
        ccount=10
    )

    ax.w_xaxis.gridlines.set_lw(0.2)
    ax.w_yaxis.gridlines.set_lw(0.2)
    ax.w_zaxis.gridlines.set_lw(0.2)

    return ax


def sphere_to_cart(magnitude, phi, theta):
    z = magnitude * math.cos(theta)
    lx = magnitude * math.sin(theta)
    x = lx * math.cos(phi)
    y = lx * math.sin(phi)

    return x, y, z


def propagate():
    """ Propagte muon in ice threw a cylindric detector

    Returns:
        (Particle) Particle representing the start position
        (Geometry) Geometry of the detector
        (list)     List of secondarys particles represeint interactions
    """
    ptype = pyPROPOSAL.ParticleType.MuMinus
    mu = pyPROPOSAL.Particle(ptype)

    med = pyPROPOSAL.Medium("ice")

    cuts_inside = pyPROPOSAL.EnergyCutSettings(500, -1)
    cuts_infront = pyPROPOSAL.EnergyCutSettings(-1, 0.05)
    cuts_behind = pyPROPOSAL.EnergyCutSettings(-1, 0.05)

    prop = pyPROPOSAL.Propagator(
        med,
        cuts_inside,
        ptype,
        "../../resources/tables"
    )
    mu = prop.particle

    mu.energy = 9e6
    mu.Z = -1e5

    mu_start = pyPROPOSAL.Particle(mu)

    geo_detector = pyPROPOSAL.Geometry()
    geo_detector.init_cylinder(
        x0=0,
        y0=0,
        z0=0,
        radius=800,
        inner_radius=0,
        z=1600
    )

    geo_outside = pyPROPOSAL.Geometry()
    geo_outside.init_box(
        x0=0,
        y0=0,
        z0=0,
        x=500000,
        y=500000,
        z=500000
    )

    collection_detector = pyPROPOSAL.ProcessCollection(mu, med, cuts_inside)
    collection_detector.location = 1
    collection_detector.geometry = geo_detector
    collection_detector.enable_randomization = True

    collection_infront = pyPROPOSAL.ProcessCollection(mu, med, cuts_infront)
    collection_infront.location = 0
    collection_infront.geometry = geo_outside
    collection_infront.enable_randomization = True

    collection_behind = pyPROPOSAL.ProcessCollection(mu, med, cuts_behind)
    collection_behind.location = 2
    collection_behind.geometry = geo_outside
    collection_behind.enable_randomization = True

    prop.collections = [
        collection_infront,
        collection_detector,
        collection_behind
    ]

    prop.detector = geo_detector
    prop.apply_options()

    secondarys = prop.propagate()

    return mu_start, geo_detector, secondarys


def plot_track(mu_start, geo_detector, secondarys):
    """ Create 3D plot of the track.

    Args:
        mu_start (Particle): particle representing the start position
        geo_detector (Geometry): Geometry of the detector
        secondarys (list): List of secondarys particles representing the
                           interactions

    Returns: None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax = plot_3D_cylinder(
        ax,
        geo_detector.radius / 100,
        geo_detector.z / 100,
        elevation=-geo_detector.z / 2.0 / 100,
        resolution=10,
        color='k',
        alpha=0.5,
        x_center=geo_detector.x0 / 100,
        y_center=geo_detector.y0 / 100
    )

    ax.scatter(
        mu_start.X/100,
        mu_start.Y/100,
        mu_start.Z/100,
        c='r'
    )

    x, y, z = sphere_to_cart(200, mu_start.phi, mu_start.theta)

    ax.quiver(
        mu_start.X/100,
        mu_start.Y/100,
        mu_start.Z/100,
        x,
        y,
        z
    )

    for interaction in secondarys:
        if geo_detector.is_particle_inside(interaction):
            color = 'b'
            if interaction.type == pyPROPOSAL.ParticleType.Brems:
                color = 'b'
            elif interaction.type == pyPROPOSAL.ParticleType.EPair:
                color = 'r'
            elif interaction.type == pyPROPOSAL.ParticleType.NuclInt:
                color = 'm'
            elif interaction.type == pyPROPOSAL.ParticleType.DeltaE:
                color = 'g'
            else:
                color = 'k'

            ax.scatter(
                interaction.X/100,
                interaction.Y/100,
                interaction.Z/100,
                c=color,
                alpha=0.2,
                s=interaction.energy/400
            )

    ax.set_xlabel('x / m')
    ax.set_ylabel('y / m')
    ax.set_zlabel('z / m')

    plt.show()

    return fig

if __name__ == "__main__":

    fig = plot_track(*propagate())
    # fig.savefig('track.pdf')
