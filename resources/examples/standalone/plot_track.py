
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

import pyPROPOSAL as pp


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

    geo_def_detector = pp.GeometryDefinition()
    geo_def_detector.shape = pp.Shape.Cylinder
    geo_def_detector.radius = 80000
    geo_def_detector.height = 160000

    geo_def_outside = pp.GeometryDefinition()
    geo_def_outside.shape = pp.Shape.Box
    geo_def_outside.depth = 50000000
    geo_def_outside.height = 50000000
    geo_def_outside.width = 50000000

    med_def = pp.MediumDefinition()
    med_def.type = pp.MediumType.Ice

    # Infront

    sec_def_infront = pp.SectorDefinition()
    sec_def_infront.medium_def = med_def
    sec_def_infront.geometry_def = geo_def_outside
    sec_def_infront.particle_location = pp.ParticleLocation.infront_detector

    sec_def_infront.scattering_model = pp.ScatteringModel.moliere

    sec_def_infront.e_cut = -1
    sec_def_infront.v_cut = 0.05

    # Inside

    sec_def_inside = pp.SectorDefinition()
    sec_def_inside.medium_def = med_def
    sec_def_inside.geometry_def = geo_def_outside
    sec_def_inside.particle_location = pp.ParticleLocation.inside_detector

    sec_def_inside.scattering_model = pp.ScatteringModel.moliere

    sec_def_inside.e_cut = 500
    sec_def_inside.v_cut = -1

    # Behind

    sec_def_behind = pp.SectorDefinition()
    sec_def_behind.medium_def = med_def
    sec_def_behind.geometry_def = geo_def_outside
    sec_def_behind.particle_location = pp.ParticleLocation.behind_detector

    sec_def_behind.scattering_model = pp.ScatteringModel.moliere

    sec_def_behind.e_cut = -1
    sec_def_behind.v_cut = 0.05

    interpolation_def = pp.InterpolationDef()
    # interpolation_def.path_to_tables = "~/.local/share/PROPOSAL/tables"

    prop = pp.Propagator(
        particle_def=pp.MuMinusDef.get(),
        sector_defs=[sec_def_infront, sec_def_inside, sec_def_behind],
        detector=pp.Cylinder(pp.Vector3D(), 0, 800, 1600),
        interpolation_def=interpolation_def
    )

    mu = prop.particle

    mu.energy = 9e6
    mu.direction = pp.Vector3D(0, 0, 1)

    pos = mu.position
    pos.set_cartesian_coordinates(0, 0, -1e5)
    mu.position = pos

    mu_start = pp.Particle(mu)

    secondarys = prop.propagate()

    return mu_start, prop.detector, secondarys


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
        geo_detector.height / 100,
        elevation=-geo_detector.height / 2.0 / 100,
        resolution=10,
        color='k',
        alpha=0.5,
        x_center=geo_detector.position.x / 100,
        y_center=geo_detector.position.y / 100
    )

    ax.scatter(
        mu_start.position.x/100,
        mu_start.position.y/100,
        mu_start.position.z/100,
        c='r'
    )

    x, y, z = sphere_to_cart(
        200,
        mu_start.direction.phi,
        mu_start.direction.theta
    )

    ax.quiver(
        mu_start.position.x/100,
        mu_start.position.y/100,
        mu_start.position.z/100,
        x,
        y,
        z
    )

    for interaction in secondarys:
        if geo_detector.is_inside(interaction.position, interaction.direction):
            color = 'b'
            if interaction.id == pp.Data.Brems:
                color = 'b'
            elif interaction.id == pp.Data.Epair:
                color = 'r'
            elif interaction.id == pp.Data.NuclInt:
                color = 'm'
            elif interaction.id == pp.Data.DeltaE:
                color = 'g'
            else:
                color = 'k'

            ax.scatter(
                interaction.position.x/100,
                interaction.position.y/100,
                interaction.position.z/100,
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
