import pyPROPOSAL
import numpy as np

ptype = pyPROPOSAL.ParticleType.MuMinus
mu = pyPROPOSAL.Particle(ptype)

med = pyPROPOSAL.Medium("ice")

cuts_inside = pyPROPOSAL.EnergyCutSettings(500, -1)
cuts_infront = pyPROPOSAL.EnergyCutSettings(-1, 0.05)
cuts_behind = pyPROPOSAL.EnergyCutSettings(-1, 0.05)

prop = pyPROPOSAL.Propagator(med, cuts_inside, ptype, "../../resources/tables")
mu = prop.particle

mu.energy = 9e6
mu.Z = -1e5

mu_start_energy = mu.energy
mu_start_pos = (mu.X, mu.Y, mu.Z)
mu_start_phi = mu.phi
mu_start_theta = mu.theta


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

# print(prop.collections)

secondarys = prop.propagate()

# # =========================================================
# # 	Plot
# # =========================================================
#
# try:
#     import matplotlib.pyplot as plt
#     from mpl_toolkits.mplot3d import Axes3D
# except ImportError:
#     print("matplotlib not installed.  no plots for you.")
#
# try:
#     import numpy as np
# except ImportError:
#     print("numpy not installed.  no plots for you.")
#
# tex_preamble = [
#     r"\usepackage{amsmath}",
#     r"\usepackage[utf8]{inputenc}",
#     r"\usepackage[T1]{fontenc}",
#     r"\usepackage{siunitx}",
# ]
#
# font_size = 10
#
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
#
# plt.rcParams.update(params)
#
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# # Cylinder
# x = np.linspace(-1, 1, 100)
# z = np.linspace(-2, 2, 100)
# Xc, Zc = np.meshgrid(x, z)
# Yc = np.sqrt(1-Xc**2)
#
# # Draw parameters
# rstride = 20
# cstride = 10
# ax.plot_surface(Xc, Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)
# ax.plot_surface(Xc, -Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)
#
# ax.set_xlabel("x / m")
# ax.set_ylabel("y / m")
# ax.set_zlabel("z / m")
# plt.show()
#
# # fig.savefig("track.pdf")

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle


def plot_3D_cylinder(
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

    # ax.plot_surface(X, Y, Z, linewidth=0, color=color, alpha=alpha)
    # ax.plot_surface(X, (2*y_center-Y), Z, linewidth=0, color=color, alpha=alpha)
    ax.plot_wireframe(
        X,
        Y,
        Z,
        linewidth=1.0,
        color=color,
        alpha=alpha,
        rcount=1,
        ccount=1
    )
    ax.plot_wireframe(
        X,
        (2*y_center-Y),
        Z,
        linewidth=1.0,
        color=color,
        alpha=alpha,
        rcount=1,
        ccount=1
    )

    ax.w_xaxis.gridlines.set_lw(0.2)
    ax.w_yaxis.gridlines.set_lw(0.2)
    ax.w_zaxis.gridlines.set_lw(0.2)
    # floor = Circle((x_center, y_center), radius, color=color, alpha=alpha)
    # ax.add_patch(floor)
    # art3d.pathpatch_2d_to_3d(floor, z=elevation, zdir="z")
    #
    # ceiling = Circle((x_center, y_center), radius, color=color, alpha=alpha)
    # ax.add_patch(ceiling)
    # art3d.pathpatch_2d_to_3d(ceiling, z=elevation+height, zdir="z")

    return ax


def sphere_to_cart(magnitude, phi, theta):
    z = magnitude * math.cos(theta)
    lx = magnitude * math.sin(theta)
    x = lx * math.cos(phi)
    y = lx * math.sin(phi)

    return x, y, z


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax = plot_3D_cylinder(
    geo_detector.radius / 100,
    geo_detector.z / 100,
    elevation=-geo_detector.z / 2.0 / 100,
    resolution=10,
    color='k',
    alpha=0.5,
    x_center=geo_detector.x0 / 100,
    y_center=geo_detector.y0 / 100
)

mu_end = prop.particle
ax.scatter(mu_start_pos[0]/100, mu_start_pos[1]/100, mu_start_pos[2]/100, c='r')

x, y, z = sphere_to_cart(200, mu_start_phi, mu_start_theta)
ax.quiver(mu_start_pos[0]/100, mu_start_pos[1]/100, mu_start_pos[2]/100, x, y, z)
# ax.scatter(mu_end.X/100, mu_end.Y/100, mu_end.Z/100, c='r')

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

ax.set_xlabel('x-axis')
ax.set_ylabel('y-axis')
ax.set_zlabel('z-axis')

plt.show()
# fig.savefig('track.pdf')
