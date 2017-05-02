
import pyPROPOSAL

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("Matplotlib not installed!")

try:
    import numpy as np
except ImportError:
    raise ImportError(
        "Numpy not installed! Needed to calculate the detector cylinder"
    )

ptype = pyPROPOSAL.ParticleType.MuMinus
mu = pyPROPOSAL.Particle(ptype)
med = pyPROPOSAL.Medium("ice")
E = pyPROPOSAL.EnergyCutSettings(500, 0.05)

photo = pyPROPOSAL.Photonuclear(mu, med, E)
photo.parametrization = 12
photo.enable_dEdx_interpolation()

dEdx = []
energy = []

for log_E in np.arange(0, 12, 0.2):
    E = mu.mass + np.power(10, log_E)
    mu.energy = E

    energy.append(E)
    dEdx.append(photo.calculate_dEdx() / E)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.semilogx(
    energy,
    dEdx,
    color='k',
    linestyle='-',
    label='ALLM 97'
)

ax.set_xlabel(r'$E$ / MeV')
ax.set_ylabel(r'')

ax.legend(loc='best')

fig.tight_layout()
fig.savefig('photo.pdf')
plt.show()
