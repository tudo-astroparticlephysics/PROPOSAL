import pyPROPOSAL
import numpy as np

ptype = pyPROPOSAL.ParticleType.MuMinus
mu = pyPROPOSAL.Particle(ptype)
med = pyPROPOSAL.Medium("standard_rock")
E = pyPROPOSAL.EnergyCutSettings(-1, -1)

photo = pyPROPOSAL.Photonuclear(mu, med, E)
photo.parametrization = 12
photo.enable_dEdx_interpolation()

print(photo)

dEdx = list()
energy = list()

for log_E in np.arange(0, 12, 0.2):
    E = mu.mass + np.power(10, log_E)
    mu.energy = E

    energy.append(E)
    dEdx.append(photo.calculate_dEdx() / E)

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import pylab

    pylab.figure()
    pylab.title("Mu Lengths")
    pylab.loglog(energy, dEdx, ls="-")
    pylab.xlabel(r'$l_{\mu}(\rm{m})$')
    pylab.savefig('Photo.pdf')

except ImportError :
    print("pylab not installed.  no plots for you.")
