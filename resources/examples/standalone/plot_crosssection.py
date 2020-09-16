import proposal as pp
from plot_utility import plot_crosssection
import matplotlib.pyplot as plt

plt.style.use("Solarize_Light2")

param = pp.parametrization.bremsstrahlung.KelnerKokoulinPetrukhin(False)
# param = pp.parametrization.pairproduction.KelnerKokoulinPetrukhin(False)
particle = pp.particle.EMinusDef()
medium = pp.medium.Air()

fig = plt.figure(figsize=(15,5))
plot_crosssection(fig, param, particle, medium, energies=[1e3, 1e9])
plt.tight_layout()
plt.show()
