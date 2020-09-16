import proposal as pp
import numpy as np
muon = pp.particle.MuMinusDef()
medium = pp.medium.Ice()
brems = pp.parametrization.bremsstrahlung.KelnerKokoulinPetrukhin(False)
cut = pp.EnergyCutSettings(np.Infinity,1, False)
cross = pp.crosssection.make_crosssection(brems, muon, medium, cut, False)
cross.calculate_dNdx(1e6)
