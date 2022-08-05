# Documentation of default CrossSections

To describe particle interactions, parametrizations for each interaction type are needed.
PROPOSAL provides default parametrizations for each particle type.
They serve as a recommendation and provide all relevant interactions for a standard simulation.

When initializing a Propagator, if no CrossSections are specified in the [configuration file](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/docs/config_docu.md),
the default CrossSections will be used automatically.

In python, the default CrossSections can be obtained by using the command:

```
cross = pp.crosssection.make_std_crosssection(particle_def, medium, 
                                              energy_cut_settings, do_interpolate)
```

## List of default CrossSections

The default CrossSections for each particle type are defined as follows.
More information about each parametrization are given [here](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/docs/default_crosssections.md).

### Electrons:

- Bremsstrahlung: `ElectronScreening` (LPM effect disabled)
- Electron-Positron pair production: `ForElectronPositron` (LPM effect disabled)
- Ionization: `BergerSeltzerBhabha`
- Nuclear interactions: `AbramowiczLevinLevyMaor97` (with the shadowing parametrization `ButkevichMikheyev`)

### Positrons:

- Bremsstrahlung: `ElectronScreening` (LPM effect disabled)
- Electron-Positron pair production: `ForElectronPositron` (LPM effect disabled)
- Ionization: `BergerSeltzerMoller`
- Nuclear interactions: `AbramowiczLevinLevyMaor97` (with the shadowing parametrization `ButkevichMikheyev`)
- Annihilation: `Heitler`

### Muons:

- Bremsstrahlung: `KelnerKokoulinPetrukhin` (LPM effect disabled)
- Electron-Positron pair production: `KelnerKokoulinPetrukhin` (LPM effect disabled)
- Ionization: `BetheBlochRossi`
- Nuclear interactions: `AbramowiczLevinLevyMaor97` (with the shadowing parametrization `ButkevichMikheyev`)

### Tau lepton:

- Bremsstrahlung: `KelnerKokoulinPetrukhin` (LPM effect disabled)
- Electron-Positron pair production: `KelnerKokoulinPetrukhin` (LPM effect disabled)
- Ionization: `BetheBlochRossi`
- Nuclear interactions: `AbramowiczLevinLevyMaor97` (with the shadowing parametrization `ButkevichMikheyev`)

### Gamma:

- Electron-positron pair production by photons: `KochMotz`
- Compton scattering: `KleinNishina`
- Photoproduction: `Rhode`
- Photoeffect: `Sauter`