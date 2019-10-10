# Documentation for Photon propagation #

For photon propagation PROPOSAL provides **Compton Scattering** and **Pair Production of electron positron pairs by photons** as possible interactions.

All other interactions should be disabled by choosing `None` as parametrization option when propagating photons.
[This configuration](config_photon.json) provides an example for valid options for photon propagation.

## Cross section parameters 

### Compton Scattering

Compton scattering describes the scattering of an ingoing photon by an atomic electron. These parametrizations are available and can be set with the `compton` parameter:
  - `"ComptonKleinNishina"` (SLAC Report number SLAC-R-730), atomic electrons are treated as free particles 

| Keyword                  | Type   | Default    | Description |
| -----------------------  | ------ | ---------- | ----------- |
| `compton_multiplier`     | Double | `1.0`      | Scales the compton parametrization |
| `compton`                | String | `None`     | Compton parametrization |

### Photo pair production

Photo pair production describes the production of an electron-positron pair by an incoming photon. These parametrizations are available and can be set with the `photopair` parameter:
  - `"PhotoPairTsai"` (Review of Modern Physics, Vol. 46, No. 4, October 1974)
  
Photo pair production interactions are always stochastic and always catastrophic: 
After an interaction, the initial photon is destroyed and the energy is divided between the electron and the positron.

To sample the angles of the created particles different models can be used. Those parametrizations are available and can be set via the `hotoangle` parameter:
  - `"PhotoAngleNoDeflection"` Both particles inherit the direction of the initial photon
  - `"PhotoAngleEGS"` The created particles are deflected by an angle of 1/energy with random azimuth
  - `"PhotoAngleTsaiIntegral"` The deflection angle is sampled for each interaction (Review of Modern Physics, Vol. 46, No. 4, October 1974)

| Keyword                  | Type   | Default                   | Description |
| -----------------------  | ------ | ------------------------- | ----------- |
| `photopair_multiplier`   | Double | `1.0`                     | Scales the photopairproduction parametrization |
| `photopair`              | String | `None`                    | Photopairproduction parametrization |
| `photoangle`             | String | `PhotoAngleNoDeflection`  | Parametrization of the deflection angle|
