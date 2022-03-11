# PROPOSAL configuration file documentation
PROPOSAL `Propagator` objects can be initialized in an easy way using a **json configuration file**, describing the propagation environment. This documentation provides an overview of the different options that can be used.

Examples for configuration files can be found in the examples folder.

The propagation environment of PROPOSAL is structured using `Sectors`. Each sector is defined by its `Geometry`, `Medium` and other individual properties. These properties will be described in this documentation.

The top-level of the json configuration file consists of two objects:

| Keyword          | Type   | Description |
| ---------------  | ------ | ----------- |
| `sectors`        | Array  | List of Sector objects. Each Sector object includes options that describe the Sector properties. |
| `global`         | Object | Defines global Sector settings. First of all, each Sector will use the settings defined in the corresponding Sector object. If the Sector object does not define a setting, PROPOSAL will look whether this setting is defined in this global object. If the option is not defined here, PROPOSAL will either use a default value or throw an error when defining this option if mandatory. |

#### Example:

```json5
{
	"global" : {
		// Describe global settings here
	},
	"sectors" : [
		{
			// Define properties of first Sector
		},
		{
			// Define properties of second Sector
		},
		// ...
	],
	"interpolation" : {
		// Define interpolation options here
	}
}
```

# Sector properties

All keywords with `-` as a default option are mandatory and need to be defined either for the Sector or in the global settings!

| Keyword            |  Type   | Default | Description |
| ------------------ | ------- | ------- | ----------- |
| `medium` | String | `-`  | Medium of the Sector |
| `cuts`   | Object | `-`  | Energy cut settings to be used in this Sector. |
| `geometries` | Array | `-` | List of geometry objects describing the geometry of the Sector. |
| `do_interpolation` | Boolean | `true`  | Defines if interpolation tables should be used for propagation. Note that not using interpolation tables will increase the runtime by several orders of magnitude! |
| `exact_time` | Boolean | `true`  | Defines if the elapsed time will be calculated exactly using the actual particle velocity or by using the approximation that all particles travel with the speed of light. |
| `scattering` | Object | No scattering  | Object to define multiple scattering and stochastic deflection behaviour. Per default, multiple scattering and stochastic deflection are disabled. |
| `density_distribution`   | Object | Homogeneous density distribution | Distribution of the mass density of the Sector. |
| `CrossSections` | Object | Standard cross sections | Cross sections that will be used in this Sector. |

## Medium

The following media are implemented and can be used with the `medium` setting:

* `Water` (always use the newest medium definition for Water)
* `WaterPDG2001`
* `WaterPDG2020`
* `Ice` (always use the newest medium definition for Ice)
* `IcePDG2001`
* `IcePDG2020`
* `Salt`
* `Standardrock`
* `Frejusrock`
* `Iron`
* `Hydrogen`
* `Lead`
* `Copper`
* `Uranium`
* `Air`
* `Paraffin`
* `AntaresWater`
* `CascadiabasinWater`
* `LiquidArgon`

## Energy cut settings

The **energy cut settings**, described by a `cuts` object, are used to differentiate between continuous losses (all energy losses with a relative energy transfer below the energy cut) and stochastic losses (all energy losses with a relative energy transfer above the energy cut). Here, the energy cut is combined from an absolute energy cut `e_cut` and a relative energy cut `v_cut`, where the minimum of both settings if used. With these options, the precision and performance of the propagation process can be adjusted. 

Furthermore, **continuous randomization** can be enabled to add random fluctuations to continuous energy losses. This can be useful when using rather high energy cuts to still consider fluctuations in energy losses.

| Keyword     |  Type            | Default | Description |
| ----------- | ---------------- | ------- | ----------- |
| `e_cut`     | Number or String | `-`     | Absolute energy cut in MeV. Pass `"inf"` or `"infinity"` to use an infinite energy cut. |
| `v_cut`     | Number           | `-`     | Relative energy cut. Pass `1` to use an infinite energy cut. |
| `cont_rand` | Boolean          | `-`     | Defines whether to use continuous randomization. |

#### Example:

```json5
"cuts" :
{
	"e_cut" : "inf",
	"v_cut" : 0.05,
	"cont_rand" : true
}
```

## Geometries

Each sector consists of at least one **geometry**. Three different geometry types can be used and need to be specified with the keyword `shape`, their `origin` as well as specific keywords depending on the geometry type. 

| Keyword     |  Type  | Default | Description                                                |
| ----------- | ------ | ------- | ---------------------------------------------------------- |
| `shape`     | String | `-`     | Type of geometry, i.e. `sphere`, `box` or `cylinder`.      |
| `origin`    | Array  | `-`     | Position of the center of the geometry, consisting out of three numbers describing the coordinates in cm. |
| `hierarchy` | Number | `0`     | Hierarchy of the geometry. For overlapping geometries, the geometry with the highest hierarchy will be used first. For overlapping Sector geometries with equal hierarchy, PROPOSAL will prefer the Sector that has been defined first in the json file. |


### Sphere

| Keyword        |  Type  | Default | Description                       |
| -------------- | ------ | ------- | --------------------------------- |
| `outer_radius` | Number | `-`     | Outer radius of the sphere in cm. |
| `inner_radius` | Number | `0`     | Inner radius of the sphere in cm. |

### Box

| Keyword  |  Type  | Default | Description                         |
| -------- | ------ | ------- | ----------------------------------- |
| `length` | Number | `-`     | Length of box in x-direction in cm. |
| `width`  | Number | `-`     | Width of box in y-direction in cm.  |
| `height` | Number | `-`     | Height of box in z-direction in cm. |

### Cylinder

| Keyword        |  Type  | Default | Description                               |
| -------------- | ------ | ------- | ----------------------------------------- |
| `height`       | Number | `-`     | Height of cylinder in z-direction in cm.  |
| `outer_radius` | Number | `-`     | Outer radius of the cylinder in cm.       |
| `inner_radius` | Number | `0`     | Inner radius of the cylinder in cm.       |

#### Example

```json5
"geometries" : [
	{
		"hierarchy" : 0,
		"shape" : "sphere",
		"origin" : [0, 0, 0],
		"inner_radius" : 0,
		"outer_radius" : 637.1e6
	},
	{
		// ...
	}
],
```

## Scattering
Particles in PROPOSAL may be deflected both during continuous propagation steps (multiple scattering) and in stochastic interactions (stochastic deflection).
The behaviour for both processes can be set individually using the `scattering` object.

| Keyword                            |  Type  | Default        | Description                                                |
| ---------------------------------- | ------ | -------------- | ---------------------------------------------------------- |
| `multiple_scattering`              | String | `NoScattering` | Parametrization used to describe multiple scattering interactions. |
| `multiple_scattering_multiplier`   | Number | `1.0`          | Factor to scale the multiple scattering effects by a constant factor. |
| `stochastic_deflection`            | Array  | `[]`           | List of parametrizations to describe stochastic deflection effects for different interaction types, or list of interaction types (indicated with the prefix `default`) where PROPOSAL will use default parametrizations to describe stochastic deflections. |
| `stochastic_deflection_multiplier` | Object | `1.0` for each interaction type | Object with multipliers to scale stochastic deflection effects for each interaction type. Names of keys must correspond to the interaction types where the multiplier should be applied. |

If no scattering object is defined, neither multiple scattering nor stochastic deflections will be considered during propagation.

### Multiple scattering
PROPOSAL provides different multiple scattering models, which can be set with the `multiple_scattering` parameter, namely:

* `NoScattering`: No multiple scattering effects.
* `HighlandIntegral`: Gaussian approximation of Molière theory, derived by [Highland](https://doi.org/10.1016/0029-554X(75)90743-0) and corrected by [Lynch/Dahl](https://doi.org/10.1016/0168-583X(91)95671-Y). 
* `Highland`: Same as `HighlandIntegral`, but with the approximation that the particle energy is constant during a propagation step. *This parametrization should only be used for small step sizes where this approximation is valid.*
* `Moliere`: Scattering based on [Molière's theory](https://zfn.mpdl.mpg.de/data/Reihe_A/3/ZNA-1948-3a-0078.pdf). *Significantly slower compared to other scattering models.*

Multiple scattering effects can be scaled with a constant factor using the option `multiple_scattering_multiplier`.
This scales the sampled deflections, which are described by their displacement in cartesian coordinates, with a constant factor for each component. 

### Stochastic deflection
PROPOSAL provides different parametrizations for stochastic deflection interactions, where an individual parametrization must be chosen for each interaction type.
These parametrizations can be set with the option `stochastic_deflection`. 
This array is expected to either include a list of parametrization names or a list of interaction types, indicated with the prefix `default`.
In the latter case, PROPOSAL will use an appropriate default parametrization for this interaction type.
For interaction types that are not covered in this Array, stochastic deflection will not be considered. 

Stochastic deflections may be scaled using the option `stochastic_deflection_multiplier`.
This object consists of multipliers, where each key should correspond to the interaction type and eac value to the scaling factor for this interaction type.
These multipliers will multiply a constant factor to the stochastic deflection, described by a radial deflection angle.
Interaction types that are not covered in this object will not be scaled (i.e. factor `1.0`).

#### Examples

Scattering object where we only want to consider multiple scattering, described by the `HighlandIntegral` parametrization:

```json5
"scattering":
{
	"multiple_scattering": "HighlandIntegral"
}
```

Scattering object where we want to consider, additional to multiple scattering, stochastic deflections.
Here, we want PROPOSAL to choose default parametrizations for the stochastic deflection:

```json5
"scattering":
{
	"multiple_scattering": "HighlandIntegral", 
	"stochastic_deflection": ["DefaultIoniz", "DefaultBrems", "DefaultEpair", "DefaultPhotonuclear"]
}
```

## Density distribution

Per default, PROPOSAL assumes a homogeneous density distribution with the standard medium density.
Additionally, inhomogeneous environments are supported using different density distribution models that can be specified with the keyword `type` as well a type-specific keywords.

| Keyword        |  Type  | Default | Description                                                |
| -------------- | ------ | ------- | ---------------------------------------------------------- |
| `type`         | String | `-`     | Type of density distribution, i.e. `homogeneous`, `exponential`, `polynomial` or `spline`.      |
| `mass_density` | Number | `-`     | Base mass density in g/cm^3. |

For all inhomogeneous density distributions, we need to define an axis along which the density distribution will be evaluated.

| Keyword     |  Type  | Default | Description                                              |
| ----------- | ------ | ------- | -------------------------------------------------------- |
| `axis_type` | String | `-`     | Axis type, either `radial` or `cartesian`.               |
| `fp0`       | Array  | `-`     | Origin of the axis (coordinates in cm).                                      |
| `fAxis`     | Array  | `-`     | Necessary for the cartesian axis, defines the direction of the axis. |

The `radial` axis starts at the origin `fp0` and the depth <img src="https://render.githubusercontent.com/render/math?math=d(\vec{x})"> will be calculated as the distance from this origin, i.e.

> <img src="https://render.githubusercontent.com/render/math?math=d(\vec{x}) = |\vec{x} -\vec{f_{p_0}}|">.

For the `cartesian` axis, the axis starts at the origin `fp0` and has a direction `fAxis`. The depth will be calculated with respect to this axis, i.e.

> <img src="https://render.githubusercontent.com/render/math?math=d(\vec{x}) = \vec{f_\mathrm{axis}} \cdot (\vec{x} -\vec{f_{p_0}})">.

### Homogeneous

Homogeneous density distribution with the density `mass_density`. No additional keywords required.

### Exponential

Exponential density distribution along an axis of the form 

> <img src="https://render.githubusercontent.com/render/math?math=\rho(x) = \rho_0 \cdot \exp{\left(\frac{d(x) - d_0}{\sigma}\right)}">

with the scaling parameter `sigma`, the shifting parameter `d_0` in cm, the base `mass_density` <img src="https://render.githubusercontent.com/render/math?math=\rho_0"> and the depth d relative to the axis.

| Keyword |  Type  | Default | Description                |
| ------- | ------ | ------- | -------------------------- |
| `sigma` | Number | `1.0`   | Scaling parameter          |
| `d0`    | Number | `0.0`   | Shifting parameter (in cm) |

#### Example

Model the earth's atmosphere by creating an exponential density distribution with an radial axis.
Here, we assume that `[0, 0, 0]` is the position of the earth's core and that the earth is a perfect sphere with a radius of `6.3781e8` cm:
 
```json5
"density_distribution":
{
	"type": "exponential",
	"mass_density" : 1.225e-3,
	"axis_type" : "radial",
	"fp0" : [0, 0, 0],
	"sigma" : -10.4e5,
	"d0" : 6.3781e8,
}
```

### Polynomial

Density distribution described by a polynomial of the form 

> <img src="https://render.githubusercontent.com/render/math?math=\rho(x) = \rho_0 \cdot \sum_{k=0}^{n} a_k d(x)^k">

with the base `mass_density` <img src="https://render.githubusercontent.com/render/math?math=\rho_0">, the `coefficients` <img src="https://render.githubusercontent.com/render/math?math=a_k"> and the depth d relative to the axis.

| Keyword        |  Type | Default | Description       |
| -------------- | ----- | ------- | ----------------- |
| `coefficients` | Array | `-`     | Array of polynomial coefficients, starting with the lowest order coefficient |

### Spline

Density distribution described by splines along an axis, using either `linear` or `cubic` spline interpolation.

| Keyword       |  Type  | Default | Description       |
| ------------- | ------ | ------- | ----------------- |
| `spline_type` | String | `-`     | Type of spline, either `linear` or `cubic`. |
| `x`           | Array  | `-`     | Coordinates (along the axis) to evaluate the density correction. |
| `y`           | Array  | `-`     | Density distribution correction evaluated at corresponding `x` coordinates. |

## Cross sections

PROPOSAL provides the option to consider additional interaction types as well as to use different physical parametrizations of interactions. 
If no `CrossSections` object is included in the json file, PROPOSAL chooses a set of interaction types and parametrizations appropriate for the particle. 
These default cross sections for each particle type are listed [here](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/docs/default_crosssections.md).

To use alternative parametrizations or to enable additional interaction types, the `CrossSections` object needs to contain one object with the appropriate name for each interaction type (e.g. `annihilation`, `brems`, etc.). If there is no object for an specific interaction type, this interaction type will be disabled.

Each of the defined objects has to contain the keyword `parametrization`, specifying the parametrization that should be used, as well as additional, type-specific keywords.
For each interaction type, a `multiplier` can be defined which scales the total cross section by a constant coefficient.

#### Example

Example where the interactions bremsstrahlung, electron-positron pair production, ionization and nuclear interactions are enabled. In contrast to the default parametrizations for muons and taus, the LPM effect for electron-positron pair production and bremsstrahlung will be enabled.

```json
"CrossSections" : {
	"brems": {
		"parametrization": "KelnerKokoulinPetrukhin",
		"lpm": true
	},
	"epair": {
		"parametrization": "KelnerKokoulinPetrukhin",
		"lpm": true
	},
	"ioniz": {
		"parametrization": "BetheBlochRossi"
	},
	"photo": {
		"parametrization": "AbramowiczLevinLevyMaor97",
		"shadow": "ButkevichMikheyev"
	}
}
```

*Note that all interaction types that are not listed in the `CrossSections` object will be disabled.*

### Annihilation: `annihilation`

| Keyword           |  Type  | Default | Description                |
| ----------------- | ------ | ------- | -------------------------- |
| `parametrization` | String | `-`     | Parametrization to be used |
| `multiplier`      | Number | `1.0`   | Multiplier to scale the cross section with a coefficient |

Annihilation of an ingoing positron with an atomic electron. Available `annihilation` parametrizations are:

* `Heitler`: Heitler formula

### Bremsstrahlung: `brems`

| Keyword           |  Type   | Default | Description                |
| ----------------- | ------- | ------- | -------------------------- |
| `parametrization` | String  | `-`     | Parametrization to be used |
| `multiplier`      | Number  | `1.0`   | Multiplier to scale the cross section with a coefficient |
| `lpm`             | Boolean | `true`  | Enabling or disabling the reduction of the cross section at high energies by the Landau-Pomeranchuk-Migdal effect.  |

Note that the LPM correction currently only supports the correct description of homogeneous media.
If the LPM effect is enabled, it will be assumed that the whole medium has the base `mass_density` <img src="https://render.githubusercontent.com/render/math?math=\rho_0">.

Available `brems` parametrizations are:

* `KelnerKokoulinPetrukhin`: ([Preprint MEPhI (1995) no. 024-95](http://cds.cern.ch/record/288828)) and (Phys. Atom. Nucl. 62 (1999), 272)
* `AndreevBezrukovBugaev`: (Phys. Atom. Nucl. 57 (1994), 2066)
* `PetrukhinShestakov`: (Canad. J. Phys. 46 (1968), 377)
* `CompleteScreening`: ([Rev. Mod. Phys. 46 (1974), 815](https://doi.org/10.1103/RevModPhys.46.815))
* `SandrockSoedingreksoRhode`: ([arXiv:1807.08475](https://arxiv.org/abs/1807.08475))
* `ElectronScreening`: For electrons ans positrons with empirical low-energy corrections (SLAC-R-730)

### Compton: `compton`

| Keyword           |  Type   | Default | Description                |
| ----------------- | ------- | ------- | -------------------------- |
| `parametrization` | String  | `-`     | Parametrization to be used |
| `multiplier`      | Number  | `1.0`   | Multiplier to scale the cross section with a coefficient |

Compton scattering of photons on free electrons. Available `compton` parametrizations are:

* `KleinNishina`: Klein-Nishina formula (SLAC Report number SLAC-R-730)

### Electron-positron pair production: `epair`

| Keyword           |  Type   | Default | Description                |
| ----------------- | ------- | ------- | -------------------------- |
| `parametrization` | String  | `-`     | Parametrization to be used |
| `multiplier`      | Number  | `1.0`   | Multiplier to scale the cross section with a coefficient |
| `lpm`             | Boolean | `true`  | Enabling or disabling the reduction of the cross section at high energies by the Landau-Pomeranchuk-Migdal effect.  |

Creation of electron-positron pairs by ingoing leptons.
Note that the LPM correction currently only supports the correct description of homogeneous media.
If the LPM effect is enabled, it will be assumed that the whole medium has the base `mass_density` <img src="https://render.githubusercontent.com/render/math?math=\rho_0">.

Available `epair` parametrizations are:

* `KelnerKokoulinPetrukhin`: (Proc. 12th ICCR (1971), 2436) with corrections for the interaction with atomic electrons (Phys. Atom. Nucl. 61 (1998), 448)
* `SandrockSoedingreksoRhode`: ([arXiv:1807.08475](https://arxiv.org/abs/1807.08475))
* `ForElectronPositron`: Adaption of (Phys. Atom. Nucl. Vol. 63, No.9 (2000), pp. 1603-1611, DOI: 10.1134/1.1312894) for electrons and positrons

### Ionization: `ioniz`

| Keyword           |  Type   | Default | Description                |
| ----------------- | ------- | ------- | -------------------------- |
| `parametrization` | String  | `-`     | Parametrization to be used |
| `multiplier`      | Number  | `1.0`   | Multiplier to scale the cross section with a coefficient |

Available `ioniz` parametrizations are:

* `BetheBlochRossi`: Ionization described by Bethe-Bloch formula with corrections for muons and taus by (B. B. Rossi. Prentice-Hall, Inc., Englewood Cliffs, NJ, 1952)
* `BergerSeltzerBhabha`: Ionization for positrons based on Berger-Seltzer formula 
* `BergerSeltzerMoller`: Ionization for electrons based on Berger-Seltzer forumla

### Muon pair production: `mupair`

| Keyword           |  Type   | Default | Description                |
| ----------------- | ------- | ------- | -------------------------- |
| `parametrization` | String  | `-`     | Parametrization to be used |
| `multiplier`      | Number  | `1.0`   | Multiplier to scale the cross section with a coefficient |

Creation of muon-antimuon pairs by ingoing leptons. Available `mupair` parametrizations:

* `KelnerKokoulinPetrukhin`:  (Physics of Atomic Nuclei, Vol. 63, No. 9, 2000, pp. 1603-1611)

### Nuclear interactions: `photo`

| Keyword           |  Type   | Default              | Description                |
| ----------------- | ------- | -------------------- | -------------------------- |
| `parametrization` | String  | `-`                  | Parametrization to be used |
| `multiplier`      | Number  | `1.0`   | Multiplier to scale the cross section with a coefficient |
| `hard_component`  | Boolean | `true`               | Enabling or disabling the hard component for parametrizations using the real photon approximation |
| `shadow`          | String  | `ButkevichMikheyev`  | Parametrization of the shadowing effect, relevant for parametrizations using the momentum integration approach |

There are two classes of `photo` parametrizations for nuclear interactions, using two different approaches:

The first approach is to use the approximation where a **real photon** scatters inelastically with a nucleus with an additional factor to make the photon virtual. Parametrizations using this approach are:

* `Zeus`: ([Eur. Phys. J. C 7 (1999), 609](https://doi.org/10.1007/s100529901084))
* `BezrukovBugaev`: (Sov. J. Nucl. Phys. 33 (1981), 635) with corrections for particles with higher mass like taus ([Phys. Rev. D67 (2003), 034027](https://doi.org/10.1103/PhysRevD.67.034027))
* `Kokoulin`
* `Rhode`

For these parametrizations, the hard component can be considered using the `hard_component` parameter.

The second approach is to use data from electron-positron scattering experiments and extrapolate to low momenta of the virtual photon transferred to the nucleus. Then, the cross section is **integrated over the photon momentum**. Parametrizations using this approach are:

* `AbramowiczLevinLevyMaor91`: ([Phys. Let. B269 (1991), 465](https://doi.org/10.1016/0370-2693(91)90202-2))
* `AbramowiczLevinLevyMaor97`: ([arXiv::hep-ph/9712415](https://arxiv.org/abs/hep-ph/9712415))
* `ButkevichMikheyev`: ([JETP 95 (2002), 11](https://doi.org/10.1134/1.1499897))
* `RenoSarcevicSu`: ([Astrop. Phys. 24 (2005), 107](https://doi.org/10.1016/j.astropartphys.2005.06.002)), uses the parametrization of ALLM97, but with corrections for spin 0 particles. Should only be used for particles with spin 0 (e.g. sTaus).

For there parametrizations, the parametrization of the shadowing effect can be chosen using the keyword `shadow`. Available `shadow` parametrizations are:

* `DuttaRenoSarcevicSeckel`: ([Phys. Rev. D63 (2001), 094020](https://doi.org/10.1103/PhysRevD.63.094020))
* `ButkevichMikheyev`: ([JETP 95 (2002), 11](https://doi.org/10.1134/1.1499897))

### Electron-positron production by photons: `photopair`

| Keyword           |  Type   | Default | Description                |
| ----------------- | ------- | ------- | -------------------------- |
| `parametrization` | String  | `-`     | Parametrization to be used |
| `multiplier`      | Number  | `1.0`   | Multiplier to scale the cross section with a coefficient |

Creation of an electron-positron pair by an ingoing photon. Available `photopair` parametrizations are:

* `Tsai`: (Review of Modern Physics, Vol. 46, No. 4, October 1974)
* `KochMotz`: (Review of Modern Physics, Vol. 61, 1959)

### Muon-antimuon production by photons: `photomupair`

| Keyword           |  Type   | Default | Description                |
| ----------------- | ------- | ------- | -------------------------- |
| `parametrization` | String  | `-`     | Parametrization to be used |
| `multiplier`      | Number  | `1.0`   | Multiplier to scale the cross section with a coefficient |

Creation of an muon-antimuon pair by an ingoing photon. Available `photomupair` parametrizations are:

* `BurkhardtKelnerKokoulin`: ([CERN-SL-2002-016-AP](https://cds.cern.ch/record/558831/files/sl-2002-016.pdf)) with atomic form factors from (Phys. Atom. Nucl. 60 (1997), 576)

### Photonuclear intractions by photons: `photoproduction`

| Keyword           |  Type   | Default | Description                |
| ----------------- | ------- | ------- | -------------------------- |
| `parametrization` | String  | `-`     | Parametrization to be used |
| `multiplier`      | Number  | `1.0`   | Multiplier to scale the cross section with a coefficient |

Photonuclear interaction of a photon with an atomic nucleus. Available `photoproduction` parametrizations are:

* `Zeus`: ([Eur. Phys. J. C 7 (1999), 609](https://doi.org/10.1007/s100529901084))
* `BezrukovBugaev`: (Sov. J. Nucl. Phys. 33 (1981), 635) with corrections for particles with higher mass like taus ([Phys. Rev. D67 (2003), 034027](https://doi.org/10.1103/PhysRevD.67.034027))
* `Caldwell`: (Phys. Rev. Let. 42 (1979), 553)
* `Kokoulin`
* `Rhode`

### Weak interaction: `weak`

| Keyword           |  Type   | Default | Description                |
| ----------------- | ------- | ------- | -------------------------- |
| `parametrization` | String  | `-`     | Parametrization to be used |
| `multiplier`      | Number  | `1.0`   | Multiplier to scale the cross section with a coefficient |

Weak interaction of an ingoing charged lepton. Available `weak` parametrizations are:

* `CooperSarkarMertsch`: ([arXiv:1106.3723](https://arxiv.org/abs/1106.3723v1))

# Global settings

Properties that should be set for all or more than one sector may be set in a `global` object at the top-level of the json configuration file.

Objects or keywords that can be defined in this global object are:

* `Medium`
* `CrossSections`
* `cuts`
* `exact_time`
* `do_interpolation`
* `scattering`

Note that these options can and will still be overwritten by options in the individual sector objects: PROPOSAL will first look if an object or keyword is defined the in sector object in the `sectors` list. Only if an option is undefined here, PROPOSAL uses the definition in the `global` setting sections.
If an option is undefined here as well, PROPOSAL will either use an appropriate default value or throw an exception if the option is mandatory.

#### Example

In this example, be define a sector which describes an earth made out of ice as well as a sector with an air atmosphere which surrounds the earth.
Since the EnergyCut settings `e_cut = 500`, `v_cut = 1`, `cont_rand = false` are defined for the air sector, these EnergyCuts will be used.
Since there are no `cuts` defined in the ice sector, PROPOSAL looks at the `global` section, where the EnergyCut settings `e_cut = inf`, `v_cut = 0.05`, `cont_rand = true` are defined and will therefore be used:

```json5
{
	"global":
	{
		"cuts":
		{
			"e_cut": "inf",
			"v_cut": 0.05,
			"cont_rand": true
		}
	},
	"sectors": [
		{
			"medium": "ice",
			"geometries": [
				{
					"hierarchy": 0,
					"shape": "sphere",
					"origin": [0, 0, 0],
					"outer_radius": 6.3781e8
				}
			]
		},
		{
			"medium": "air",
			"geometries": [
				{
					"hierarchy": 0,
					"shape": "sphere",
					"origin": [0, 0, 0],
					"inner_radius": 6.3781e8,
					"inner_radius": 1e20,
				}
			],
			"cuts":
			{
				"e_cut": "500",
				"v_cut": 1,
				"cont_rand": false
			}
		}
	]
}
```