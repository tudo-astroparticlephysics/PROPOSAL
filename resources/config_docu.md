
# Documentation of the Configuration file #

In the configuration file you can specify all parameters for your propagation.
It's written in a .json file. That's why we have this (nice) Docu. ;)

The configuration file is structured into three parts:
- In the "global" part general parameters are defined (path for interpolation tables, seed, etc.)
- In the "sectors" part the different media with their geometry are defined. The energy-cuts, which means how accurate the propagation process should be calculated, can be set for each sector, which then overrides the global configurations.
- In the "detector" part the geometry of the detector is defined.

An example for a configuration file can be found [here](config_docu.md).

## The `global` configurations ##

### General parameters ###

There is the option to set a `seed`, if the internal random number generator is used.
The seed has no effect, if you use an external random number generator.

If the Output of the Secondaries should not only include the stochastic energy losses (and the particles produced in an interaction or decay), but also the continuous energy losses, this can be set by the `do_continuous_energy_loss_output` parameter.

When the Output should just contain the secondaries (energy losses or particles produced in an interaction or decay), that occurred inside the detector volume, and not the ones outside of the detector, this can be set with the `only_loss_inside_detector` parameter.

| Keyword                            | Type    | Default   | Description |
| ---------------------------------- | ------- | --------- | ----------- |
| `seed`                             | Integer | `0`       | seed for the internal random number generator|
| `do_continuous_energy_loss_output` | Bool    | `False`   | Decides whether continuous losses should be emitted in the Output of Secondaries|
| `only_loss_inside_detector`        | Bool    | `False`   | Decides whether only secondaries created inside the detector should be included in the Output of Secondaries|

### Interpolation parameters ###
The `interpolation` parameter is an own json-object.
This object can contain multiple parameters dealing with the interpolation tables and described in the following section.

By default, PROPOSAL uses interpolation tables to store the calculated results.
It is also possible to do the calculations with integrations by setting the parameter `do_interpolation` to False
Be aware that this increases the amount of calculation time needed by around four orders of magnitude ! This should just be used for tests and comparisons. 
The interpolation tables have an accuracy of at least 1.e-5 for most energies. 
Note: The reason why integration is that much slower is that the interpolation used for propagation is near the surface, using already calculated numbers from 'lower' interpolations. When integrating, these integrals going 'deep' (up to 4 layers).

There are two kind of paths that can be set: A path, where the program only requires reading permission (`path_to_tables_readonly`) as well as a path here it requires writing permission (`path_to_tables`).
The argument can either be a string or a list of strings.
If the given string or the first string in the list fulfills the readable/writable requirement, it is set. If the first item in the list does not satisfy the requirement the program is looking for the next item and so on until it reaches the end of the list.
If none of the given elements in the list of strings or the given string is a valid path, then the program writes the tables in the memory.

If both the `path_to_tables_readonly` and `path_to_tables` are valid, PROPOSAL first looks at the readonly path.
If the interpolation table file given by the readonly path doesn't exist or is empty (e.g. if two instances are called in parallel and one process is already building the table but hasn't finished yet), PROPOSAL uses the writing path.
There it again looks, if the interpolation table file already exist or is not empty.
If the tables given by the path have already been built PROPOSAL just uses them.
If there are no tables corresponding to the needed propagation properties PROPOSAL builds the corresponding tables in the folder given by the `path_to_tables`.
If the string is empty, the folder doesn't exist or PROPOSAL has no permission to write, the tables that are needed are stored in the memory.
Note: The tables differ in the parameters given below. These information are stored in the file name. For not too long file names, these values are hashed.

There is the option that just the readonly path should be used (`just_use_readonly_path`). So if there is not the required tables prebuild in the readonly path the Initialization/program wil break and not try to look or write at the `path_to_tables` or in the memory.
When this parameter is enabled but the required tables are not prebuilt in the `path_to_tables_readonly` PROPOSAL will neither look at the `path_to_tables`, nor write the tables in this path nor write the tables in the memory. Instead, the program will stop!

The parameter `do_binary_tables` decides whether the tables are stored as binary files or as a (human readable) text files.

The upper energy limit can be modified (`max_node_energy`) up to the maximum possible primary particle energy, 
to prevent values for particles with energies greater than the maximum energy from being extrapolated.
If particles are propagated with primary energies greater than `max_node_energy`, the interpolation error increases rapidly. 
This should be avoided.

If the error of the interpolation becomes too large, the number of sampling points can be increased by changing the properties `nodes_cross_section`, `nodes_continous_randomization` and `nodes_propagate`. 
This however increases the runtime of PROPOSAL.

| Keyword                         | Type   | Default | Description |
| ------------------------------- | ------ | ------- | ----------- |
| `do_interpolation`              | Bool   | `True`  | Decides, whether to calculate with interpolation tables or integrations |
| `path_to_tables`                | String | `""`    | Path pointing to the folder with the interpolation tables |
| `path_to_tables_readonly`       | String | `""`    | Path pointing to the folder with the interpolation tables with reading permissions only |
| `just_use_readonly_path`        | Bool   | `False` | Decides, if only the readonly path should be used |
| `do_binary_tables`              | Bool   | `True`  | Decides, whether the tables are stored in binary format or in a human readable text format |
| `max_node_energy`               | Double | `1.e14` | Energy in MeV up to which the interpolation tables are built |
| `nodes_cross_section`           | Integer| `100`   | Number of interpolation points for the interpolation of the crosssection integral |
| `nodes_continous_randomization` | Integer| `200`   | Number of interpolation points for the interpolation of the continous randomization integral |
| `nodes_propagate`               | Integer| `1000`  | Number of interpolation points for the interpolation of the propagation integral |

### Accuracy parameters and Scattering ###
There are several parameters with which the precision or speed for advancing the particles can be adjusted.

There are two ways to calculate the elapsed time for the propagated particles chosen with the parameter `exact_time`.
If enabled, the time is calculated taking the energy dependence into account (using the dE/dx information).
Otherwise, the time is calculated by simply divinding the propagated distance by the speed of light. 

When `stopping_decay` is enabled, every particle that gets below a threshold energy (and hasn't decayed yet) is forced to stop and to decay.
This threshold energy is a property of the particle object and is defined in the Particle Definitions. It is by default the same as the particle mass.

There are four different multiple scattering parametrizations describing the deviation of the propagation direction. Those can be set with the `scattering` parameter.
Choosing a different scattering parametrization can significantly alter the computing time. These parametrizations are available:

  - `"Moliere"` Parametrization of [Moliere](http://zfn.mpdl.mpg.de/data/Reihe_A/3/ZNA-1948-3a-0078.pdf)
  - `"Highland"` First Order approximation (a gaussian function) of Moliere's theory derived by [Highland](https://doi.org/10.1016/0029-554X(75)90743-0) and corrected by [Lynch/Dahl](https://doi.org/10.1016/0168-583X(91)95671-Y).
  - `"HighlandIntegral"` From the old PROPOSAL version (which had just this scattering mode). Also using Highland approximation (corrected by Lynch/Dahl), but taking into account the energy dependence of the propagated distance
  - `"NoScattering"` Here the particle always propagate in the initial direction without deviation from the propagation axes.


| Keyword          | Type   | Default     | Description |
| ---------------- | ------ | ----------- | ----------- |
| `exact_time`     | Bool   | `True`      | Decides, whether the energy dependence is considered when calculating the elapsed time |
| `stopping_decay` | Bool   | `True`      | Decides, whether the particle gets stopped and forced to decay, if its energy get below a pre-defined threshold. |
| `scattering`     | String | `"HighlandIntegral"` | Multiple scattering parametrization describing the displacement from the initial propagation direction |


### Cross section parameters ###
There are several parametrizations defining the cross sections and further options which simply scales them.
For most of the interactions, it's possible to choose between multiple parametrizations.

The cross section `multiplier`, available for each cross section, scales this cross section by its factor.

Choosing `"None"` as the option for a parametrizations disables this interaction completely (default option for some non-standard interactions)

For **Ionization**, there are crosssections for electrons, positrons or heavier charged leptons:
  - `"IonizBetheBlochRossi"`
  - `"IonizBergerSeltzerBhabha"` (for positrons)
  - `"IonizBergerSeltzerMoller"` (for electrons)

The **electron pair production** parametrizations are:
  - `"EpairKelnerKokoulinPetrukhin"` (Proc. 12th ICCR (1971), 2436) with corrections for the interaction with atomic electrons (Phys. Atom. Nucl. 61 (1998), 448)
  - `"EpairSandrockSoedingreksoRhode"` 

The **bremsstrahlung** parametrizations are:
  - `"BremsKelnerKokoulinPetrukhin"` ([Preprint MEPhI (1995) no. 024-95](http://cds.cern.ch/record/288828)) and (Phys. Atom. Nucl. 62 (1999), 272)
  - `"BremsAndreevBezrukovBugaev"` (Phys. Atom. Nucl. 57 (1994), 2066)
  - `"BremsPetrukhinShestakov"` (Canad. J. Phys. 46 (1968), 377)
  - `"BremsCompleteScreening"` (Tsai [Rev. Mod. Phys. 46 (1974), 815](https://doi.org/10.1103/RevModPhys.46.815))
  - `"BremsSandrockSoedingreksoRhode"`
  - `"BremsElectronScreening"` (For electrons and positrons with empirical low-energy corrections, SLAC-R-730)

There are two different approaches to parametrise the **nuclear interaction**:

The first way is to work with the approximation where a real photon scatters inelastically with a nucleus with an additional factor to make the photon virtual.
The available parametrizations using this approach are:
- `"PhotoKokoulin"`
- `"PhotoRhode"`
- `"PhotoBezrukovBugaev"` (Sov. J. Nucl. Phys. 33 (1981), 635) with corrections for particles with higher mass like taus ([Phys. Rev. D67 (2003), 034027](https://doi.org/10.1103/PhysRevD.67.034027))
- `"PhotoZeus"` ([Eur. Phys. J. C 7 (1999), 609](https://doi.org/10.1007/s100529901084))

For these parametrizations the hard component can additionally be considered (only affecting the parametrizations using the real photon approximation).

The second way is to use data from electron proton scattering experiments and extrapolate to low momentum of the virtual photon transfered to the nucleus (Q**2).
There the cross section is first integrated over the photon momentum.
The available parametrizations using this approach are:
- `"PhotoAbramowiczLevinLevyMaor91"` by Abramowicz Levin Levy Maor [Phys. Let. B269 (1991), 465](https://doi.org/10.1016/0370-2693(91)90202-2)
- `"PhotoAbramowiczLevinLevyMaor97"` by Abramowicz Levin Levy Maor [arXiv::hep-ph/9712415](https://arxiv.org/abs/hep-ph/9712415)
- `"PhotoButkevichMikhailov"` [JETP 95 (2002), 11](https://doi.org/10.1134/1.1499897)
- `"PhotoRenoSarcevicSu"` [Astrop. Phys. 24 (2005), 107](https://doi.org/10.1016/j.astropartphys.2005.06.002) (which uses the parametrization of ALLM97, but with corrections for spin 0 particles and should therefore only be selected for spin 0 particles like sTau)

For these parametrizations the parametrization of the shadowing factor can be chosen from one of the following parametrizations (only affecting the nuclear interaction parametrizations with momentum integration):
- `"ShadowButkevichMikhailov"` from their calculation of nuclear interaction
- `"ShadowDuttaRenoSarcevicSeckel"` by Dutta, Reno, Sarcevic, Seckel [Phys. Rev. D63 (2001), 094020](https://doi.org/10.1103/PhysRevD.63.094020)

The **LPM effect** (Landau-Pomeranschuk-Migdal), suppressing the bremsstrahlung and the pair production at high energies and the Ter-Mikaelian effect, suppress low bremsstrahlung energy losses, can also be incorporated.

**Annihilation** of positrons with electrons is per default disabled and can be enabled using the parametrization:
  - `"AnnihilationHeitler"` 

The **muon pair production** (which is an optional process and per default disabled) parametrizations are:
  - `"MupairKelnerKokoulinPetrukhin"` (Phys. Atom. Nucl. Vol. 63, No.9 (2000),  pp. 1603-1611, DOI: 10.1134/1.1312894)
  
The **weak interaction** (charged current interaction of a charged lepton, per default disabled) parametrizations are:
  - `"CooperSarkarMertsch"` by Cooper-Sarkar, Mertsch, Sarkar [arXiv:1106.3723 ](https://arxiv.org/abs/1106.3723v1)

| Keyword                  | Type   | Default    | Description |
| -----------------------  | ------ | ---------- | ----------- |
| `brems_multiplier`       | Double | `1.0`        | Scales the bremsstrahlung |
| `epair_multiplier`       | Double | `1.0`        | Scales the electron pair production |
| `ioniz_multiplier`       | Double | `1.0`        | Scales the ionization |
| `photo_multiplier`       | Double | `1.0`        | Scales the nuclear interaction |
| `annihilation_multiplier`| Double | `1.0`        | Scales the annihilation |
| `mupair_multiplier`      | Double | `1.0`        | Scales the muon pair production |
| `weak_multiplier`        | Double | `1.0`        | Scales the weak interaction |
| `brems`                  | String | `"BremsKelnerKokoulinPetrukhin"` | Bremsstrahlung parametrization |
| `epair`                  | String | `"EpairKelnerKokoulinPetrukhin"` | Electron pair production parametrization |
| `ioniz`                  | String | `"IonizBetheBlochRossi"` | Ionization parametrization |
| `photo`                  | String | `"PhotoAbramowiczLevinLevyMaor97"` | Nuclear interaction parametrization |
| `photo_hard_component`   | Bool   | `True`     | Including the hard components |
| `photo_shadow`           | String | `"ShadowButkevichMikhailov"` | Shadowing parametrization |
| `lpm`                    | Bool   | `True`     | Incorporate the LPM-effect and TM-effect |
| `annihilation`           | String | `"None"` | Annihilation parametrization |
| `mupair`                 | String | `"None"` | Muon pair production parametrization |
| `mupair_particle_output` | Bool   | `True`     | Produced muon pairs are treated as particles with corresponding energies in the Output of Secondaries (and not as DynamicData objects) |
| `weak`                   | String | `"None"` | Weak interaction parametrization |

There are also parametrizations that can be used for **Photon propagation**.
[Here](config_photon.md) they are described in detail. 
All parametrizations in connection with photon propagation are per default disabled.

### Energy-cut parameters ###
The energy cut settings and the continuous randomization option are separated between `cuts_inside`, `cuts_infront` and `cuts_behind` the detector, which are again own json-objects.
The parameter `e_cut` describes the cut in the total energy loss in MeV and the parameter `v_cut` describes the cut in the energy loss relative to the particle energy.
Above this cut, the energy losses are calculated stochastically and below this cut, the energy loss is calculated via the average energy loss (continuously).
- If both values are set, the minimum of the total energy cut and the relative cut times the particle energy is applied.
- If just one of the energy cuts should be used, the other value should be set to a negative value (usually it's -1), because negative values for the energy cuts are neglected.
- If both values are set to a negative value, only continuous losses (average energy loss dE/dx) are taking into account.

The continuous randomization randomizes the continuous losses, which affects, that the continuous losses are not always the same for the same particle energy, because they are randomised a little bit.
This behaviour can be enabled with the `cont_rand` parameter.

Note: The energy cuts and the continuous randomization settings can also be specified for each sector.
Then the global settings will be overwritten.

For the `cuts_inside` option, the default values are

| Keyword     | Type   | Default   | Description |
| ----------- | ------ | --------- | ----------- |
| `e_cut`     | Double | `500.0`   | Total energy loss cut inside the detector |
| `v_cut`     | Double | `-1.0`    | Relative energy loss cut inside the detector |
| `cont_rand` | Bool   | `True`    | Includes the continuous randomization inside the detector |

For the `cuts_infront` option, the default values are

| Keyword     | Type   | Default   | Description |
| ----------- | ------ | --------- | ----------- |
| `e_cut`     | Double | `-1.0`    | Total energy loss cut in front the detector |
| `v_cut`     | Double | `0.001`   | Relative energy loss cut in front the detector |
| `cont_rand` | Bool   | `True`    | Includes the continuous randomization in front the detector |

For the `cuts_behind` option, the default values are

| Keyword     | Type   | Default   | Description |
| ----------- | ------ | --------- | ----------- |
| `e_cut`     | Double | `-1.0`    | Total energy loss cut behind the detector |
| `v_cut`     | Double | `-1.0`    | Relative energy loss cut behind the detector |
| `cont_rand` | Bool   | `False`   | Includes the continuous randomization behind the detector |


## The `sectors` configurations ##

Each sector needs a medium, a geometry and a priority. Furthermore, the energy cut parameters of a sector can be set to override the global settings.

### Hierarchy ###

The priority of a sector compared to other sectors can be set with the `hierarchy` option.
If multiple sectors are defined and if these sectors overlap, the decision which sector should be preferred is made with the hierarchy option.
The bigger hierarchy is preferred. If the overlapping sectors have the same hierarchy, then the one with the more dense medium is chosen. 
If both the hierarchy and the mass density are equal, then the first sector in the list is chosen.

| Keyword     | Type    | Default | Description |
| ----------- | ------- | ------- | ----------- |
| `hierarchy` | Integer | `1`     | Decides, which sector is favored, when they overlap |

### Medium parameters ###
The following media are implemented:
- `"Hydrogen"`
- `"Iron"`
- `"Copper"`
- `"Lead"`
- `"Uranium"`
- `"Air"`
- `"Water"`
- `"Ice"`
- `"StrandardRock"`
- `"FrejusRock"`
- `"Salt"`
- `"Paraffin"`
- `"AntaresWater"`

A density correction factor for the medium, specifying the density of the medium for the experiment, can optionally be set.

| Keyword              | Type   | Default | Description |
| -------------------- | ------ | ------- | ----------- |
| `medium`             | String | `"Water"` | Medium in this sector |
| `density_correction` | Double | `1.0`   | Density correction factor of the medium |

### Geometry parameters ###

Every sector needs a geometry, described by an vector that defines its `origin` and a `shape`.
All lengths PROPOSAL uses are in **centimeter**!
(Unit of length in the setting below are in **meter**!)

The following three shapes are available:
- `"Sphere"`
- `"Cylinder"`
- `"Box"`

| Keyword  | Type                     | Default           | Description |
| -------- | ------------------------ | ----------------- | ----------- |
| `shape`  | String                   |  `-`              | Shape of the Geometry |
| `origin` | [Double, Double, Double] | `[0.0, 0.0 ,0.0]` | Center of the Geometry |

Every shape has its own specifications:
- A sphere needs an inner and an outer radius between which the sector is defined
- A cylinder needs also an inner and an outer radius between which the sector is defined and a height (there is just one height, not an inner and an outer one)
- A box needs a length, a width and a height. Unlike the other two geometries, where the sector is defined between two shapes, the box has just one shape and the sector is inside this.

| Keyword        | Type   | Default | Description |
| -------------- | ------ | ------- | ----------- |
| `outer_radius` | Double | `-`     | Outer radius of the sphere or the cylinder |
| `inner_radius` | Double | `-`     | Inner radius of the sphere or the cylinder |
| `height`       | Double | `-`     | Height of the cylinder or the box (length in z-direction) |
| `length`       | Double | `-`     | Length of the box (length in x-direction) |
| `width`        | Double | `-`     | Width of th box (length in y-direction) |

### Energy cut parameters ###

The energy cut parameters can be specified for every sector, which then overwrite the globally defined cut settings.
These cuts use the same specifications like the global ones - For a description see the previous section.

## The `detector` configurations ##

The `detector` is defined by a geometry, its parameters can be specified like the geometry of a sector described in the previous section.
