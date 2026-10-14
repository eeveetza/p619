# MATLAB/Octave Implementation of Recommendation ITU-R P.619

<!--[![DOI](https://zenodo.org/badge/459560677.svg)](https://zenodo.org/badge/latestdoi/459560677) -->

This development code repository contains a MATLAB/Octave software implementation of [Recommendation ITU-R P.619-6](https://www.itu.int/rec/R-REC-P.619/en) (09/2025) "Propagation data required for the evaluation of interference between stations in space and those on the surface of the Earth".

The implementation covers the basic transmission loss for single-entry interference (§ 3.1, equation (14)) and for multiple-entry interference (§ 3.2, equation (15)), together with the propagation mechanisms of § 2 and Attachments A to F. It is validated against a set of validation examples (see [Validation](#validation)). The implementation and validation are work in progress and has not yet been formally reviewed; see [Status and limitations](#status-and-limitations).

## Contents

- [MATLAB/Octave Implementation of Recommendation ITU-R P.619](#matlaboctave-implementation-of-recommendation-itu-r-p619)
  - [Contents](#contents)
  - [Repository structure](#repository-structure)
  - [Integrating ITU digital products](#integrating-itu-digital-products)
  - [Basic transmission loss for single-entry interference (§ 3.1)](#basic-transmission-loss-for-single-entry-interference--31)
    - [Function call](#function-call)
    - [Required input arguments](#required-input-arguments)
    - [Optional input arguments](#optional-input-arguments)
    - [Outputs](#outputs)
    - [Examples](#examples)
  - [Basic transmission loss for multiple-entry interference (§ 3.2)](#basic-transmission-loss-for-multiple-entry-interference--32)
    - [Function call](#function-call-1)
    - [Optional input arguments](#optional-input-arguments-1)
    - [Outputs](#outputs-1)
    - [Example](#example)
  - [Methods by section](#methods-by-section)
  - [Implementation notes](#implementation-notes)
  - [Validation](#validation)
  - [Status and limitations](#status-and-limitations)
  - [Software versions](#software-versions)
  - [References](#references)

## Repository structure

| File/Folder | Description |
|---|---|
| `matlab/P619.m` | MATLAB class with methods implementing Recommendation ITU-R P.619-6 |
| `matlab/validate_p619.m` | Runs the validation vectors in `validation/vectors` against `P619.m` |
| `matlab/test_p619.m` | Development tests and plots reproducing figures of the Recommendation |
| `matlab/initiate_digital_maps.m` | Converts ITU-R digital maps into MATLAB functions (see below) |
| `matlab/private/` | Map functions/files generated from the ITU-R digital maps (P.453 N<sub>wet</sub> maps are included) and helper functions |
| `matlab/validation_examples/` | ITU-R Study Group 3 validation examples (Excel workbook) |
| `validation/vectors/` | Implementation-independent validation vectors (CSV) and tolerances |
<!-- | `validation/reference/` | Independent Python reference implementation and vector generator | -->
| `validation/README.md` | Description of the validation suite and of all conventions |

## Integrating ITU digital products

This software uses ITU digital products that are an integral part of Recommendations. These products must not be reproduced or distributed without explicit written permission from the ITU.

| Map | Recommendation | Used for | Status |
|---|---|---|---|
| N<sub>wet</sub> (annual statistics, 18 percentages) | ITU-R P.453-14 | Median N<sub>wet</sub> for tropospheric scintillation (Attachment D) in `tl_p619_single`, unless `'Nwet'` is given | included in `matlab/private` (`DigitalMaps_NWET_Annual_xx.m`); can also be generated from `NWET_Annual_xx.TXT` |
| `DN50.TXT` | ITU-R P.452-18 | Median effective Earth radius for the ducting term of § 2.6, only when a terrain profile is given and `'ae'` is not | to be downloaded |
| `RHO_xx`, `VSCH_xx` (annual, 18 percentages), `TOPO_0DOT5` | ITU-R P.836-6 | Surface water-vapour density for A<sub>g</sub>(p<sub>1</sub>), only with `'waterVapour', 'P.836'` | to be downloaded |

To use these maps:

1. Download the digital maps of [ITU-R P.452-18](https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.452-18-202310-I!!ZIP-E.zip) and extract `DN50.TXT` to `matlab/private/maps`.
2. Download the digital maps of [ITU-R P.836-6](https://www.itu.int/rec/R-REC-P.836-6-201712-I/en) and extract the annual surface water-vapour density maps (`RHO_01` … `RHO_99`), water-vapour scale height maps (`VSCH_01` … `VSCH_99`) and `TOPO_0DOT5` to `matlab/private/maps` (sub-folders are allowed; names are matched case-insensitively by their prefix).
3. Optional: download the digital maps of [ITU-R P.453-14](https://www.itu.int/rec/R-REC-P.453-14-201908-I/en) and extract `NWET_Annual_01.TXT` … `NWET_Annual_99.TXT` to `matlab/private/maps`. This is only needed if the `DigitalMaps_NWET_Annual_xx.m` functions are not present.
4. Run `initiate_digital_maps.m`. It creates `matlab/private/DigitalMaps_DN50.m`, `matlab/private/DigitalMaps_P836.mat` and `matlab/private/DigitalMaps_P453.mat`; maps that are not present are skipped with a warning.

The N<sub>wet</sub> maps are read once per MATLAB session (`private/p453_nwet_maps.m`), from `DigitalMaps_P453.mat` if it exists and otherwise from the `DigitalMaps_NWET_Annual_xx.m` functions.

## Basic transmission loss for single-entry interference (§ 3.1)

Equation (14):

L<sub>b</sub> = L<sub>bfs</sub> + A<sub>xp</sub> + A<sub>g</sub>(p<sub>1</sub>) + A<sub>bs</sub> + A<sub>s</sub>(p<sub>2</sub>) + L<sub>dtb</sub>(p<sub>1</sub>)

### Function call

```matlab
p = P619;
L = p.tl_p619_single(f, dkm, e2sflag, p1, He, Hs, phi_e, phi_s, Dphi, atm_type, ...
                     lat_s, lat_e, lon_s, lon_e, p2, Ga);
L = p.tl_p619_single(..., p2, Ga, d, h);                 % with a terrain profile (§ 2.6)
L = p.tl_p619_single(..., p2, Ga, Name, Value);          % name-value pairs, no profile
L = p.tl_p619_single(..., p2, Ga, d, h, Name, Value);
```

### Required input arguments

| Variable | Type | Units | Limits | Description |
|---|---|---|---|---|
| `f` | scalar double | GHz | 0.1 ≤ `f` ≤ 100 | Frequency |
| `dkm` | scalar double or `[]` | km | | Path length. `[]`: straight-line distance from Attachment A |
| `e2sflag` | logical | | | `true`: Earth-to-space (ascending ray, Att. C.5); `false`: space-to-Earth (descending ray, Att. C.4) |
| `p1` | scalar double | % | 0.001 ≤ `p1` ≤ 50 | Time percentage for A<sub>g</sub> and L<sub>dtb</sub> |
| `He` | scalar double | km | | Altitude of the earth-based station antenna above sea level |
| `Hs` | scalar double | km | `Hs` > `He` | Altitude of the space station above sea level (may be above 100 km) |
| `phi_e` | scalar double or `[]` | deg | | Elevation of the main beam of the earth-based station. `[]`: pointing along the path |
| `phi_s` | scalar double or `[]` | deg | | Elevation of the main beam of the space station. `[]`: pointing along the path |
| `Dphi` | scalar double | deg | | Half-power beamwidth of the space station (`e2sflag` true) or the earth-based station (false). 360: no beamwidth check |
| `atm_type` | integer | | 1, 2, 31, 32, 41, 42 | ITU-R P.835 reference atmosphere: 1 mean annual global, 2 low-latitude, 31/32 mid-latitude summer/winter, 41/42 high-latitude summer/winter |
| `lat_s`, `lon_s` | scalar double | deg | | Latitude and longitude of the sub-satellite point |
| `lat_e`, `lon_e` | scalar double | deg | | Latitude and longitude of the earth-based station |
| `p2` | scalar double | % | 0 < `p2` < 100 | Time percentage for scintillation (50 %: A<sub>s</sub> ≈ 0) |
| `Ga` | scalar double | dBi | | Earth-based antenna gain in the direction of the path (Attachment D) |

### Optional input arguments

| Variable | Type | Units | Default | Description |
|---|---|---|---|---|
| `d` | vector double | km | `[]` | Terrain profile distances from the earth-based station (`d(1)` = 0) along the azimuth of the path, extending beyond the station's horizon |
| `h` | vector double | m | `[]` | Terrain heights above sea level at `d` |
| `'xpdMode'` | string | | `'none'` | Depolarization (§ 2.2): `'none'` (A<sub>xp</sub> = 0), `'xpd'`, `'faraday'` or `'hydro'` |
| `'XPDValue'` | scalar double | dB | | Cross-polarization discrimination for `'xpd'` |
| `'Bav'`, `'NT'` | scalar double | T, el/m² | | Earth magnetic field and total electron content for `'faraday'` |
| `'Tp'`, `'Ap'`, `'tau'`, `'theta'` | scalar double | %, dB, deg, deg | | Inputs of the ITU-R P.618 hydrometeor XPD for `'hydro'` |
| `'ae'` | scalar double | km | DN50 map | Median effective Earth radius for the horizon parameters of § 2.6 |
| `'Ld'` | scalar double | dB | `[]` | Diffraction loss of a specific obstruction (e.g. from `finite_width_screen_loss`); if given, L<sub>dtb</sub> = L<sub>d</sub> |
| `'terrainMethod'` | string | | `'slant'` | Terrain diffraction: `'slant'` (ITU-R P.526-16 § 4.6) or `'knife-edge'` (single knife edge) |
| `'scintillation'` | logical | | `true` | `false`: A<sub>s</sub> is not evaluated |
| `'waterVapour'` | string | | `'reference'` | Water-vapour density for A<sub>g</sub>(p<sub>1</sub>) (`atm_type` = 1): `'reference'` (ρ<sub>0</sub> = 7.5 g/m³) or `'P.836'` (surface density exceeded for 100 − p<sub>1</sub> % at the earth-based station, ITU-R P.836-6) |
| `'rho0'` | scalar double | g/m³ | `[]` | Sea-level surface water-vapour density of the profile ρ(h) = ρ<sub>0</sub> exp(−h/2) (`atm_type` = 1); overrides `'waterVapour'` |
| `'Nwet'` | scalar double | ppm | P.453 map | Wet term of the surface refractivity for tropospheric scintillation; default: median (50 %) value of the ITU-R P.453-14 maps at the earth-based station (from local temperature and humidity: `p453_Nwet`) |

### Outputs

`L` is a structure with the fields:

| Field | Units | Description |
|---|---|---|
| `Lb` | dB | Basic transmission loss, equation (14) |
| `Lbfs` | dB | Free-space basic transmission loss, § 2.1 |
| `Ax` | dB | Attenuation due to depolarization, § 2.2 |
| `Ag` | dB | Attenuation due to atmospheric gases, Attachment C |
| `Abs` | dB | Beam-spreading loss, § 2.4.2 |
| `Ast` | dB | Tropospheric scintillation, Attachment D (negative: enhancement) |
| `Ldtb`, `Ld` | dB | Diffraction/ducting loss and diffraction loss, § 2.6 |
| `th_hoz`, `d_hoz`, `nu_max` | mrad, km, – | Horizon elevation and distance, diffraction parameter (§ 2.6) |
| `dkm` | km | Path length |
| `theta_0`, `theta` | deg | Free-space and apparent elevation angle (Attachments A, B) |
| `rho0`, `rho_s`, `q_rho` | g/m³, g/m³, % | Sea-level water-vapour density used; P.836 surface density at the station and its exceedance percentage (NaN if not used) |
| `Nwet` | ppm | Wet term of the surface refractivity used for A<sub>s</sub> (NaN if scintillation is not evaluated) |

### Examples

```matlab
p = P619;

% Bern earth station (0.55 km asl) to a GSO satellite at 0°E, 20 GHz, Earth-to-space,
% antennas pointing along the path, mean annual global atmosphere
L = p.tl_p619_single(20, [], true, 50, 0.55, 35786, [], [], 360, 1, ...
                     0, 46.95, 0, 7.44, 50, 40);
% L.Lb = 210.44 dB, L.Lbfs = 210.10 dB, L.Ag = 0.3404 dB, L.dkm = 38134 km, L.theta_0 = 35.536°

% Tromsø (0.1 km asl) to a GSO satellite at 20°W, p1 = 0.01 %, with a ridge 8 km away
d = (0:0.05:30)';                          % km
h = 1150*exp(-((d - 8)/1.5).^2);           % m asl
L = p.tl_p619_single(20, [], true, 0.01, 0.1, 35786, [], [], 360, 1, ...
                     0, 69.65, -20, 18.96, 50, 40, d, h, 'ae', 6371*4/3);
% L.theta = 7.195°, L.nu_max = 5.50, L.Ld = 41.60 dB, L.Lb = 254.27 dB
% (with 'terrainMethod', 'knife-edge': L.Ld = 27.74 dB)

% Time-percentage dependence of the gaseous attenuation (P.836-6 maps required):
% Ag not exceeded for p1 = 10 % uses the density exceeded for 90 % of the time
L = p.tl_p619_single(20, [], true, 10, 0.55, 35786, [], [], 360, 1, ...
                     0, 46.95, 0, 7.44, 50, 40, 'waterVapour', 'P.836');
% L.q_rho = 90 %, L.rho_s = 4.189 g/m³, L.rho0 = 5.515 g/m³, L.Ag = 0.2738 dB (0.3404 dB with 7.5 g/m³)
```

## Basic transmission loss for multiple-entry interference (§ 3.2)

Equation (15), for one path of a multiple-entry (e.g. Monte Carlo) simulation:

L<sub>b</sub> = L<sub>bfs</sub> + A<sub>xp</sub> + A<sub>g</sub>(p) + A<sub>bs</sub> + L<sub>c</sub>(p<sub>Lc</sub>) + L<sub>be</sub>(p<sub>Lbe</sub>) + L<sub>dtb</sub>(p)

Scintillation is not included. Unless `'pLc'` and `'pLbe'` are given, new uniformly distributed random percentages of locations are drawn on each call (set the state of `rand` for reproducible simulations). The aggregate interference over all paths is the sum in linear units, `aggregate_power`.

### Function call

```matlab
L = p.tl_p619_multiple(f, dkm, e2sflag, p, He, Hs, phi_e, phi_s, Dphi, atm_type, ...
                       lat_s, lat_e, lon_s, lon_e);
L = p.tl_p619_multiple(..., lon_e, d, h);
L = p.tl_p619_multiple(..., lon_e, Name, Value);
L = p.tl_p619_multiple(..., lon_e, d, h, Name, Value);
Pagg = p.aggregate_power(P);       % P: received powers per path (dBW)
```

The required inputs are those of `tl_p619_single` without `p2` and `Ga`; `p` (%) is the time percentage for A<sub>g</sub> and L<sub>dtb</sub>.

### Optional input arguments

| Variable | Type | Units | Default | Description |
|---|---|---|---|---|
| `d`, `h` | vector double | km, m | `[]` | Terrain profile, as for `tl_p619_single` |
| `'Axp'` | scalar double | dB | 0 | Attenuation due to polarization (3 dB for arbitrary polarization, § 2.2) |
| `'clutter'` | logical | | `false` | Earth-based station within urban/suburban clutter (`false`: L<sub>c</sub> = 0) |
| `'clutterModel'` | string | | `'P.2108-2-draft'` | `'P.2108-2-draft'`: draft revision ITU-R P.2108-2 § 3.3 (0.5–100 GHz); `'P.2108-1'`: ITU-R P.2108-1 § 3.3 (10–100 GHz) |
| `'hs'` | scalar double | m | | Ground station height above ground (≥ 1 m); required by `'P.2108-2-draft'` |
| `'hm'` | scalar double | m | | Median clutter height (low-rise ≤ 8, mid-rise 8–20, high-rise > 20); required by `'P.2108-2-draft'` |
| `'pLc'` | scalar double | % | random | Percentage of locations for L<sub>c</sub> |
| `'indoor'` | logical | | `false` | Indoor earth-based station (`false`: L<sub>be</sub> = 0) |
| `'buildingType'` | integer | | 1 | ITU-R P.2109: 1 traditional, 2 thermally efficient |
| `'pLbe'` | scalar double | % | random | Probability for L<sub>be</sub> |
| `'ae'`, `'Ld'`, `'terrainMethod'` | | | | As for `tl_p619_single` (§ 2.6) |
| `'waterVapour'`, `'rho0'` | | | | As for `tl_p619_single` (A<sub>g</sub>(p), with p in place of p<sub>1</sub>) |

### Outputs

`L` is a structure with the fields `Lb`, `Lbfs`, `Axp`, `Ag`, `Abs`, `Lc`, `Lbe`, `Ldtb` (dB), `pLc`, `pLbe` (%, the percentages used; NaN if not applicable), `dkm` (km), `theta_0` and `theta` (deg), `rho0`, `rho_s` (g/m³) and `q_rho` (%).

### Example

```matlab
% One draw of a Monte Carlo simulation: indoor station in mid-rise clutter
rand('seed', 1);
L = p.tl_p619_multiple(20, [], true, 50, 0.55, 35786, [], [], 360, 1, ...
                       0, 46.95, 0, 7.44, 'Axp', 3, 'clutter', true, 'hs', 10, 'hm', 15, ...
                       'indoor', true, 'buildingType', 1);
% L.Lb = 282.97 dB (L.Lc = 26.33 dB at L.pLc = 86.68 %, L.Lbe = 43.20 dB at L.pLbe = 91.27 %)
```

## Methods by section

All methods are called on an instance of the class, e.g. `p = P619; p.tl_free_space(20, 38134)`.

| Method | Reference | Description |
|---|---|---|
| `tl_p619_single` | § 3.1, eq. (14) | Basic transmission loss for single-entry interference |
| `tl_p619_multiple` | § 3.2, eq. (15) | Basic transmission loss for multiple-entry interference |
| `aggregate_power` | § 3.2 | Aggregate interference power (linear sum) |
| `tl_free_space` | § 2.1, eq. (1) | Free-space basic transmission loss |
| `co_and_cross_polar_attenuation` | § 2.2.1, eqs. (2a), (2b) | Co- and cross-polar attenuation from XPD |
| `faraday_co_and_cross_polar_attenuation` | § 2.2.2, eqs. (3a), (3b), (4) | Faraday rotation angle and attenuations |
| `p618_hydrometeor_xpd`, `p618_xpd_scaling`, `hydrometeor_depolarization_attenuation` | § 2.2.3, eq. (6); ITU-R P.618 § 4.1 | Hydrometeor depolarization |
| `beam_spreading_loss` | § 2.4.2; ITU-R P.834 § 5 | Beam-spreading loss (θ<sub>0</sub> ≤ 10°, h ≤ 5 km; 0 outside) |
| `tropospheric_scintillation` | Attachment D, eqs. (48)–(50); ITU-R P.618-12 § 2.4.1 | Tropospheric scintillation A<sub>st</sub> and intensity σ |
| `fresnel_parameters` | § 2.6, eqs. (12a), (12b) | Diffraction parameter ν and first Fresnel zone radius |
| `knife_edge_loss` | § 2.6; ITU-R P.526-16 eq. (31) | Single knife-edge loss J(ν) |
| `finite_width_screen_loss` | § 2.6; ITU-R P.526-16 § 5.1, eqs. (71), (72) | Finite-width screen (e.g. building), J<sub>min</sub> and J<sub>av</sub> |
| `slant_path_terrain_loss` | § 2.6; ITU-R P.526-16 § 4.6 | General Earth–space slant-path terrain loss (Bullington, effective plane-Earth, blending) |
| `terrain_obstruction_loss` | § 2.6, Attachment E | Terrain diffraction loss along the refracted ray (`'slant'` or `'knife-edge'`) |
| `diffraction_ducting_loss` | § 2.6, eqs. (13)–(13e); ITU-R P.452-18 | Terrain L<sub>d</sub>, horizon parameters and ducting term L<sub>dtb</sub> |
| `tl_dtb`, `p452_path_profile` | § 2.6, eq. (13); ITU-R P.452-18 eqs. (137), (141) | Ducting-enhanced diffraction loss; horizon angle and distance |
| `cl_p2108_3` | § 2.7; draft ITU-R P.2108-2 § 3.3 | Earth–space clutter loss (default) |
| `cl_p2108_1` | § 2.7; ITU-R P.2108-1 § 3.3 | Earth–space clutter loss (P.2108-1) |
| `bel_p2109` | § 2.8; ITU-R P.2109-2 | Building entry loss |
| `straight_line_e2s` | Attachment A, eqs. (18)–(24) | Distance, free-space elevation and azimuth of a straight-line Earth–space path |
| `fs2apparent`, `apparent2fs` | Attachment B, eqs. (25)–(28) | Conversion between free-space and apparent elevation angles |
| `atm_attenuation_E2s` | Attachment C.5, eqs. (39)–(46) | Gaseous attenuation, ascending ray (cases φ<sub>e</sub> ≥ 0 and φ<sub>e</sub> < 0) |
| `atm_attenuation_s2E` | Attachment C.4, eqs. (33)–(38) | Gaseous attenuation, descending ray |
| `solve45` | Attachment C.5, eqs. (45)–(47) | Altitude H<sub>min</sub> by Newton iteration |
| `p676d11_ga` | ITU-R P.676 Annex 1, eq. (1) | Specific attenuation due to oxygen and water vapour |
| `p676_slant_path_geometry16` | ITU-R P.676 Annex 1, eq. (16) | Layer boundaries between two heights |
| `p835_std_atm_profiles` | ITU-R P.835 | T, P, ρ and refractive index of the reference atmospheres |
| `p453_n`, `p453_Nwet` | ITU-R P.453 | Refractive index, wet term of refractivity |
| `p836_rho_surface` | § 2.3; ITU-R P.836-6 Annex 1, ITU-R P.1144 | Surface water-vapour density exceeded for p % at a location and altitude |
| `beam_clearance` | Attachment E, eqs. (51)–(59) | Ray height profile with atmospheric refraction |
| `test_precipitation_scatter` | Attachment F, eqs. (60a)–(67) | Test whether a full precipitation-scatter calculation is needed, with all intermediate values |
| `specific_rain_attenuation_p838` | ITU-R P.838-3 | Specific rain attenuation (path inclination in radians) |

## Implementation notes

The points below are choices made where the text of the Recommendation is ambiguous or appears to contain an error. All of them are documented in more detail in [validation/README.md](validation/README.md).

* **Attachment A, equations (23) and (24).** The elevation angle is computed as arctan(Z<sub>2</sub>/G<sub>ts</sub>) and the bearing from arctan(Y<sub>2</sub>/X<sub>2</sub>). P.619-6 defines atan2(x, y) = arctan(x/y), which applied literally to equations (23) and (24) gives the zenith angle and a wrong bearing.
* **Attachment C.** Steps 4 to 6 are evaluated for all non-nadir paths. Layers follow ITU-R P.676 equation (16) up to 100 km; a space station above 100 km is treated as being in free space (n<sub>s</sub> = 1). The P.835 mixing-ratio floor (e/P ≥ 2·10<sup>−6</sup>) is applied to the water-vapour density. Case 2 (φ<sub>e</sub> < 0) uses the Newton iteration (46)–(47); rays that would pass below sea level return 10<sup>20</sup> dB.
* **§ 2.6.** The ray from Attachment E is the direct ray of ITU-R P.526. For the slant-path method (P.526-16 § 4.6, default), heights relative to the ray are converted to along-ray and perpendicular coordinates with the path elevation. The single knife-edge method (largest ν along the ray) remains available. Equation (13b) uses the coefficient 0.198. Equation (13) is applied to terrain obstruction only.
* **Attachment F, equation (60b).** `2*r1/sind(theta)` is used instead of the printed 2r<sub>1</sub> sin θ (suspected typographical error, reported to WP 3M): only the division keeps equations (65a) and (65b) consistent for beams of equal radius.
* **§ 3.2.** L<sub>c</sub> and L<sub>be</sub> use the apparent elevation angle of the path. The draft revision P.2108-2 § 3.3 is the default clutter model.
* **Tropospheric scintillation** (Attachment D, method of ITU-R P.618-12 § 2.4.1) uses the median (50 %) N<sub>wet</sub> of the ITU-R P.453-14 annual maps at the earth-based station (bilinear interpolation on the 0.75° grid; between map percentages, linear in log p), as in ITU-R P.618-14 and the ITU-R SG3 validation examples. The time percentage p<sub>2</sub> applies only to the fade/enhancement, not to N<sub>wet</sub>.
* **Beam-spreading loss** is set to 0 outside the validity range of ITU-R P.834 § 5 (θ<sub>0</sub> > 10° or h > 5 km).
* **Gaseous attenuation** uses the line-by-line method of ITU-R P.676 Annex 1 with the spectroscopic data of P.676-11, as referenced in P.619-6.
* **Time-percentage dependence of A<sub>g</sub>(p<sub>1</sub>)** (§ 2.3, § 4). Attachment C itself has no time percentage. Following § 2.3 (P.836 maps exceeded for a percentage of time; scaling ρ = ρ<sub>0</sub> exp(−h/2)) and § 4 (ρ<sub>0</sub> as input to Attachment C), A<sub>g</sub>(p<sub>1</sub>), the attenuation *not exceeded* for p<sub>1</sub> % of the time, is computed with the water-vapour density *exceeded* for q = 100 − p<sub>1</sub> % from ITU-R P.836-6 at the earth-based station (altitude H<sub>e</sub>), converted to the sea-level value ρ<sub>0</sub> = ρ<sub>s</sub> exp(H<sub>e</sub>/2) of the P.835 mean annual global reference profile. The P.836 maps cover 0.1 ≤ q ≤ 99 %; for p<sub>1</sub> < 1 % the 99 % value is used (with a warning). Only the water-vapour profile varies; temperature and dry pressure remain those of the reference atmosphere.

## Validation

`validation/` contains implementation-independent test vectors (CSV) with inputs, expected outputs and tolerances. Where available, the expected values are taken from, or checked against, the ITU-R Study Group 3 validation examples (ITU-R P.676-13, P.618-14, P.838-3) and the SG3 Clutter and BEL workbook (draft P.2108-2).

```matlab
cd matlab
results = validate_p619();
```

All 2791 values in the 29 vector files pass (with the P.836-6 maps processed; without them, the 162 values that need the maps are reported as not provided). To validate another implementation, only the adapter functions (`run_*`) at the end of `validate_p619.m` need to be rewritten. See [validation/README.md](validation/README.md) for details.

## Status and limitations

Not yet implemented:

* Ionospheric scintillation, § 2.5.1 (Recommendation ITU-R P.531)
* Differential rain attenuation, § 2.10 (Recommendation ITU-R P.1815)
* Multiple-entry aggregation methods of § 5 (Monte Carlo and analytical); `aggregate_power` performs only the linear summation of § 3.2

Other remarks:

* Case 2 of Attachment C (φ<sub>e</sub> < 0) is implemented for the reference standard atmospheres only.
* The clutter model of the draft revision ITU-R P.2108-2 is not yet published; P.2108-1 can be selected with `'clutterModel', 'P.2108-1'`.

## Software versions

The code was tested and runs on:

* MATLAB version 2022a
* Octave version 8.4.0 (validation suite)

<!-- The validation reference implementation requires Python ≥ 3.9 with `numpy`, `scipy`, `itur` (ITU-Rpy) and `openpyxl`. -->

## References

* [Recommendation ITU-R P.619](https://www.itu.int/rec/R-REC-P.619/en)
* [Recommendation ITU-R P.526](https://www.itu.int/rec/R-REC-P.526/en)
* [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)
* [Recommendation ITU-R P.453](https://www.itu.int/rec/R-REC-P.453/en)
* [Recommendation ITU-R P.618](https://www.itu.int/rec/R-REC-P.618/en)
* [Recommendation ITU-R P.676](https://www.itu.int/rec/R-REC-P.676/en)
* [Recommendation ITU-R P.834](https://www.itu.int/rec/R-REC-P.834/en)
* [Recommendation ITU-R P.835](https://www.itu.int/rec/R-REC-P.835/en)
* [Recommendation ITU-R P.836](https://www.itu.int/rec/R-REC-P.836/en)
* [Recommendation ITU-R P.838](https://www.itu.int/rec/R-REC-P.838/en)
* [Recommendation ITU-R P.1144](https://www.itu.int/rec/R-REC-P.1144/en)
* [Recommendation ITU-R P.2108](https://www.itu.int/rec/R-REC-P.2108/en), draft revision P.2108-2: [Document 3/41](https://www.itu.int/md/meetingdoc.asp?lang=en&parent=R23-SG03-C-0041), implementation: [eeveetza/p2108](https://github.com/eeveetza/p2108/tree/dev)
* [Recommendation ITU-R P.2109](https://www.itu.int/rec/R-REC-P.2109/en)
