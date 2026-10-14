# Validation suite for Recommendation ITU-R P.619-6

Implementation-neutral test vectors for Recommendation ITU-R P.619-6 (09/2025), intended to
check this MATLAB/Octave implementation and any other implementation (Python, C++, …) against
the same reference values.

```
matlab/
└── validate_p619.m     runner for P619.m (adapters per vector file)
validation/
├── README.md           this file
├── vectors/            CSV test vectors (inputs + expected outputs) and tolerances.csv
└── reference/          independent Python reference implementation + vector generator
```

## Running

**MATLAB / Octave** (from `matlab/`):

```matlab
results = validate_p619();
```

The runner prints, per file and quantity, the number of values passed, the number the
implementation does not provide (`n/a`), the largest deviation and the row id where it occurs.

**Another implementation:** read the CSVs, evaluate the listed inputs, and compare each output
column with `|impl − ref| ≤ abs_tol + rel_tol·|ref|` using `vectors/tolerances.csv`. In
`validate_p619.m` everything implementation-specific is in the adapter functions (`run_*`), so a
port of the runner only needs new adapters.

**Regenerating the vectors** (Python ≥ 3.9 with `numpy scipy itur openpyxl`):

```
cd validation/reference
python3 generate_vectors.py [path/to/CG-3M3J-13-ValEx-Rev8.1.1.xlsx]
```

## Where the expected values come from

| Source | Used for |
|---|---|
| **ITU-R SG 3 validation examples** (`CG-3M3J-13-ValEx-Rev8.1.1.xlsx`, sheets *P.676-13 SpAtt*, *A_Gas_A1_2.2.1a/b*, *P.618-14 A_Scint*) | specific attenuation, three layered slant paths, 48 scintillation cases |
| **ITU-Rpy** (`itur` 0.4, P.676-12 line tables) | specific attenuation at further p, ρ, T, f |
| **`reference/p619_ref.py`** – written from the Recommendation texts, not ported from `P619.m` | everything else |

The reference implementation was itself checked before generating vectors:

* slant-path Ag vs the three official P.676-13 examples: ≤ 5·10⁻⁷ dB (the residual is the official
  sheet using the eq. (13)–(15) layer grid rather than eq. (16));
* γo, γw vs the official sheet: all printed digits;
* MAGRA T, P vs ITU-Rpy (P.835-6): exact;
* Attachment A vs an independent ECEF/ENU look-angle computation: ≤ 2·10⁻¹² km / deg;
* Attachment D vs the 48 official P.618-14 rows (fade side, with D_eff = √η·D): ≤ 10⁻¹³ relative;
* eq. (45) root by Brent's method, i.e. a different algorithm from the Newton iteration (46)–(47).

Agreement with the reference therefore means agreement between two independently written
implementations plus the official ITU examples – not a proof of correctness. Where the
Recommendation text is ambiguous, the choice made is listed under *Conventions*.

## Vector files

| File | P.619-6 | Inputs → outputs |
|---|---|---|
| `01_free_space.csv` | §2.1 eq. (1) | f, d → Lbfs |
| `02_depolarisation.csv` | §2.2.1, §2.2.3 eqs. (2a), (2b), (6) | XPD → Ax, Ac, Axq |
| `03_faraday.csv` | §2.2.2 eqs. (3a), (3b), (4) | f, Bav, NT → θF, AxF, AcF |
| `04_beam_spreading.csv` | §2.4.2 (P.834 §5) | θ0, h → Abs |
| `05_geometry_attA.csv` | Att. A eqs. (18)–(24) | Hs, Ht, φs, φt, δ → Dts, θ0, ψ |
| `06_attB_fs_to_apparent.csv` | Att. B eqs. (25)–(26) | Ht, θ0 → θ |
| `07_attB_apparent_to_fs.csv` | Att. B eqs. (27)–(28) | Ht, θ → θ0 |
| `08_specific_attenuation.csv` | Att. C step 5 (P.676 Annex 1 eq. (1)) | f, p_dry, ρ, T → γo, γw |
| `09_atmosphere_magra.csv` | Att. C (P.835-6 MAGRA, P.453) | h → T, P, ρ, e, n |
| `10_ag_ascending.csv` | Att. C.5 Case 1, eqs. (39)–(40) | f, He, Hs, φe → Ag |
| `11_ag_descending.csv` | Att. C.4, eqs. (33)–(34) | f, He, Hs, φe (arrival), φs (launch) → Ag |
| `12_hmin.csv` | Att. C.5 Case 2, eq. (45) | He, φe → Hmin |
| `13_ag_negative_elevation.csv` | Att. C.5 Case 2, eq. (44) | f, He, Hs, φe → Hmin, Ag |
| `14_scintillation_attD.csv` | Att. D eqs. (48)–(50), P.618-12 §2.4.1 | f, θ, Nwet, Ga, p → σ, Ast |
| `15_ray_height_attE.csv` | Att. E eqs. (51)–(59) | Ht, θ, Dc → Hr |
| `16_ducting_ldtb.csv` | §2.6 eqs. (13)–(13e) | f, p, φ, θhoz, dhoz, Ld → Ldtb |
| `17_end_to_end_eq14.csv` | §3.1 eq. (14) | scenario → Dts, θ0, θ, Lbfs, Ag, Abs, As, Ldtb, Lb |
| `18_precip_scatter_attF.csv` | Att. F eqs. (60a)–(67) | P_int, f, r1, r2, Ptx, G, dtx, drx, γg, cylinder, θ, θpath, R, τ → L1, L2, S, Pin, γr, Pout, Pscat, Peirps, Fnis, Ptxs, decision |
| `19_fresnel_knife_edge.csv` | §2.6 eqs. (12a), (12b); P.526-16 eq. (31) | h, d1, f → ν, R1, J(ν) |
| `20_finite_width_screen.csv` | §2.6; P.526-16 §5.1 eqs. (71), (72) | νtop, νleft, νright → Jmin, Jav |
| `21_terrain_ldtb.csv` | §2.6, Att. E, eq. (13); P.452-18 (137), (141); single knife-edge method | f, p, φ, Ht, θ, ae, ridge profile → νmax, dmax, hmax, Ld, θhoz, dhoz, Ldtb |
| `22_terrain_slant_p526.csv` | §2.6, Att. E, eq. (13); P.526-16 §4.6 slant-path method | f, p, φ, Ht, θ, ae, two-ridge profile → ν, h, d1, J, Le, kb, Ld, θhoz, dhoz, Ldtb |
| `23_bel_p2109.csv` | §2.8; P.2109-2 eqs. (1)–(10) | f, p, building type, θ → Lbe |
| `24_multiple_entry_eq15.csv` | §3.2 eq. (15) | scenario, p, Axp, clutter/hs/hm/pLc, indoor/building/pLbe → Dts, θ0, θ, Lbfs, Ag, Abs, Lc, Lbe, Ldtb, Lb |
| `25_clutter_p2108_2_draft.csv` | §2.7; draft P.2108-2 §3.3 (Doc. 3/41(Rev.2)) | f, θ, p, h, hm → Lc |
| `26_p836_rho.csv` | §2.3; P.836-6 Annex 1 (annual), P.1144 | lat, lon, alt, p → ρ |
| `27_gas_ag_p1.csv` | §3.1 Ag(p1), §2.3, §4, Att. C | scenario, p1, ρ0 given (mode 1) or from P.836 (mode 2) → q, ρs, ρ0, θ, Ag |
| `28_nwet_p453.csv` | Att. D; P.453-14 §2.2 (annual maps), P.1144 | lat, lon, p → Nwet |
| `29_end_to_end_scint.csv` | §3.1 eq. (14) with As(p2), Att. D, P.453-14 | scenario, p2, Ga → θ0, Nwet, As, Lb |

Every file has an `id` column (1-based) and a `source` column; `17_end_to_end_eq14.csv` also has
`scenario` and `direction` (`E2s`/`s2E`).

## Conventions

* Units: frequency GHz, heights and distances km, angles degrees, losses dB; `log` = log10.
  Earth radius Re = 6371 km.
* **Attachment A.** θ0 is the elevation *above the horizontal*, i.e. arctan(Z2/Gts), and the
  bearing before step 8 is arctan(Y2/X2). P.619-6 defines atan2(x, y) = arctan(x/y), which,
  applied literally to eqs. (23)–(24), returns the zenith angle and a wrong bearing; the
  equations are only consistent with atan2(x, y) = arctan(y/x). δ is wrapped to (−180°, 180°].
  ψ is east of true north in [0°, 360°) and NaN for a vertical path.
* **Atmosphere.** P.835-6 mean annual global reference atmosphere. The mixing-ratio floor
  e/P ≥ 2·10⁻⁶ is applied to **both e and ρ** (as in the official P.676-13 slant-path sheets).
  P is total pressure; γ uses the dry pressure P − e. n from P.453 eqs. (1), (6). The 86 km
  boundary (where P.835-6 T is discontinuous by 0.08 K) is deliberately not tested.
* **Layering.** Layers per P.676 Annex 1 eq. (16) between He and min(Hs, 100 km); T, P, ρ, n
  and γ are evaluated at layer mid-points; ne is n of the first layer; ns is n of the top layer if
  Hs ≤ 100 km and 1 (free space) above. The Step 2/3 beam-width and ducting checks are not part
  of these vectors (antennas are taken to point along the ray).
* **Case 2.** Hmin is the root of eq. (45); Ag = Ag(Hmin→He, 0°) + Ag(Hmin→Hs, 0°). Rays for
  which Hmin < 0 (intersecting the Earth) are not included.
* **Beam spreading.** P.834 §5 formula, reported as a positive loss; in the end-to-end file
  Abs = 0 outside its validity range (θ0 > 10° or h > 5 km).
* **End-to-end (eq. 14).** Antennas point along the path. θ = apparent elevation from Att. B when
  Ht ≤ 3 km and −1° ≤ θ0 ≤ 10°, otherwise θ = θ0. Ag is computed for a ray launched at θ
  (Case 2 if θ < 0); s2E rows give the same Ag by reciprocity. Axp = 3 dB (§2.2 default),
  p2 = 50 % (As ≈ 0, tolerance 10⁻³ dB to absorb the Nwet-dependent residual), no terrain
  (Ldtb = 0). File 29 checks As with Nwet from the P.453 maps.
* **Faraday.** θF from eq. (4); AxF, AcF from eqs. (3a), (3b) *as printed*. Note that (3a) gives
  AxF = 0 for θF = 0 while (2a) with XPD = −20 log tan θF gives Ax → ∞ – the two are only
  consistent if (3a) and (3b) are swapped. Test cases keep 0 < θF < π/2.

* **Attachment F – suspected typo in eq. (60b).** P.619-6 prints
  L2 = max(L1·cos θ, 2·r1·**sin** θ). The vectors use **2·r1/sin θ**: it is the length of cylinder 2
  crossed by cylinder 1 (the same construction as L1 = 2·r2/sin θ in (60a)), and only with the
  division do (65a) and (65b) give the same scattered e.i.r.p. for equal beams (r1 = r2), as they
  must when the two cylinders are interchangeable. With the printed form, Peirps exceeds Pscat
  for r1 = r2 (by 1.5 dB at θ = 45°). Rows 3/4 and 7/8 of `18_precip_scatter_attF.csv` are such
  equal-beam pairs and must give identical Ptxs. Reported to ITU-R WP 3M.
* **Attachment F – other conventions.** θ is the scatter angle (10° ≤ θ ≤ 90°). γr from
  P.838-3 with the path inclination θpath and polarisation tilt τ given as inputs (0° horizontal,
  90° vertical, 45° circular). γg is an input and multiplies dtx and drx as written in (61) and
  (67). Ptxs is computed as written in (67), i.e. without a receive-antenna gain. The rain-cell
  validity condition (common-volume extent ≤ 6600·R^−0.08 m) is not evaluated. `full_calc` = 1
  when Pint − Ptxs < 20 dB; Pint is set 5 dB either side of that threshold.

* **§2.6 diffraction/ducting.**
  * ν and R1 from eqs. (12a), (12b) (h > 0: obstruction above the ray). J(ν) from P.526-16
    eq. (31), J = 0 for ν ≤ −0.78; this also covers “Ld = 0 if no obstacle enters the first Fresnel
    zone” (ν < −√2). The P.619-6 worked example for Fig. 10 quotes ν ≈ 4.9 for h = 60.3 m,
    d1 = 24 km, 30 GHz; eq. (12a) gives 5.51 (4.9 results from using R1 at 30 km instead of 24 km).
  * Finite-width screen: Jmin = −20 log(Σ 10^(−J(νk)/20)), Jav = −10 log(Σ 10^(−J(νk)/10)) over the
    top and the two side edges (P.526-16 eqs. (71), (72); (67), (68) in P.526-15).
  * Terrain: the ray is traced from Ht at the apparent elevation θ with Attachment E (1 km steps,
    linear interpolation between steps; eq. (59) evaluated directly for θ > 5°; points beyond the
    trace, where the ray is above 10 km, are ignored). This traced ray is the direct ray of P.526.
  * **Slant-path method (default; P.526-16 §4.6, file 22).** Vertical heights Δh relative to the
    ray at horizontal distance d are converted to the along-ray distance s = d/cos θ + Δh·sin θ and
    the perpendicular height n = Δh·cos θ. §4.6.1: if any n > 0, the Bullington apex is formed by
    the tangent ray from the antenna (α1 = max atan(n/s)) and the ray parallel to the direct ray
    tangent to the terrain (h = max n), d1 = h/tan α1; otherwise h and d1 are taken at the point of
    minimum Fresnel clearance (largest n/√s). ν from eq. (26) with d2 omitted (λ = c/f,
    c = 299 792 458 m/s). If ν ≤ −0.78, Ld = 0. §4.6.2: terrain raised by R1/π = √(λs/π),
    α1 and h of the raised terrain, Le from eq. (68) (satellite at great distance, h ≈ α2·d);
    Le = 0 if the raised terrain stays below the ray (cannot occur when ν > −0.78). §4.6.3:
    kb = max((Le − 26)/70, 0), Ld = kb·J(ν) + (1 − kb)·Le, used as printed (kb is not limited to
    1; P.526 states validity up to about 50 dB). Note that for Le < 26 dB, Ld = Le.
  * **Single knife-edge method (file 21).** The terrain is represented by a single knife edge at
    the profile point with the largest ν (vertical height relative to the ray, horizontal distance,
    eq. (12a)), Ld = J(νmax).
  * Horizon: θhoz and dhoz over the intermediate profile points (i = 2 … n−1) with
    θi = 1000·atan((hi − hts)/(1000·di) − di/(2ae)) mrad, as in the validated P.452 implementation;
    ae is an input (6371·4/3 km in the vectors; P.452 eqs. (5), (6a) with ΔN from the DN50 map in
    `tl_p619_single` when not given).
  * Synthetic profiles: di = i·step for i = 0 … round(dmax/step), hi = h0 + A·exp(−((di − dc)/w)²)
    m asl (file 22: sum of two such ridges with A1, dc1, w1 and A2, dc2, w2). The cases cover: unobstructed, obstructed with Ldtb = Ld (high horizon), partial ducting
    reduction, clamped at 0 (also |φ| > 70°), p ≥ β, and θ > 5°.
  * Eq. (13) is applied only to terrain obstruction; an Ld given for a specific obstruction (e.g. a
    building, finite-width screen) is used as Ldtb unchanged.

* **§3.2 multiple-entry, eq. (15).** Lb = Lbfs + Axp + Ag(p) + Abs + Lc(pLc) + Lbe(pLbe) + Ldtb(p);
    no scintillation. Lc by default from the draft revision P.2108-2 §3.3 (inputs hs, hm; 0.5–100 GHz),
    Q⁻¹ = inverse complementary normal; 0–90°. P.2108-1 §3.3 (`cl_p2108_1`, 'clutterModel',
    'P.2108-1'; P.619-6 cites P.2108-0) remains selectable in the code but is not part of the
    validation vectors. Lbe from P.2109-2 with
    p in % (P = p/100) and θ = path elevation at the façade. Both use the apparent elevation θ of
    the path. Lc = 0 for a station outside clutter, Lbe = 0 for an outdoor station. In a simulation
    pLc and pLbe are new uniform random values in (0, 100) for each loss; the vectors fix them.
    Aggregate power is the sum in linear units (`aggregate_power`).
  * Draft revision P.2108-2 §3.3 (Document 3/41(Rev.2); `cl_p2108_3`, selectable with
    'clutterModel', 'P.2108-2-draft', inputs h and hm): 419 values in file 25 are taken from the
    SG3 Clutter and BEL workbook values distributed with github.com/eeveetza/p2108 (branch dev):
    400 sampled rows and all rows at the boundary p = p_FcLoS, where eq. (16) applies eq. (15).
    The independent Python implementation reproduces all 282 800 workbook values within 5·10⁻⁹ dB.
    For k → 0 in eqs. (8a)/(10a) (mid-rise, h = 2.6 m) the limit (θ + Ce)/(90 + Ce) is used
    (3 rows, reference only; the draft does not state it).

* **Ag(p1) (files 26, 27).** P.836-6 annual surface water-vapour density exceeded for p %:
    four surrounding 1.125° grid points, each scaled to the station altitude with
    exp(−(alt − alt_i)/VSCH_i), grid-point altitudes from TOPO_0DOT5 (bicubic, P.1144 kernel
    a = −0.5; longitude wrapped to [0, 360), rows clamped at the poles), bilinear interpolation,
    linear interpolation in ln p between the 18 map percentages. Ag(p1) uses the density exceeded
    for q = 100 − p1 % (limited to 0.1 … 99 %) at the earth-based station, ρ0 = ρs·exp(He/2), in
    the P.835 mean annual global profile ρ(h) = ρ0·exp(−h/2) (mixing-ratio floor applied).
    Mode 1 rows give ρ0 directly and need no maps. The 32 ITU validation examples (ValEx sheet
    P.836-6 WV) are reproduced to 5·10⁻¹⁴ g/m³. The other reference values are computed with the
    map arrays distributed with ITU-Rpy (the ITU maps are not part of this repository); points
    north of 88.875° N are avoided because those arrays contain NaN cells in that row. Without the
    maps processed by `initiate_digital_maps.m`, the P.836 values are reported as not provided.

* **Nwet (files 28, 29).** P.453-14 annual maps of Nwet exceeded for p % (0.75° grid, rows
    −90 … 90°, columns −180 … 180°), bilinear interpolation of the four surrounding grid points
    (P.1144), longitude wrapped to [−180, 180) (the columns at ±180° are identical), linear
    interpolation in log p between the 18 map percentages (independent of the base of the log).
    The 8 ITU validation examples (ValEx sheet P.453-14 Nwet, p = 50 %) are reproduced to
    4·10⁻¹³; the other values are computed with the map arrays distributed with ITU-Rpy, which are
    identical to the P.453-14 files and to `matlab/private/DigitalMaps_NWET_Annual_xx.m`.
    File 29: Attachment D refers to the method of P.618-12 §2.4.1, in which Nwet is the long-term
    (median) value; as in P.618-14 and the ValEx sheet P.618-14 A_Scint, Nwet is the median (50 %)
    value at the earth-based station, whatever p2. As uses the free-space elevation θ0 (≥ 5° in all
    rows), Ga = 30 or 45 dBi, p2 = 0.01 … 99.99 %; the other terms are as in file 17 (MAGRA,
    p1 = 50 %, Axp = 3 dB, no terrain). One scenario gives the longitudes in [0, 360).

## Not covered (yet)

* §2.5.1 ionospheric scintillation (P.531), §2.10, §5.
* ΔN digital map (P.452): vectors take ae as an input instead.
* Beam-width and ducting checks of Att. C Steps 2–3.

## Status of P619.m (branch `dev`)

All 2791 values pass (with the P.836-6 maps; without them 162 values are not provided). (θF and σ, 96 values, are returned as additional optional outputs of `faraday_co_and_cross_polar_attenuation` and `tropospheric_scintillation` since this revision.)

Earlier failures and their fixes in `P619.m`:

| File | Deviation before | Cause | Fix |
|---|---|---|---|
| `09` ρ, e (h ≥ 25 km) | ≤ 2.2·10⁻⁵ g/m³ | P.835 mixing-ratio floor applied only to e inside `p453_n`, not to ρ; γ (P.676) used the unfloored ρ | floor applied to ρ in `p835_reference_atmosphere` (also reduces Ag deviations from 8·10⁻⁶ to 5·10⁻⁷ dB) |
| `16` Ldtb | ≤ 4.9·10⁻⁶ dB | coefficient 0.195 in eq. (13b) | 0.198 as in P.619-6 (and P.452) |
| `17` Abs, Lb (θ0 > 10°) | ≤ 0.022 dB | P.834 formula applied outside θ0 ≤ 10°, h ≤ 5 km | `beam_spreading_loss` returns 0 outside its validity range |
| `18` γr and everything after it | up to 8 % in γr | `test_precipitation_scatter` passed θpath in degrees to `specific_rain_attenuation_p838`, which expects radians | converted to radians |
| `18` L2 (5 of 12 rows) | – | eq. (60b) as printed (2·r1·sin θ) | 2·r1/sin θ, suspected typo (see *Conventions*) |
| (none) `bel_p2109` | could not run | `A` was never assigned (second assignment overwrote `B`); `cl_p2108_3`, `bel_p2109` lacked the `obj` argument | fixed |
| `29` Nwet, As (p2 ≠ 50 %) | – | `tl_p619_single` took Nwet exceeded for p2 % from the P.453 maps | median (50 %) Nwet; maps cached, any longitude range |
| `25` Lc (draft P.2108-2) | NaN in 16 of 282 800 workbook rows, 3 rows 0.002 dB | older copy of the model: no k → 0 limit; eq. (14) used at p = p_FcLoS (0/0 when p_LoS = p_FcLoS) | aligned with github.com/eeveetza/p2108 (dev) |
