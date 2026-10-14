"""
Independent reference implementation of selected parts of Recommendation
ITU-R P.619-6 (09/2025), used to generate the validation vectors in
../vectors/.

It was written from the text of the Recommendations (P.619-6, P.834,
P.835-6, P.453, P.618-12 §2.4.1, P.676 Annex 1) and NOT ported from
P619.m, so that agreement between the two is meaningful. Where a step is
iterative, a different numerical method is used on purpose (e.g. eq. (45)
is solved with a bracketing root finder instead of the Newton iteration
(46)-(47)).

Specific attenuation (P.676 Annex 1, eq. 1) is delegated to the
third-party ITU-Rpy package (itur.models.itu676.gamma0_exact /
gammaw_exact), so it is also independent of P619.m.

Conventions (see ../README.md for details):
  * Earth radius Re = 6371 km
  * angles in degrees, heights and distances in km, frequency in GHz
  * log = log10
"""

import numpy as np
from scipy.optimize import brentq

import itur.models.itu676 as itu676

RE = 6371.0

# ---------------------------------------------------------------------------
# §2.1  Free-space basic transmission loss
# ---------------------------------------------------------------------------

def lbfs(f_ghz, d_km):
    """Eq. (1)."""
    return 92.45 + 20.0 * np.log10(f_ghz * d_km)


# ---------------------------------------------------------------------------
# §2.2  Depolarisation
# ---------------------------------------------------------------------------

def xpd_attenuation(xpd_db):
    """Eqs. (2a), (2b): cross-polar Ax and co-polar Ac attenuation from XPD (dB)."""
    ax = 10.0 * np.log10(1.0 + 10.0 ** (0.1 * xpd_db))
    ac = 10.0 * np.log10(1.0 + 10.0 ** (-0.1 * xpd_db))
    return ax, ac


def faraday_rotation(f_ghz, bav_t, nt_m2):
    """Eq. (4): Faraday rotation angle (rad); f in GHz, Bav in T, NT in el/m^2."""
    return 2.36e-14 * bav_t * nt_m2 / f_ghz ** 2


def faraday_attenuation_literal(theta_f):
    """Eqs. (3a), (3b) as printed in P.619-6 (see README: possibly swapped)."""
    return -20.0 * np.log10(np.cos(theta_f)), -20.0 * np.log10(np.sin(theta_f))


def hydrometeor_depol_attenuation(xpd_db):
    """Eq. (6): Axq from the hydrometeor XPD."""
    return -20.0 * np.log10(np.cos(np.arctan(10.0 ** (-xpd_db / 20.0))))


# ---------------------------------------------------------------------------
# §2.4.2  Beam-spreading loss (Rec. ITU-R P.834 §5)
# ---------------------------------------------------------------------------

def beam_spreading_loss(theta0_deg, h_km):
    """Abs (dB, reported as a positive loss) for theta0 <= 10 deg, h <= 5 km."""
    t, h = theta0_deg, h_km
    num = 0.5411 + 0.07446 * t + h * (0.06272 + 0.0276 * t) + h ** 2 * 0.008288
    den = (1.728 + 0.5411 * t + 0.03723 * t ** 2
           + h * (0.1815 + 0.06272 * t + 0.0138 * t ** 2)
           + h ** 2 * (0.01727 + 0.008288 * t)) ** 2
    return abs(10.0 * np.log10(1.0 - num / den))


# ---------------------------------------------------------------------------
# Attachment A  Geometry of a straight-line Earth-space path
# ---------------------------------------------------------------------------

def straight_line_geometry(hs, ht, lat_s, lat_t, delta):
    """Eqs. (18)-(24). Returns Dts (km), theta0 (deg), psi (deg E of N, [0,360)).

    atan2 is used in the sense that makes theta0 the angle above the
    horizontal: theta0 = arctan(Z2/Gts), psi_south = arctan(Y2/X2).
    """
    rs, rt = RE + hs, RE + ht
    lat_s, lat_t, delta = map(np.radians, (lat_s, lat_t, delta))
    x1 = rs * np.cos(lat_s) * np.cos(delta)
    y1 = rs * np.cos(lat_s) * np.sin(delta)
    z1 = rs * np.sin(lat_s)
    x2 = x1 * np.sin(lat_t) - z1 * np.cos(lat_t)
    y2 = y1
    z2 = z1 * np.sin(lat_t) + x1 * np.cos(lat_t) - rt
    dts = np.sqrt(x2 ** 2 + y2 ** 2 + z2 ** 2)
    gts = np.sqrt(x2 ** 2 + y2 ** 2)
    theta0 = np.degrees(np.arctan2(z2, gts))
    psi = np.mod(180.0 - np.degrees(np.arctan2(y2, x2)), 360.0)
    return dts, theta0, psi


def look_angles_ecef(hs, ht, lat_s, lon_s, lat_t, lon_t):
    """Independent cross-check of Attachment A using spherical-Earth ECEF
    vectors and a local ENU frame (does not use eqs. 19-24)."""
    def ecef(r, lat, lon):
        lat, lon = np.radians(lat), np.radians(lon)
        return r * np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    s = ecef(RE + hs, lat_s, lon_s)
    t = ecef(RE + ht, lat_t, lon_t)
    v = s - t
    la, lo = np.radians(lat_t), np.radians(lon_t)
    east = np.array([-np.sin(lo), np.cos(lo), 0.0])
    north = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    up = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    e, n, u = v @ east, v @ north, v @ up
    return np.linalg.norm(v), np.degrees(np.arctan2(u, np.hypot(e, n))), np.mod(np.degrees(np.arctan2(e, n)), 360.0)


# ---------------------------------------------------------------------------
# Attachment B  Apparent <-> free-space elevation
# ---------------------------------------------------------------------------

def fs_to_apparent(theta0, ht):
    """Eqs. (25)-(26), Ht <= 3 km, -1 <= theta0 <= 10 deg."""
    t1 = 1.728 + 0.5411 * theta0 + 0.03723 * theta0 ** 2
    t2 = 0.1815 + 0.06272 * theta0 + 0.01380 * theta0 ** 2
    t3 = 0.01727 + 0.008288 * theta0
    return theta0 + 1.0 / (t1 + ht * t2 + ht ** 2 * t3)


def apparent_to_fs(theta, ht):
    """Eqs. (27)-(28)."""
    t1 = 1.314 + 0.6437 * theta + 0.02869 * theta ** 2
    t2 = 0.2305 + 0.09428 * theta + 0.01096 * theta ** 2
    t3 = 0.008583
    return theta - 1.0 / (t1 + ht * t2 + ht ** 2 * t3)


# ---------------------------------------------------------------------------
# Rec. ITU-R P.835-6 mean annual global reference atmosphere (MAGRA)
# ---------------------------------------------------------------------------

def _geopotential(h):
    return 6356.766 * h / (6356.766 + h)


def magra_temperature(h):
    """Temperature (K) at geometric height h (km), 0 <= h <= 100."""
    hp = _geopotential(h)
    if h > 86.0 or np.isclose(h, 86.0):
        if h <= 91.0:
            return 186.8673
        return 263.1905 - 76.3232 * np.sqrt(1.0 - ((h - 91.0) / 19.9429) ** 2)
    if hp <= 11:
        return 288.15 - 6.5 * hp
    if hp <= 20:
        return 216.65
    if hp <= 32:
        return 216.65 + (hp - 20.0)
    if hp <= 47:
        return 228.65 + 2.8 * (hp - 32.0)
    if hp <= 51:
        return 270.65
    if hp <= 71:
        return 270.65 - 2.8 * (hp - 51.0)
    return 214.65 - 2.0 * (hp - 71.0)


def magra_pressure(h):
    """Total pressure (hPa) at geometric height h (km)."""
    hp = _geopotential(h)
    if h > 86.0 or np.isclose(h, 86.0):
        a = (95.571899, -4.011801, 6.424731e-2, -4.789660e-4, 1.340543e-6)
        return np.exp(sum(ai * h ** i for i, ai in enumerate(a)))
    if hp <= 11:
        return 1013.25 * (288.15 / (288.15 - 6.5 * hp)) ** (-34.1632 / 6.5)
    if hp <= 20:
        return 226.3226 * np.exp(-34.1632 * (hp - 11.0) / 216.65)
    if hp <= 32:
        return 54.74980 * (216.65 / (216.65 + (hp - 20.0))) ** 34.1632
    if hp <= 47:
        return 8.680422 * (228.65 / (228.65 + 2.8 * (hp - 32.0))) ** (34.1632 / 2.8)
    if hp <= 51:
        return 1.109106 * np.exp(-34.1632 * (hp - 47.0) / 270.65)
    if hp <= 71:
        return 0.6694167 * (270.65 / (270.65 - 2.8 * (hp - 51.0))) ** (-34.1632 / 2.8)
    return 0.03956649 * (214.65 / (214.65 - 2.0 * (hp - 71.0))) ** (-34.1632 / 2.0)


def magra(h, rho0=7.5):
    """T (K), total P (hPa), rho (g/m^3), e (hPa) with the P.835 mixing-ratio
    floor e/P >= 2e-6 applied to both e and rho."""
    t = magra_temperature(h)
    p = magra_pressure(h)
    rho = rho0 * np.exp(-h / 2.0)
    e = rho * t / 216.7
    if e < 2e-6 * p:
        e = 2e-6 * p
        rho = e * 216.7 / t
    return t, p, rho, e


def refractive_index(t, p, e):
    """Rec. ITU-R P.453 eqs. (1), (6): n from total pressure P, e, T."""
    pd = p - e
    n_units = 77.6 * pd / t + 72.0 * e / t + 3.75e5 * e / t ** 2
    return 1.0 + 1e-6 * n_units


def magra_n(h):
    t, p, _, e = magra(h)
    return refractive_index(t, p, e)


# ---------------------------------------------------------------------------
# Attachment C  Attenuation by atmospheric gases
# ---------------------------------------------------------------------------

def layer_boundaries(h_lower, h_upper):
    """Rec. ITU-R P.676 Annex 1 eqs. (14)-(16): boundaries of the exponentially
    thickening layers between h_lower and h_upper (0 <= h_lower < h_upper <= 100)."""
    i_l = int(np.floor(100.0 * np.log(1e4 * h_lower * (np.exp(0.01) - 1.0) + 1.0) + 1.0))
    i_u = int(np.ceil(100.0 * np.log(1e4 * h_upper * (np.exp(0.01) - 1.0) + 1.0) + 1.0))
    m = (np.exp(0.02) - np.exp(0.01)) / (np.exp(i_u / 100.0) - np.exp(i_l / 100.0)) * (h_upper - h_lower)
    i = np.arange(i_l, i_u + 1)
    return h_lower + m * (np.exp(0.01 * (i - 1)) - np.exp(0.01 * (i_l - 1))) / (np.exp(0.01) - 1.0)


def _layer_profiles(hb, rho0=7.5):
    hm = 0.5 * (hb[:-1] + hb[1:])
    prof = np.array([magra(h, rho0) for h in hm])
    t, p, rho, e = prof.T
    n = refractive_index(t, p, e)
    return t, p, rho, e, n


def _gamma(f, t, p, rho, e):
    """Specific attenuation (dB/km), P.676 Annex 1 eq. (1), via ITU-Rpy."""
    pd = p - e
    g0 = np.atleast_1d(itu676.gamma0_exact(f, pd, rho, t).value)
    gw = np.atleast_1d(itu676.gammaw_exact(f, pd, rho, t).value)
    return g0 + gw


def specific_attenuation(f, pd, rho, t):
    """γo, γw (dB/km) at dry pressure pd (hPa), rho (g/m^3), T (K)."""
    return (float(itu676.gamma0_exact(f, pd, rho, t).value),
            float(itu676.gammaw_exact(f, pd, rho, t).value))


def ag_ascending(f, he, hs, phi_e, rho0=7.5):
    """C.5 Case 1 (phi_e >= 0), eqs. (39)-(40), layered MAGRA atmosphere from
    he to min(hs, 100) km; above 100 km the space station is in free space.
    The beam-width / ducting checks of Steps 2-3 are not applied."""
    hb = layer_boundaries(he, min(hs, 100.0))
    t, p, rho, e, n = _layer_profiles(hb, rho0)
    r = RE + hb
    ne = n[0]
    re = RE + he
    c2 = (ne / n * re * np.cos(np.radians(phi_e))) ** 2
    l = np.sqrt(np.maximum(r[1:] ** 2 - c2, 0)) - np.sqrt(np.maximum(r[:-1] ** 2 - c2, 0))
    return float(np.sum(_gamma(f, t, p, rho, e) * l))


def ag_descending(f, he, hs, phi_s):
    """C.4, eqs. (33)-(34): ray launched from the space station at phi_s (<0)."""
    hb = layer_boundaries(he, min(hs, 100.0))
    t, p, rho, e, n = _layer_profiles(hb)
    r = RE + hb
    ns = 1.0 if hs > hb[-1] + 1e-9 else n[-1]
    rs = RE + hs
    c2 = (ns / n * rs * np.cos(np.radians(phi_s))) ** 2
    l = np.sqrt(np.maximum(r[1:] ** 2 - c2, 0)) - np.sqrt(np.maximum(r[:-1] ** 2 - c2, 0))
    return float(np.sum(_gamma(f, t, p, rho, e) * l))


def space_launch_angle(he, hs, theta_e):
    """Elevation (negative, deg) at which a ray must leave the space station to
    arrive at the earth station with elevation theta_e (Snell, eq. 35)."""
    hb = layer_boundaries(he, min(hs, 100.0))
    _, _, _, _, n = _layer_profiles(hb)
    ns = 1.0 if hs > hb[-1] + 1e-9 else n[-1]
    return -np.degrees(np.arccos(min(1.0, (RE + he) * n[0] / ((RE + hs) * ns) * np.cos(np.radians(theta_e)))))


def hmin_eq45(he, phi_e):
    """Root of eq. (45) by Brent's method (independent of the Newton iteration)."""
    c = (RE + he) * magra_n(he) * np.cos(np.radians(phi_e))
    g = lambda h: (RE + h) * magra_n(h) - c
    return brentq(g, 0.0, he, xtol=1e-13, rtol=1e-15, maxiter=500)


def ag_negative_elevation(f, he, hs, phi_e):
    """C.5 Case 2, eq. (44): He -> Hmin plus Hmin -> space, both at 0 deg at Hmin."""
    hmin = hmin_eq45(he, phi_e)
    return hmin, ag_ascending(f, hmin, he, 0.0) + ag_ascending(f, hmin, hs, 0.0)


# ---------------------------------------------------------------------------
# Attachment D  Tropospheric scintillation
# ---------------------------------------------------------------------------

def scintillation_sigma(f, theta, nwet, ga):
    """σ (dB) per Rec. ITU-R P.618-12 §2.4.1 Steps 3-7, with the effective
    antenna diameter from eq. (48). theta = elevation (deg), f in GHz."""
    sigma_ref = 3.6e-3 + 1e-4 * nwet
    s = np.sin(np.radians(theta))
    l = 2.0 * 1000.0 / (np.sqrt(s ** 2 + 2.35e-4) + s)
    deff = 0.3 * 10.0 ** (0.05 * ga) / (np.pi * f)
    x = 1.22 * deff ** 2 * f / l
    if x >= 7.0:
        return 0.0
    g = np.sqrt(3.86 * (x ** 2 + 1.0) ** (11.0 / 12.0) * np.sin(11.0 / 6.0 * np.arctan(1.0 / x)) - 7.08 * x ** (5.0 / 6.0))
    return sigma_ref * f ** (7.0 / 12.0) * g / s ** 1.2


def scintillation_ast(f, p, theta, nwet, ga):
    """Eqs. (49a), (49b), (50): Ast(p) (dB), negative = enhancement."""
    sigma = scintillation_sigma(f, theta, nwet, ga)
    if p <= 50.0:
        lp = np.log10(p)
        return -sigma * (2.672 - 1.258 * lp - 0.0835 * lp ** 2 - 0.0597 * lp ** 3)
    lq = np.log10(100.0 - p)
    return sigma * (3.0 - 1.71 * lq + 0.072 * lq ** 2 - 0.061 * lq ** 3)


# ---------------------------------------------------------------------------
# Attachment E  Beam clearance
# ---------------------------------------------------------------------------

def ray_height(ht, theta, dc):
    """Ray altitude Hr (km) at horizontal distance dc (km, integer for theta <= 5)."""
    if theta <= 5.0:
        hr, eps = ht, np.radians(theta)
        for _ in range(int(round(dc))):
            d_eps = 1.0 * (1.0 / RE - 4.28715e-5 * np.exp(-hr / 7.348))
            hr = hr + 1.0 * eps
            eps = eps + d_eps
        return hr
    return ht + dc * np.tan(np.radians(theta)) + dc ** 2 / (2.0 * RE)


# ---------------------------------------------------------------------------
# §2.6  Diffraction/ducting loss
# ---------------------------------------------------------------------------

def ldtb(f, p, lat, theta_hoz, d_hoz, ld):
    """Eqs. (13)-(13e)."""
    beta = 10.0 ** (-0.015 * abs(lat) + 1.67) if abs(lat) <= 70.0 else 4.17
    dn = 600.0 / (1.0 + f)
    lb = np.log10(beta)
    gam = 1.076 / (2.0058 - lb) ** 1.012 * np.exp(-(9.51 - 4.8 * lb + 0.198 * lb ** 2) * 1e-6 * dn ** 1.13)
    if p >= beta:
        return ld
    a = (1.2 + 3.7e-3 * dn) * np.log10(p / beta) + 12.0 * ((p / beta) ** gam - 1.0)
    th = theta_hoz - 0.1 * d_hoz
    ads = 20.0 * np.log10(1.0 + 0.361 * th * np.sqrt(f * d_hoz)) + 0.264 * th * f ** (1.0 / 3.0) if th > 0 else 0.0
    return max(ld + min(a + ads, 0.0), 0.0)


# ---------------------------------------------------------------------------
# Attachment F  Test for whether a precipitation-scatter calculation is needed
# ---------------------------------------------------------------------------

def precipitation_scatter_test(p_int, f, r1, r2, p_tx, g, d_tx, d_rx, gamma_g, cyl, theta, theta_path, r_rain, tau):
    """Eqs. (60a)-(67). r1 <= r2 in m, d_tx/d_rx in km, theta = scatter angle (deg),
    theta_path = path inclination for P.838 (deg), tau = polarisation tilt (deg).

    NOTE eq. (60b): P.619-6 prints 2*r1*sin(theta); 2*r1/sin(theta) is used here
    (suspected typo, see ../README.md)."""
    import itur.models.itu838 as itu838
    s, c = np.sin(np.radians(theta)), np.cos(np.radians(theta))
    l1 = 2.0 * r2 / s                                           # (60a)
    l2 = max(l1 * c, 2.0 * r1 / s)                              # (60b), division
    p_eirp = p_tx + g                                           # (61a)
    s_pfd = p_eirp - 20.0 * np.log10(d_tx) - gamma_g * d_tx - 71.0   # (61)
    r_beam = r1 if cyl == 1 else r2
    p_in = s_pfd + 10.0 * np.log10(np.pi * r_beam ** 2)         # (62a), (62b)
    gamma_r = float(itu838.rain_specific_attenuation(r_rain, f, theta_path, tau).value)  # (63c)
    p_out = p_in - 0.001 * gamma_r * (l1 if cyl == 1 else l2)  # (63a), (63b)
    p_scat = 10.0 * np.log10(10.0 ** (0.1 * p_in) - 10.0 ** (0.1 * p_out))  # (64)
    if cyl == 1:
        p_eirps = p_scat                                        # (65a)
    else:
        p_eirps = p_scat - 10.0 * np.log10(r2 ** 2 * l2 / (r1 ** 2 * l1))  # (65b)
    if f > 10.0:
        f_nis = 1e-3 * r_rain ** 0.4 * c * (2.0 * (f - 10.0) ** 1.6 - 2.5 * (f - 10.0) ** 1.7)  # (66)
    else:
        f_nis = 0.0
    p_txs = p_eirps + f_nis - 20.0 * np.log10(d_rx * f) - gamma_g * d_rx - 92.4  # (67)
    return dict(L1_m=l1, L2_m=l2, S_dBWm2=s_pfd, Pin_dBW=p_in, gamma_r_dBkm=gamma_r, Pout_dBW=p_out,
                Pscat_dBW=p_scat, Peirps_dBW=p_eirps, Fnis_dB=f_nis, Ptxs_dBW=p_txs,
                full_calc=1.0 if p_int - p_txs < 20.0 else 0.0)


# ---------------------------------------------------------------------------
# §2.6  Diffraction/ducting loss due to specific obstruction
# ---------------------------------------------------------------------------

def fresnel_parameters(h_m, d1_km, f):
    """Eqs. (12a), (12b): diffraction parameter nu and first Fresnel radius R1 (m)."""
    return 0.08168 * h_m * np.sqrt(f / d1_km), 17.314 * np.sqrt(d1_km / f)


def knife_edge_loss(nu):
    """Rec. ITU-R P.526-16 eq. (31); 0 for nu <= -0.78."""
    if nu <= -0.78:
        return 0.0
    return 6.9 + 20.0 * np.log10(np.sqrt((nu - 0.1) ** 2 + 1.0) + nu - 0.1)


def finite_width_screen_loss(nu_top, nu_left, nu_right):
    """Rec. ITU-R P.526-16 §5.1 eqs. (71), (72): (Jmin, Jav)."""
    inv = [10.0 ** (-knife_edge_loss(n) / 20.0) for n in (nu_top, nu_left, nu_right)]
    return -20.0 * np.log10(sum(inv)), -10.0 * np.log10(sum(x * x for x in inv))


def ray_trace(ht, theta):
    """Attachment E, eqs. (51)-(58): (Dc, Hr) in km at 1 km steps for theta <= 5 deg,
    up to and including the first point with Hr >= 10 km (at most 2000 points)."""
    dc, hr, eps = [0.0], [ht], np.radians(theta)
    while len(hr) < 2000 and hr[-1] < 10.0:
        d_eps = 1.0 / RE - 4.28715e-5 * np.exp(-hr[-1] / 7.348)
        hr.append(hr[-1] + eps)
        eps += d_eps
        dc.append(dc[-1] + 1.0)
    return np.array(dc), np.array(hr)


def ray_height_at(ht, theta, d):
    """Ray height asl (km) at profile distances d (km); NaN beyond the trace."""
    d = np.asarray(d, dtype=float)
    if theta > 5.0:
        return ht + d * np.tan(np.radians(theta)) + d ** 2 / (2.0 * RE)   # (59)
    dc, hr = ray_trace(ht, theta)
    out = np.interp(d, dc, hr)
    out[d > dc[-1]] = np.nan
    return out


def ridge_profile(h0, a, dc, w, dmax, step):
    """Synthetic profile used by the vectors: d_i = i*step, i = 0..round(dmax/step),
    h_i = h0 + a*exp(-((d_i - dc)/w)^2) (m asl)."""
    n = int(round(dmax / step))
    d = np.arange(n + 1) * step
    return d, h0 + a * np.exp(-((d - dc) / w) ** 2)


def terrain_obstruction(f, ht, theta, d, h):
    """Ld with the single knife edge at max nu along the ray, and its parameters."""
    d, h = np.asarray(d, float), np.asarray(h, float)
    sel = d > 0
    hray = ray_height_at(ht, theta, d[sel]) * 1000.0
    ok = ~np.isnan(hray)
    if not ok.any():
        return 0.0, dict(nu_max=-np.inf, d_max=np.nan, h_max=np.nan)
    ds, hob = d[sel][ok], h[sel][ok] - hray[ok]
    nu = 0.08168 * hob * np.sqrt(f / ds)
    i = int(np.argmax(nu))
    return knife_edge_loss(nu[i]), dict(nu_max=nu[i], d_max=ds[i], h_max=hob[i])


def horizon_p452(d, h, hts_m, ae):
    """Horizon elevation angle (mrad) and distance (km) over the intermediate profile
    points i = 2..n-1 (P.452-18 eqs. (137), (141)), in the atan form used by the
    validated P.452 implementation."""
    d, h = np.asarray(d, float)[1:-1], np.asarray(h, float)[1:-1]
    th = 1000.0 * np.arctan((h - hts_m) / (1000.0 * d) - d / (2.0 * ae))
    i = int(np.argmax(th))
    return th[i], d[i]


def diffraction_ducting(f, p, lat, ht, theta, d, h, ae):
    ld, o = terrain_obstruction(f, ht, theta, d, h)
    th_hoz, d_hoz = horizon_p452(d, h, ht * 1000.0, ae)
    o.update(Ld=ld, th_hoz=th_hoz, d_hoz=d_hoz, Ldtb=ldtb(f, p, lat, th_hoz, d_hoz, ld))
    return o


def slant_path_p526(f, s_km, n_m):
    """Rec. ITU-R P.526-16 §4.6 (general Earth-space slant path).
    s: distances along the direct ray (km, > 0); n: heights perpendicular to it (m).
    Returns dict with nu, h (m), d1 (km), J, Le, kb, Leb (dB)."""
    lam = 299792458.0 / (f * 1e9)
    s = np.asarray(s_km, float) * 1000.0
    n = np.asarray(n_m, float)
    if n.max() > 0.0:                                   # 4.6.1, obstructed: Bullington apex
        a1 = np.max(np.arctan2(n, s))
        h = n.max()
        d1 = h / np.tan(a1)
    else:                                               # line of sight: minimum Fresnel clearance
        i = int(np.argmax(n / np.sqrt(s)))
        h, d1 = n[i], s[i]
    nu = h * np.sqrt(2.0 / (lam * d1))                  # (26), d2 -> infinity
    res = dict(nu=nu, h=h, d1=d1 / 1000.0, J=0.0, Le=0.0, kb=0.0, Leb=0.0)
    if nu <= -0.78:
        return res
    J = knife_edge_loss(nu)                             # (31)
    ne = n + np.sqrt(lam * s / np.pi)                   # 4.6.2, terrain raised by R_1/pi
    he = ne.max()
    Le = 0.0
    if he > 0.0:
        a1e = np.max(np.arctan2(ne, s))
        x = np.pi * a1e * he / lam
        Le = -40.0 * np.log10(np.sqrt(x + 1.0) - np.sqrt(x))   # (68)
    kb = max((Le - 26.0) / 70.0, 0.0)                   # (70)
    res.update(J=J, Le=Le, kb=kb, Leb=kb * J + (1.0 - kb) * Le)  # (69)
    return res


def ridges_profile(h0, ridges, dmax, step):
    """d_i = i*step; h_i = h0 + sum_k A_k exp(-((d_i - dc_k)/w_k)^2), ridges = [(A, dc, w), ...]."""
    n = int(round(dmax / step))
    d = np.arange(n + 1) * step
    h = np.full_like(d, float(h0))
    for a, dc, w in ridges:
        h = h + a * np.exp(-((d - dc) / w) ** 2)
    return d, h


def terrain_slant(f, ht, theta, d, h):
    """Terrain heights relative to the Attachment E ray, converted to along-ray
    distance s and perpendicular height n with the path elevation theta, then §4.6."""
    d, h = np.asarray(d, float), np.asarray(h, float)
    sel = d > 0
    hray = ray_height_at(ht, theta, d[sel]) * 1000.0
    ok = ~np.isnan(hray)
    if not ok.any():
        return dict(nu=-np.inf, h=np.nan, d1=np.nan, J=0.0, Le=0.0, kb=0.0, Leb=0.0)
    dv, dh = d[sel][ok], (h[sel][ok] - hray[ok])
    c, sn = np.cos(np.radians(theta)), np.sin(np.radians(theta))
    return slant_path_p526(f, dv / c + dh / 1000.0 * sn, dh * c)


def diffraction_ducting_slant(f, p, lat, ht, theta, d, h, ae):
    o = terrain_slant(f, ht, theta, d, h)
    th_hoz, d_hoz = horizon_p452(d, h, ht * 1000.0, ae)
    o.update(Ld=o['Leb'], th_hoz=th_hoz, d_hoz=d_hoz, Ldtb=ldtb(f, p, lat, th_hoz, d_hoz, o['Leb']))
    return o


# ---------------------------------------------------------------------------
# §3.2  Multiple-entry interference: building entry loss (P.2109-2)
# ---------------------------------------------------------------------------

BEL_COEF = {1: dict(r=12.64, s=3.72, t=0.96, u=9.6, v=2.0, w=9.1, x=-3.0, y=4.5, z=-2.0),
            2: dict(r=28.19, s=-3.00, t=8.48, u=13.5, v=3.8, w=27.8, x=-2.9, y=9.4, z=-2.1)}


def bel_p2109(f, p, building, theta):
    """Rec. ITU-R P.2109-2 eqs. (1)-(10): building entry loss (dB) not exceeded for
    p % (0..100) of locations; building 1 = traditional, 2 = thermally efficient."""
    from scipy.stats import norm
    c = BEL_COEF[building]
    lf = np.log10(f)
    mu1 = c['r'] + c['s'] * lf + c['t'] * lf ** 2 + 0.212 * abs(theta)
    mu2 = c['w'] + c['x'] * lf
    s1 = c['u'] + c['v'] * lf
    s2 = c['y'] + c['z'] * lf
    fi = norm.ppf(p / 100.0)
    a, b, cc = fi * s1 + mu1, fi * s2 + mu2, -3.0
    return 10.0 * np.log10(10.0 ** (0.1 * a) + 10.0 ** (0.1 * b) + 10.0 ** (0.1 * cc))


# Draft revision P.2108-2 §3.3 (Document 3/41(Rev.2), 2025): Earth-space clutter loss as a
# function of f, theta, p, ground station height h (m) and median clutter height hm (m)
P2108_2_T7 = {'low': (4.9, 6.7, 2.6, 0.19, 0.0, 1.4, 74.0), 'mid': (-2.6, 6.6, 2.0, 0.42, -6.7, 0.15, 97.0),
              'high': (2.4, 7.0, 1.0, 0.19, -2.7, 0.15, 98.0)}           # ak bk ck aC bC aV bV
P2108_2_T8 = {'low': (6.0, 0.07, 0.15, 5.4, -0.3, 3.2, 0.07, -27.0, 1.6, -17.0),
              'mid': (3.6, 0.05, 0.17, 13.0, -0.2, 3.7, 0.05, -41.0, 1.0, -21.0),
              'high': (5.0, 0.003, 0.17, 32.6, 0.012, -23.9, -0.07, -41.0, 1.0, -18.0)}
#              a'k  b'k   a'c   b'c1  b'c2  b'c3  b'c4  c'c   a'V  b'V


def _p_state(k, theta, ce, vmax):
    """Eqs. (7)/(9) shape; for k -> 0 the limit (theta + Ce)/(90 + Ce) is used."""
    if abs(k) < 1e-6:
        return max(0.0, vmax * (theta + ce) / (90.0 + ce))
    return max(0.0, vmax * (1.0 - np.exp(-k * (theta + ce) / 90.0)) / (1.0 - np.exp(-k * (90.0 + ce) / 90.0)))


def clutter_p2108_2(f, theta, p, h, hm):
    from scipy.stats import norm
    sc = 'low' if hm <= 8 else ('mid' if hm <= 20 else 'high')
    ak, bk, ck, ac, bc, av, bv = P2108_2_T7[sc]
    k = ((h + ak) / bk) ** ck                                   # (8a)
    p_los = _p_state(k, theta, ac * h + bc, min(av * h + bv, 100.0))        # (7), (8b), (8c)
    akp, bkp, acp, bc1, bc2, bc3, bc4, ccp, avp, bvp = P2108_2_T8[sc]
    kp = akp * np.exp(bkp * h)                                  # (10a)
    cep = (f * 1e9) ** acp + bc1 * np.exp(bc2 * h) + bc3 * np.exp(bc4 * h) + ccp   # (10b)
    vmp = min(avp * h + bvp, 0.0) * f ** -0.55 + 100.0          # (10c)
    p_fc = p_los * _p_state(kp, theta, cep, vmp) / 100.0        # (9), (11)
    if p > p_los:                                               # NLoS, (12)-(13)
        pp = (p - p_los) / (100.0 - p_los)
        mu = 8.54 + 17.57 * np.log(1.0 + (90.0 - theta) / 90.0) + f ** 0.63
        sg = 0.056 + 6.32 * np.log(1.0 + (90.0 - theta) / 90.0) + f ** 0.19
        return max(mu + sg * norm.isf(1.0 - pp), 6.0)           # Q^-1(1 - P)
    if p <= p_fc:                                               # FcLoS, (15)
        return -3.542 * np.exp(-155.1 * p / p_fc) - 1.06 * np.exp(-6.342 * p / p_fc)
    return 6.0 * (p - p_fc) / (p_los - p_fc)                    # FoLoS, (14)


# ---------------------------------------------------------------------------
# Rec. ITU-R P.836-6 Annex 1 §1: surface water-vapour density exceeded for p %
# of an average year (annual maps), with Rec. ITU-R P.1144 interpolation.
# Map arrays: rho[p], vsch[p] on the 1.125 deg grid (161 x 321, lat 90 -> -90,
# lon 0 -> 360) and topo on the 0.5 deg grid with a border of one row/column
# (363 x 723, lat 90.5 -> -90.5, lon -0.5 -> 360.5); heights and scale heights in km.
# ---------------------------------------------------------------------------

P836_P = (0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0, 60.0, 70.0, 80.0, 90.0, 95.0, 99.0)


def _k1144(d):
    """Rec. ITU-R P.1144 bicubic kernel (a = -0.5)."""
    d = abs(d)
    if d <= 1.0:
        return 1.5 * d ** 3 - 2.5 * d ** 2 + 1.0
    if d <= 2.0:
        return -0.5 * d ** 3 + 2.5 * d ** 2 - 4.0 * d + 2.0
    return 0.0


def _bicubic(grid, r, c):
    """Bicubic interpolation at fractional (row, column) indices."""
    R, C = int(np.floor(r)), int(np.floor(c))
    nr, nc = grid.shape
    val = 0.0
    for i in range(R - 1, R + 3):
        ii = min(max(i, 0), nr - 1)                    # clamp at the poles
        row = sum(grid[ii, min(max(j, 0), nc - 1)] * _k1144(c - j) for j in range(C - 1, C + 3))
        val += row * _k1144(r - i)
    return val


def p836_rho(maps, lat, lon, p, alt_km):
    """Surface water-vapour density (g/m^3) exceeded for p % (0.1..99) at (lat, lon), alt (km)."""
    lon = lon % 360.0
    ps = np.array(P836_P)
    if p in P836_P:
        pb = pa = p
    else:
        k = int(np.searchsorted(ps, p)) - 1
        pb, pa = ps[k], ps[k + 1]
    r, c = (90.0 - lat) / 1.125, lon / 1.125
    R, C = int(np.floor(r)), int(np.floor(c))
    R = min(R, 159)
    corners = [(R, C), (R + 1, C), (R, C + 1), (R + 1, C + 1)]
    w = [(R + 1 - r) * (C + 1 - c), (r - R) * (C + 1 - c), (R + 1 - r) * (c - C), (r - R) * (c - C)]
    alt_i = []
    for (i, j) in corners:
        la, lo = 90.0 - i * 1.125, (j * 1.125) % 360.0
        alt_i.append(_bicubic(maps['topo'], (90.5 - la) / 0.5, (lo + 0.5) / 0.5))

    def at(pp):
        tot = 0.0
        for (i, j), wk, ak in zip(corners, w, alt_i):
            tot += wk * maps['rho'][pp][i, j] * np.exp(-(alt_km - ak) / maps['vsch'][pp][i, j])
        return tot
    if pa == pb:
        return at(pa)
    ra, rb = at(pa), at(pb)
    return rb + (ra - rb) * (np.log(p) - np.log(pb)) / (np.log(pa) - np.log(pb))


def p836_maps_from_itur():
    """Map arrays as distributed with ITU-Rpy (used for the reference vectors only)."""
    import os, itur
    d = os.path.join(os.path.dirname(itur.__file__), 'data', '836')
    ld = lambda f: (lambda z: z[z.files[0]])(np.load(os.path.join(d, f)))
    tag = lambda p: str(p).replace('.0', '').replace('.', '') if p < 1 else str(int(p))
    return {'rho': {p: ld('v6_rho_%s.npz' % tag(p)) for p in P836_P},
            'vsch': {p: ld('v6_vsch_%s.npz' % tag(p)) for p in P836_P},
            'topo': ld('v6_topo_0dot5.npz')}


def rho0_from_p836(maps, p1, lat, lon, he):
    """Sea-level density rho0 for Ag(p1): P.836 density exceeded for q = 100 - p1 %
    (limited to 0.1..99 %) at the station altitude he, scaled with exp(he/2)."""
    q = min(max(100.0 - p1, 0.1), 99.0)
    rho_s = p836_rho(maps, lat, lon, q, he)
    return q, rho_s, rho_s * np.exp(he / 2.0)


# ---------------------------------------------------------------------------
# Rec. ITU-R P.453-14 §2.2: wet term of the surface refractivity Nwet exceeded
# for p % of an average year (annual maps), bilinear interpolation (P.1144).
# Map arrays on the 0.75 deg grid (241 x 481, rows lat -90 -> 90, columns
# lon -180 -> 180); between the 18 map percentages Nwet is interpolated
# linearly in log(p) (the result does not depend on the base of the log).
# ---------------------------------------------------------------------------

P453_P = P836_P


def _bilinear_p453(grid, lat, lon):
    r, c = (lat + 90.0) / 0.75, (lon + 180.0) / 0.75
    R, C = min(int(np.floor(r)), 239), min(int(np.floor(c)), 479)
    fr, fc = r - R, c - C
    return ((1 - fr) * (1 - fc) * grid[R, C] + (1 - fr) * fc * grid[R, C + 1]
            + fr * (1 - fc) * grid[R + 1, C] + fr * fc * grid[R + 1, C + 1])


def nwet_p453(maps, lat, lon, p=50.0):
    """Nwet (ppm) exceeded for p % (0.1..99) at (lat, lon); lon in any range."""
    lon = (lon + 180.0) % 360.0 - 180.0         # columns -180 and 180 are identical
    ps = np.array(P453_P)
    if p in P453_P:
        return _bilinear_p453(maps[p], lat, lon)
    k = int(np.searchsorted(ps, p)) - 1
    pb, pa = ps[k], ps[k + 1]
    nb, na = _bilinear_p453(maps[pb], lat, lon), _bilinear_p453(maps[pa], lat, lon)
    return nb + (na - nb) * (np.log(p) - np.log(pb)) / (np.log(pa) - np.log(pb))


def p453_maps_from_itur():
    """Nwet annual maps as distributed with ITU-Rpy (identical to the P.453-14 files)."""
    import os, itur
    d = os.path.join(os.path.dirname(itur.__file__), 'data', '453')
    tag = lambda p: str(p).replace('.0', '').replace('.', '') if p < 1 else str(int(p))
    return {p: (lambda z: z[z.files[0]])(np.load(os.path.join(d, 'v13_nwet_annual_%s.npz' % tag(p))))
            for p in P453_P}
