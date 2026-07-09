"""
Analytical traveling-wave MZM electro-optic frequency response.

Physics
-------
The electrode is a uniform transmission line with frequency-dependent
propagation constant gamma(f) = alpha + j*beta and characteristic impedance
Z0(f), driven at z = 0 by a source of impedance Zs and terminated at z = L
with Zl.  The line voltage is the superposition of a forward and a backward
(reflected) wave; an optical wavefront co-propagating at the group velocity
c/n_g accumulates phase proportional to the local *junction* voltage.

    gamma(f) = alpha_Np(f) + j*2*pi*f*n_rf(f)/c          (HFSS)
    beta_o(f) = 2*pi*f*n_g/c                             (MODE)
    H_RC(f) = 1/(1 + j*2*pi*f*R*C)                       (CHARGE)

    Gl = (Zl - Z0)/(Zl + Z0),   Gs = (Zs - Z0)/(Zs + Z0)
    V+/Vs = Z0 / [(Z0 + Zs) * (1 - Gs*Gl*exp(-2*gamma*L))]

    F(x) = (1 - exp(-x))/x
    x_co = (gamma - j*beta_o)*L        co-propagating walk-off/loss term
    x_ct = (gamma + j*beta_o)*L        counter-propagating (reflected wave)

    m(f) = H_RC * (V+/Vs) * [ F(x_co) + Gl*exp(-x_co)*F(x_ct) ]

    EO S21(f) = 20*log10 |m(f)/m(0)| ,   m(0) = Zl/(Zs + Zl)

Conventions
-----------
* EO S21 is in ELECTRICAL dB: 20*log10 of the normalized modulation
  amplitude.  The "-3 dB EO bandwidth" is where this crosses -3 dB
  (equivalently 1.5 dB optical); the -6 dB electrical point equals the
  3 dB OPTICAL bandwidth and is also reported.
* If RFLine.includes_junction is False, the junction admittance
  Y_j = j*w*C/(1 + j*w*R*C) is first added per unit length to the HFSS
  (metal-only) line and the loaded gamma/Z0 are computed (load_line).
* The junction RC divider H_RC applies in BOTH cases: it converts electrode
  voltage into depletion-capacitance voltage and is not part of gamma/Z0.
* Electrode topology: the engine models ONE quasi-TEM mode on ONE line, so
  GSG vs GSGSG enters only through the inputs (fixed x2 / x1/2 topology
  factors move Vpi, which is out of scope, and cancel in the normalized
  response).  Single junction: per-junction R, C.  Series push-pull (GSG
  or GSGSG) and differential GSGSG: build the effective load with
  Junction.series_push_pull and use differential-mode line and drive
  quantities (odd-mode n_rf/alpha, differential Z0, e.g. 100 Ohm system).
  Two independent electrodes: model each arm, then combine_dual_drive().
  Differential modeling assumes a symmetric layout and purely differential
  excitation (no common mode / mode conversion).

Main entry point: eo_response(rf, junction, optical, drive) -> EOResponse.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.constants import c as C0

from inputs import RFLine, Junction, Optical, Drive, alpha_to_np_per_m


# ---------------------------------------------------------------------------
# numerics
# ---------------------------------------------------------------------------

def _one_minus_exp_over(x):
    """(1 - exp(-x))/x for complex x, series-protected near x = 0."""
    x = np.asarray(x, dtype=complex)
    small = np.abs(x) < 1e-6
    safe = np.where(small, 1.0, x)
    full = (1.0 - np.exp(-safe)) / safe
    series = 1.0 - x / 2.0 + x * x / 6.0 - x * x * x / 24.0
    return np.where(small, series, full)


# ---------------------------------------------------------------------------
# junction loading of a metal-only line (mode B)
# ---------------------------------------------------------------------------

def load_line(rf: RFLine, junction: Junction) -> RFLine:
    """Add the distributed junction admittance to a metal-only HFSS line.

    Telegrapher synthesis:  Z' = gamma*Z0,  Y' = gamma/Z0 + Y_junction,
    then gamma_loaded = sqrt(Z'*Y') and Z0_loaded = sqrt(Z'/Y').
    """
    if rf.includes_junction:
        raise ValueError("this RFLine already includes the junction loading")
    z = rf.z_series
    y = rf.y_shunt + junction.y_per_m(rf.f)
    gamma = np.sqrt(z * y)
    gamma = np.where(gamma.real < 0.0, -gamma, gamma)   # forward-decaying wave
    z0 = np.sqrt(z / y)
    z0 = np.where(z0.real < 0.0, -z0, z0)               # passive impedance
    n_rf = gamma.imag * C0 / (2.0 * np.pi * rf.f)
    return RFLine(rf.f, n_rf, z0, gamma.real, includes_junction=True)


# ---------------------------------------------------------------------------
# electro-optic response
# ---------------------------------------------------------------------------

@dataclass
class EOResponse:
    """Result container returned by eo_response()."""
    f: np.ndarray        # frequency grid [Hz]
    m: np.ndarray        # complex EO response, normalized to its DC value
    s21_db: np.ndarray   # 20*log10 |m|  (electrical dB)
    f3db: float          # -3 dB electrical crossing [Hz] (nan if not reached)
    f6db: float          # -6 dB electrical = -3 dB optical crossing [Hz]
    line: RFLine | None  # loaded line used (None for combined dual-drive)


def eo_response(rf: RFLine, junction: Junction, optical: Optical,
                drive: Drive) -> EOResponse:
    """Small-signal EO frequency response of the TW-MZM.

    If rf.includes_junction is False the junction loading is added first
    (see load_line); the RC voltage divider is applied in both cases.
    """
    line = rf if rf.includes_junction else load_line(rf, junction)
    f = line.f
    L = drive.length
    g = line.gamma
    z0 = line.z0
    beta_o = 2.0 * np.pi * f * optical.n_g / C0

    gl = (drive.z_load - z0) / (drive.z_load + z0)
    gs = (drive.z_src - z0) / (drive.z_src + z0)
    vfwd = z0 / ((z0 + drive.z_src) * (1.0 - gs * gl * np.exp(-2.0 * g * L)))

    x_co = (g - 1j * beta_o) * L
    x_ct = (g + 1j * beta_o) * L
    m = junction.h_rc(f) * vfwd * (_one_minus_exp_over(x_co)
                                   + gl * np.exp(-x_co) * _one_minus_exp_over(x_ct))

    m0 = drive.z_load / (drive.z_src + drive.z_load)    # exact DC limit
    m = m / m0
    s21_db = 20.0 * np.log10(np.abs(m))
    return EOResponse(f, m, s21_db,
                      find_crossing(f, s21_db, 3.0),
                      find_crossing(f, s21_db, 6.0),
                      line)


def find_crossing(f, s21_db, drop_db=3.0):
    """First frequency where s21_db falls below -|drop_db| (linear interp).

    Returns nan when the response never drops that far in the simulated band.
    """
    target = -abs(drop_db)
    below = np.asarray(s21_db) <= target
    if not below.any():
        return float("nan")
    i = int(np.argmax(below))
    if i == 0:
        return float(f[0])
    f1, f2 = f[i - 1], f[i]
    s1, s2 = s21_db[i - 1], s21_db[i]
    return float(f1 + (target - s1) * (f2 - f1) / (s2 - s1))


def combine_dual_drive(res_arm1: EOResponse, res_arm2: EOResponse) -> EOResponse:
    """MZM response when the two arms have independent electrodes (dual drive).

    Assumes equal-amplitude complementary (push-pull) drive and identical
    DC voltage division in both drive networks, so the normalized MZM
    response is the average of the two normalized arm responses,
    m = (m1 + m2)/2 (identical arms reduce to either arm).  If the arms'
    DC networks differ, weight each m by its Zl/(Zs + Zl) before averaging.

    The combined result has line=None; run the per-line diagnostics
    (line_sparams, voltage_profile, response_breakdown) on each arm.
    """
    if (res_arm1.f.shape != res_arm2.f.shape
            or not np.allclose(res_arm1.f, res_arm2.f)):
        raise ValueError("arm responses must be on the same frequency grid")
    f = res_arm1.f
    m = 0.5 * (res_arm1.m + res_arm2.m)
    s21_db = 20.0 * np.log10(np.abs(m))
    return EOResponse(f, m, s21_db,
                      find_crossing(f, s21_db, 3.0),
                      find_crossing(f, s21_db, 6.0),
                      None)


# ---------------------------------------------------------------------------
# diagnostics
# ---------------------------------------------------------------------------

def line_sparams(line: RFLine, length, z_ref=50.0):
    """(S11, S21) of the bare electrode two-port in a z_ref reference.

    Useful to compare with a full-structure HFSS/VNA measurement.  Rule of
    thumb: for a matched, loss-limited modulator the -3 dB EO point falls
    where the electrical line loss reaches about 6.4 dB.
    """
    gL = line.gamma * length
    ch, sh = np.cosh(gL), np.sinh(gL)
    a, b, c, d = ch, line.z0 * sh, sh / line.z0, ch
    den = a + b / z_ref + c * z_ref + d
    s11 = (a + b / z_ref - c * z_ref - d) / den
    s21 = 2.0 / den
    return s11, s21


def voltage_profile(line: RFLine, drive: Drive, f0, n_z=401):
    """Electrode voltage V(z)/Vs at the grid frequency nearest f0 [Hz].

    Returns (z [m], V(z)/Vs complex, actual frequency used [Hz]).
    """
    i = int(np.argmin(np.abs(line.f - f0)))
    g, z0 = line.gamma[i], line.z0[i]
    L = drive.length
    gl = (drive.z_load - z0) / (drive.z_load + z0)
    gs = (drive.z_src - z0) / (drive.z_src + z0)
    vfwd = z0 / ((z0 + drive.z_src) * (1.0 - gs * gl * np.exp(-2.0 * g * L)))
    z = np.linspace(0.0, L, n_z)
    v = vfwd * (np.exp(-g * z) + gl * np.exp(-g * (2.0 * L - z)))
    return z, v, float(line.f[i])


def response_breakdown(rf: RFLine, junction: Junction, optical: Optical,
                       drive: Drive):
    """EO S21 with each bandwidth limiter isolated (all curves in dB).

    Returns an ordered dict-like mapping:
      'Full model'      - everything on
      'Junction RC only'- lossless, velocity-matched, matched terminations
      'RF loss only'    - alpha kept; walk-off, mismatch, RC off
      'Walk-off only'   - n_rf kept; loss, mismatch, RC off
      'Mismatch only'   - real terminations kept; loss, walk-off, RC off
    The dB curves add up approximately to the full model (cross terms are
    small), which makes the dominant limiter obvious at a glance.
    """
    line = rf if rf.includes_junction else load_line(rf, junction)
    zeros = np.zeros_like(line.alpha_np_m)
    ng = np.full_like(line.n_rf, optical.n_g)

    lossless = replace(line, alpha_np_m=zeros)
    matched_v = replace(line, n_rf=ng)
    ideal = replace(line, alpha_np_m=zeros, n_rf=ng)
    matched_z = replace(drive, z_src=line.z0, z_load=line.z0)
    no_rc = Junction(0.0, junction.c_f_m)

    return {
        "Full model": eo_response(line, junction, optical, drive),
        "Junction RC only": eo_response(ideal, junction, optical, matched_z),
        "RF loss only": eo_response(matched_v, no_rc, optical, matched_z),
        "Walk-off only": eo_response(lossless, no_rc, optical, matched_z),
        "Mismatch only": eo_response(ideal, no_rc, optical, drive),
    }


def summary(res: EOResponse, junction: Junction, optical: Optical,
            drive: Drive) -> str:
    """Human-readable result summary (plain ASCII)."""
    line = res.line

    def fmt_f(x):
        return "beyond simulated band" if np.isnan(x) else f"{x / 1e9:6.2f} GHz"

    rows = [
        "TW-MZM EO response summary",
        "-" * 54,
        f"electrode length             : {drive.length * 1e3:.2f} mm",
        f"junction R*C                 : {junction.tau * 1e12:.2f} ps"
        f"   (f_RC = {junction.f_rc / 1e9:.1f} GHz)",
        f"EO bandwidth, -3 dB electrical: {fmt_f(res.f3db)}",
        f"EO bandwidth, -6 dB electrical: {fmt_f(res.f6db)}  (= -3 dB optical)",
    ]
    if not np.isnan(res.f3db) and res.line is not None:
        f3 = res.f3db
        a3 = np.interp(f3, line.f, line.alpha_db_cm)
        n3 = np.interp(f3, line.f, line.n_rf)
        z3 = np.interp(f3, line.f, line.z0.real) \
            + 1j * np.interp(f3, line.f, line.z0.imag)
        gl = abs((drive.z_load - z3) / (drive.z_load + z3))
        gs = abs((drive.z_src - z3) / (drive.z_src + z3))
        rows += [
            "at f = f3dB:",
            f"  RF line loss               : {a3:.2f} dB/cm"
            f"   (alpha*L = {a3 * drive.length * 100:.2f} dB)",
            f"  index walk-off             : n_rf = {n3:.2f}"
            f" vs n_g = {optical.n_g:.2f}",
            f"  impedance                  : |Z0| = {abs(z3):.1f} Ohm"
            f"   |Gamma_load| = {gl:.2f}, |Gamma_src| = {gs:.2f}",
        ]
    return "\n".join(rows)


# ---------------------------------------------------------------------------
# self-tests against known limiting cases
# ---------------------------------------------------------------------------

def self_test():
    """Check the engine against closed-form limiting cases.

    Returns a list of (name, passed, detail) tuples.
    """
    out = []
    f = np.linspace(0.1e9, 60e9, 240)
    ones = np.ones_like(f)
    ng = 3.8
    opt = Optical(ng)
    jn0 = Junction(0.0, 0.0)
    drv = Drive(2.0e-3, 45.0, 45.0)

    # 1) matched, lossless, velocity-matched -> exactly flat 0 dB
    line = RFLine(f, ng * ones, 45.0 * ones, 0.0 * ones, includes_junction=True)
    r = eo_response(line, jn0, opt, drv)
    dev = float(np.max(np.abs(r.s21_db)))
    out.append(("flat response when matched/lossless/velocity-matched",
                dev < 1e-9, f"max deviation {dev:.1e} dB"))

    # 2) loss only -> |m| = (1 - exp(-aL))/(aL)
    a = 400.0  # Np/m
    line = RFLine(f, ng * ones, 45.0 * ones, a * ones, includes_junction=True)
    r = eo_response(line, jn0, opt, drv)
    expect = (1.0 - np.exp(-a * drv.length)) / (a * drv.length)
    dev = float(np.max(np.abs(np.abs(r.m) - expect)))
    out.append(("loss-limited response equals (1-e^-aL)/aL",
                dev < 1e-12, f"max deviation {dev:.1e}"))

    # 3) walk-off only -> |m| = |2 sin(theta/2)/theta|
    nrf = 5.2
    line = RFLine(f, nrf * ones, 45.0 * ones, 0.0 * ones, includes_junction=True)
    r = eo_response(line, jn0, opt, drv)
    theta = 2.0 * np.pi * f * (nrf - ng) * drv.length / C0
    expect = np.abs(2.0 * np.sin(theta / 2.0) / theta)
    dev = float(np.max(np.abs(np.abs(r.m) - expect)))
    out.append(("walk-off-limited response equals |sinc|",
                dev < 1e-12, f"max deviation {dev:.1e}"))

    # 4) DC limit -> 0 dB even with mismatched terminations and junction RC
    f4 = np.geomspace(1e5, 60e9, 240)
    a4 = alpha_to_np_per_m(0.5 * np.sqrt(f4 / 1e9), "dB/cm")
    line = RFLine(f4, 4.1 * np.ones_like(f4), (35.0 - 2.0j) * np.ones_like(f4),
                  a4, includes_junction=True)
    jn = Junction(8.0e-3, 2.4e-10)
    r = eo_response(line, jn, opt, Drive(2.5e-3, 50.0, 50.0))
    dev = float(abs(r.s21_db[0]))
    out.append(("DC limit normalizes to 0 dB with mismatch + RC",
                dev < 0.01, f"|S21(f->0)| = {dev:.1e} dB"))

    # 5) junction loading trends: n_rf up, |Z0| down, loss up
    from inputs import example_inputs
    rf_u, jn_e, opt_e, drv_e = example_inputs(hfss_includes_junction=False)
    rl = load_line(rf_u, jn_e)
    ok = (np.all(rl.n_rf > rf_u.n_rf)
          and np.all(np.abs(rl.z0) < np.abs(rf_u.z0))
          and np.all(rl.alpha_np_m >= rf_u.alpha_np_m))
    out.append(("loading raises n_rf, lowers |Z0|, adds loss", bool(ok),
                f"n_rf {rf_u.n_rf[-1]:.2f}->{rl.n_rf[-1]:.2f}, "
                f"|Z0| {abs(rf_u.z0[-1]):.0f}->{abs(rl.z0[-1]):.0f} Ohm"))

    # 6) mode A (pre-loaded data) == mode B (loaded inside the model)
    ra = eo_response(rl, jn_e, opt_e, drv_e)
    rb = eo_response(rf_u, jn_e, opt_e, drv_e)
    dev = float(np.max(np.abs(ra.m - rb.m)))
    out.append(("HFSS-loaded and metal-only+synthesis modes agree",
                dev < 1e-12, f"max deviation {dev:.1e}"))

    # 7) passivity of the electrode two-port
    s11, s21 = line_sparams(rl, drv_e.length)
    mx = float(max(np.max(np.abs(s11)), np.max(np.abs(s21))))
    out.append(("electrode two-port is passive (|S| <= 1)",
                mx <= 1.0 + 1e-9, f"max |S| = {mx:.6f}"))

    # 8) series push-pull: same RC pole, half the loading admittance
    jn_1 = Junction(8.0e-3, 2.4e-10)
    jn_pp = Junction.series_push_pull(8.0e-3, 2.4e-10)
    fchk = np.linspace(1e9, 70e9, 60)
    ok = (np.isclose(jn_pp.tau, jn_1.tau, rtol=1e-15)
          and np.allclose(jn_pp.h_rc(fchk), jn_1.h_rc(fchk), rtol=1e-14)
          and np.allclose(jn_pp.y_per_m(fchk), 0.5 * jn_1.y_per_m(fchk),
                          rtol=1e-14))
    out.append(("series push-pull keeps the RC pole, halves the loading",
                bool(ok), f"tau = {jn_pp.tau * 1e12:.2f} ps in both"))

    # 9) differential (GSGSG) example: the two HFSS modes agree
    from inputs import example_inputs_differential
    rf_du, jn_d, opt_d, drv_d = example_inputs_differential(
        hfss_includes_junction=False)
    ra = eo_response(load_line(rf_du, jn_d), jn_d, opt_d, drv_d)
    rb = eo_response(rf_du, jn_d, opt_d, drv_d)
    dev = float(np.max(np.abs(ra.m - rb.m)))
    out.append(("differential example: loaded and metal-only modes agree",
                dev < 1e-12, f"max deviation {dev:.1e}"))

    # 10) dual-drive combination of identical arms equals a single arm
    rc = combine_dual_drive(ra, rb)
    dev = float(np.max(np.abs(rc.m - ra.m)))
    ok = dev < 1e-12 and bool(np.isclose(rc.f3db, ra.f3db, equal_nan=True))
    out.append(("dual-drive combine of identical arms is the identity",
                ok, f"max deviation {dev:.1e}"))
    return out
