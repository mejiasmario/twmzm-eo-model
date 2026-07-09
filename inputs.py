"""
Simulation inputs for the traveling-wave Mach-Zehnder modulator (TW-MZM)
electro-optic response model.

Data sources
------------
====================  ==========  ============================================
Tool                  Class       Quantities
====================  ==========  ============================================
Lumerical CHARGE      Junction    R [Ohm*m], C [F/m] at the operating bias
Lumerical MODE        Optical     optical group index n_g
Ansys HFSS            RFLine      f [Hz], n_rf(f), Z0(f) [Ohm], alpha(f)
(design choice)       Drive       length L [m], source/termination impedances
====================  ==========  ============================================

Unit conventions (used throughout the package)
----------------------------------------------
frequency : Hz
length    : m
R         : Ohm*m   resistance x length      (10 Ohm*mm  -> 1.0e-2 Ohm*m)
C         : F/m     capacitance per length   (250 fF/mm  -> 2.5e-10 F/m)
alpha     : stored as Np/m; constructors accept 'dB/cm', 'dB/mm', 'dB/m', 'Np/m'
Z0        : Ohm, complex allowed
phasors   : exp(+j*w*t); a forward wave varies as exp(-gamma*z)

Note on R: the junction access resistance halves when the device length
doubles, so its per-unit-length form is resistance*length (Ohm*m).  Only
then is the product R*C the length-independent junction time constant.

Replace the placeholder arrays in example_inputs() with your simulation
exports, or construct RFLine / Junction / Optical / Drive directly.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.constants import c as C0
from scipy.interpolate import PchipInterpolator

# ---------------------------------------------------------------------------
# attenuation unit handling
# ---------------------------------------------------------------------------

_NP_PER_DB = np.log(10.0) / 20.0  # 1 dB = 0.115129... Np

ALPHA_UNIT_TO_NP_PER_M = {
    "Np/m": 1.0,
    "dB/m": _NP_PER_DB,
    "dB/cm": 100.0 * _NP_PER_DB,
    "dB/mm": 1000.0 * _NP_PER_DB,
}


def alpha_to_np_per_m(alpha, unit="dB/cm"):
    """Convert an attenuation array to Np/m.

    `unit` is one of 'Np/m', 'dB/m', 'dB/cm', 'dB/mm'.
    """
    try:
        scale = ALPHA_UNIT_TO_NP_PER_M[unit]
    except KeyError:
        raise ValueError(
            f"unknown alpha unit {unit!r}; use one of {list(ALPHA_UNIT_TO_NP_PER_M)}"
        ) from None
    return np.asarray(alpha, dtype=float) * scale


# ---------------------------------------------------------------------------
# input containers
# ---------------------------------------------------------------------------

@dataclass
class RFLine:
    """Frequency-dependent transmission-line description of the electrode (HFSS).

    Attributes
    ----------
    f : 1-D array [Hz], strictly positive and ascending
    n_rf : microwave effective index (beta = 2*pi*f*n_rf/c)
    z0 : characteristic impedance [Ohm], complex allowed
    alpha_np_m : attenuation [Np/m]
    includes_junction :
        True  -> the data describe the junction-loaded line; the model only
                 applies the junction RC voltage divider on top.
        False -> the data describe the metal-only electrode; the model first
                 adds the junction admittance (model.load_line) and computes
                 the loaded gamma and Z0.
    """
    f: np.ndarray
    n_rf: np.ndarray
    z0: np.ndarray
    alpha_np_m: np.ndarray
    includes_junction: bool = True

    def __post_init__(self):
        self.f = np.asarray(self.f, dtype=float)
        self.n_rf = np.asarray(self.n_rf, dtype=float)
        self.z0 = np.asarray(self.z0, dtype=complex)
        self.alpha_np_m = np.asarray(self.alpha_np_m, dtype=float)
        if self.f.ndim != 1 or self.f.size < 2:
            raise ValueError("f must be a 1-D array with at least two points")
        if self.f[0] <= 0.0 or np.any(np.diff(self.f) <= 0.0):
            raise ValueError("f must be strictly positive and ascending")
        for name in ("n_rf", "z0", "alpha_np_m"):
            if getattr(self, name).shape != self.f.shape:
                raise ValueError(f"{name} must have the same shape as f")

    @classmethod
    def from_hfss(cls, f, n_rf, z0, alpha, alpha_unit="dB/cm",
                  includes_junction=True):
        """Build from HFSS exports, converting alpha from `alpha_unit`."""
        return cls(f, n_rf, z0, alpha_to_np_per_m(alpha, alpha_unit),
                   includes_junction)

    # ---- derived quantities ----------------------------------------------
    @property
    def gamma(self):
        """Complex propagation constant alpha + j*beta [1/m]."""
        return self.alpha_np_m + 1j * (2.0 * np.pi * self.f) * self.n_rf / C0

    @property
    def z_series(self):
        """Per-unit-length series impedance Z' = gamma*Z0 [Ohm/m]."""
        return self.gamma * self.z0

    @property
    def y_shunt(self):
        """Per-unit-length shunt admittance Y' = gamma/Z0 [S/m]."""
        return self.gamma / self.z0

    @property
    def alpha_db_cm(self):
        """Attenuation in dB/cm (for reporting and plots)."""
        return self.alpha_np_m / ALPHA_UNIT_TO_NP_PER_M["dB/cm"]

    def resampled(self, f_new):
        """Monotone (PCHIP) resampling onto a new frequency grid [Hz]."""
        f_new = np.asarray(f_new, dtype=float)

        def _interp(y):
            return PchipInterpolator(self.f, y, extrapolate=False)(f_new)

        z0 = _interp(self.z0.real) + 1j * _interp(self.z0.imag)
        return RFLine(f_new, _interp(self.n_rf), z0, _interp(self.alpha_np_m),
                      self.includes_junction)


@dataclass
class Junction:
    """Small-signal pn-junction loading (CHARGE) at the operating bias.

    r_ohm_m : access resistance x length [Ohm*m]   (10 Ohm*mm  -> 1.0e-2)
    c_f_m   : depletion capacitance / length [F/m] (250 fF/mm  -> 2.5e-10)
    """
    r_ohm_m: float
    c_f_m: float

    @classmethod
    def from_totals(cls, r_ohm, c_farad, length_m):
        """From total R [Ohm] and C [F] of a device/slice of known length [m]."""
        return cls(r_ohm * length_m, c_farad / length_m)

    @classmethod
    def series_push_pull(cls, r_j_ohm_m, c_j_f_m):
        """Effective line loading for two identical junctions electrically in
        series across the line -- covers single-drive series push-pull (GSG)
        and GSGSG driven differentially.

        Pass PER-JUNCTION values from CHARGE.  The driven mode sees the
        series combination R_eff = 2*R_j, C_eff = C_j/2: the time constant
        (hence the RC divider pole) is unchanged, but the loading admittance
        is halved -- which matters when the model adds the loading itself
        (RFLine.includes_junction=False).  With loaded HFSS data (True) the
        distinction is immaterial because only R*C enters.
        """
        return cls(2.0 * r_j_ohm_m, 0.5 * c_j_f_m)

    @property
    def tau(self):
        """Junction time constant R*C [s]."""
        return self.r_ohm_m * self.c_f_m

    @property
    def f_rc(self):
        """Intrinsic junction cutoff 1/(2*pi*R*C) [Hz] (inf if tau == 0)."""
        return np.inf if self.tau == 0.0 else 1.0 / (2.0 * np.pi * self.tau)

    def h_rc(self, f):
        """Junction voltage divider V_junction/V_electrode = 1/(1 + j*w*R*C)."""
        return 1.0 / (1.0 + 2j * np.pi * np.asarray(f, dtype=float) * self.tau)

    def y_per_m(self, f):
        """Distributed loading admittance of the series R-C branch [S/m]."""
        w = 2.0 * np.pi * np.asarray(f, dtype=float)
        return 1j * w * self.c_f_m / (1.0 + 1j * w * self.tau)


@dataclass
class Optical:
    """Optical waveguide data (MODE)."""
    n_g: float


@dataclass
class Drive:
    """Electrode drive configuration (design choice, not a simulation output).

    length : electrode (phase-shifter) length [m]
    z_src  : source impedance [Ohm], complex allowed
    z_load : termination impedance [Ohm], complex allowed (e.g. 1e9 ~ open)

    z_src/z_load may also be arrays on the RF frequency grid, which is used
    internally for the matched-termination diagnostic cases.
    """
    length: float
    z_src: complex = 50.0 + 0.0j
    z_load: complex = 50.0 + 0.0j


# ===========================================================================
# PLACEHOLDER example data -- REPLACE with your CHARGE / MODE / HFSS results
# ===========================================================================

def example_inputs(hfss_includes_junction=True):
    """Return (rf, junction, optical, drive) for a fictitious but physically
    plausible silicon depletion-mode TW-MZM -- single-ended GSG topology
    with one junction across the line per unit cell.

    The metal-only electrode curves below imitate a CPW/CPS HFSS export.
    With hfss_includes_junction=True, the returned line is *synthesized*
    from them via model.load_line, imitating an HFSS run whose stack-up
    contains the doped silicon.  Both settings therefore describe the same
    fictitious device and must yield the same EO response.
    """
    # ---- HFSS placeholder: metal-only electrode ---------------------------
    f = np.linspace(0.05e9, 70e9, 700)                    # Hz
    fghz = f / 1e9
    n_rf_u = 2.45 + 0.10 * np.exp(-fghz / 25.0)           # unloaded index
    z0_u = 61.0 - 2.0j * np.exp(-fghz / 8.0)              # Ohm
    alpha_u_db_cm = 0.35 * np.sqrt(fghz) + 0.020 * fghz   # conductor + diel.
    rf = RFLine.from_hfss(f, n_rf_u, z0_u, alpha_u_db_cm, "dB/cm",
                          includes_junction=False)

    # ---- CHARGE placeholder ------------------------------------------------
    junction = Junction(r_ohm_m=8.0e-3,   # 8 Ohm*mm
                        c_f_m=2.4e-10)    # 240 fF/mm  ->  f_RC ~ 83 GHz

    # ---- MODE placeholder ---------------------------------------------------
    optical = Optical(n_g=3.85)

    # ---- drive configuration (design choice) --------------------------------
    drive = Drive(length=2.5e-3, z_src=50.0, z_load=50.0)

    if hfss_includes_junction:
        from model import load_line  # lazy import avoids a circular import
        rf = load_line(rf, junction)
    return rf, junction, optical, drive


def example_inputs_differential(hfss_includes_junction=True):
    """Return (rf, junction, optical, drive) for a fictitious GSGSG series
    push-pull device driven DIFFERENTIALLY.

    All line quantities are differential-mode referenced: n_rf and alpha of
    the odd mode, Z0 is the differential impedance (= 2 x odd-mode Z0 of
    each half; mixed-mode HFSS extraction reports it directly), and the
    drive impedances are differential (100 Ohm system).  The junction is
    the EFFECTIVE differential load built from per-junction CHARGE values
    via Junction.series_push_pull (R_eff = 2*R_j, C_eff = C_j/2).

    Valid under the model's differential assumptions: symmetric layout and
    purely differential excitation (no common mode / mode conversion).
    """
    # ---- HFSS placeholder: metal-only GSGSG, differential (odd) mode ------
    f = np.linspace(0.05e9, 70e9, 700)                    # Hz
    fghz = f / 1e9
    n_rf_u = 2.40 + 0.10 * np.exp(-fghz / 25.0)           # odd-mode index
    z0_u = 110.0 - 3.0j * np.exp(-fghz / 8.0)             # differential Ohm
    alpha_u_db_cm = 0.30 * np.sqrt(fghz) + 0.018 * fghz
    rf = RFLine.from_hfss(f, n_rf_u, z0_u, alpha_u_db_cm, "dB/cm",
                          includes_junction=False)

    # ---- CHARGE placeholder: PER-JUNCTION values -> effective SPP load ----
    junction = Junction.series_push_pull(r_j_ohm_m=8.0e-3,   # 8 Ohm*mm
                                         c_j_f_m=2.4e-10)    # 240 fF/mm

    # ---- MODE placeholder ----------------------------------------------------
    optical = Optical(n_g=3.85)

    # ---- drive configuration: differential system -----------------------------
    drive = Drive(length=2.5e-3, z_src=100.0, z_load=100.0)

    if hfss_includes_junction:
        from model import load_line  # lazy import avoids a circular import
        rf = load_line(rf, junction)
    return rf, junction, optical, drive
