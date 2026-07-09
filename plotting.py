"""
Plotting routines for the TW-MZM EO model.

Each fig_* function returns a matplotlib Figure; save_figure() writes it to
disk.  Call apply_style() once before building figures.

Styling follows a validated colorblind-safe categorical palette with fixed
slot order (blue, aqua, yellow, green, violet, ...) plus a sequential blue
ramp for ordered families of curves; chart chrome is recessive so the data
carry the figure.
"""
from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np

import model

# validated categorical palette -- fixed slot order, never cycled/reordered
PALETTE = ["#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7",
           "#e34948", "#e87ba4", "#eb6834"]
SEQ_BLUES = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]  # ordinal, light->dark
SURFACE = "#fcfcfb"
INK = "#0b0b0b"       # primary text
INK2 = "#52514e"      # secondary text (axis labels, direct labels)
MUTED = "#898781"     # ticks, reference lines
GRID = "#e1e0d9"      # hairline grid
BASELINE = "#c3c2b7"  # axis spines


def apply_style():
    """Recessive engineering-plot style; call once before making figures."""
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "DejaVu Sans", "Arial"],
        "font.size": 10.5,
        "axes.edgecolor": BASELINE,
        "axes.linewidth": 1.0,
        "axes.labelcolor": INK2,
        "axes.titlecolor": INK,
        "axes.titlesize": 11.5,
        "axes.titlepad": 10.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK2,
        "ytick.labelcolor": INK2,
        "legend.frameon": False,
        "lines.linewidth": 2.0,
    })


def save_figure(fig, name, outdir):
    """Save a figure as PNG into outdir and return the path."""
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{name}.png"
    fig.savefig(path, dpi=170)
    return path


def _spread(values, min_gap):
    """Nudge label y-positions apart (ascending sweep) to avoid collisions."""
    values = np.asarray(values, dtype=float)
    out = values.copy()
    order = np.argsort(values)
    prev = -np.inf
    for idx in order:
        out[idx] = max(out[idx], prev + min_gap)
        prev = out[idx]
    return out


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def fig_line_parameters(line_loaded, line_metal=None, optical=None, label=None):
    """Attenuation, effective index and Z0 of the electrode vs frequency."""
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.0), constrained_layout=True)
    if label:
        fig.suptitle(label, fontsize=12.5, color=INK)
    fg = line_loaded.f / 1e9

    ax = axes[0]
    ax.plot(fg, line_loaded.alpha_db_cm, color=PALETTE[0], label="loaded line")
    if line_metal is not None:
        ax.plot(line_metal.f / 1e9, line_metal.alpha_db_cm,
                color=PALETTE[1], label="metal-only")
        ax.legend(loc="upper left")
    ax.set_title("RF attenuation")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("α (dB/cm)")

    ax = axes[1]
    ax.plot(fg, line_loaded.n_rf, color=PALETTE[0])
    if line_metal is not None:
        ax.plot(line_metal.f / 1e9, line_metal.n_rf, color=PALETTE[1])
    if optical is not None:
        ax.axhline(optical.n_g, color=MUTED, lw=1.2, ls=(0, (4, 3)))
        ax.text(0.98, optical.n_g, " n_g (optical)", color=INK2, fontsize=9.5,
                va="bottom", ha="right",
                transform=ax.get_yaxis_transform())
    ax.set_title("Effective microwave index")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("n_RF")

    ax = axes[2]
    ax.plot(fg, line_loaded.z0.real, color=PALETTE[0])
    ax.plot(fg, line_loaded.z0.imag, color=PALETTE[0], ls=(0, (4, 3)), lw=1.6)
    if line_metal is not None:
        ax.plot(line_metal.f / 1e9, line_metal.z0.real, color=PALETTE[1])
        ax.plot(line_metal.f / 1e9, line_metal.z0.imag, color=PALETTE[1],
                ls=(0, (4, 3)), lw=1.6)
    for line, dy in ((line_loaded, 0.0),):
        ax.text(fg[-1], line.z0.real[-1] + dy, " Re", color=INK2,
                fontsize=9.5, va="center", ha="left", clip_on=False)
        ax.text(fg[-1], line.z0.imag[-1] + dy, " Im", color=INK2,
                fontsize=9.5, va="center", ha="left", clip_on=False)
    ax.set_title("Characteristic impedance")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Z0 (Ω)")
    ax.set_xlim(0, fg[-1] * 1.06)
    return fig


def fig_eo_s21(res, label=None):
    """Normalized EO S21 with bandwidth markers."""
    fig, ax = plt.subplots(figsize=(7.8, 4.7), constrained_layout=True)
    fg = res.f / 1e9
    ax.plot(fg, res.s21_db, color=PALETTE[0], lw=2.3)

    for drop, fbw, dy in ((3.0, res.f3db, 10), (6.0, res.f6db, -16)):
        ax.axhline(-drop, color=MUTED, lw=0.9, ls=(0, (4, 3)))
        if np.isfinite(fbw):
            ax.plot(fbw / 1e9, -drop, marker="o", ms=8, mfc=SURFACE,
                    mec=PALETTE[0], mew=2.0, zorder=5)
            ax.annotate(f"−{drop:.0f} dB el.:  {fbw / 1e9:.1f} GHz",
                        xy=(fbw / 1e9, -drop), xytext=(10, dy),
                        textcoords="offset points", color=INK2, fontsize=10)

    ax.set_title("Electro-optic frequency response (normalized)"
                 + (f" — {label}" if label else ""))
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("EO S21 (dB, electrical)")
    ax.set_xlim(0, fg[-1])
    ax.set_ylim(min(res.s21_db.min() - 1.5, -10.0), 1.0)
    return fig


def fig_breakdown(breakdown, label=None):
    """Full EO response with each bandwidth limiter isolated."""
    fig, ax = plt.subplots(figsize=(8.6, 5.0), constrained_layout=True)
    names = list(breakdown)
    fg = breakdown[names[0]].f / 1e9

    ends = []
    for i, name in enumerate(names):
        res = breakdown[name]
        if name == "Full model":
            ax.plot(fg, res.s21_db, color=PALETTE[0], lw=2.6)
        else:
            ax.plot(fg, res.s21_db, color=PALETTE[i], lw=1.7, ls=(0, (5, 3)))
        ends.append(res.s21_db[-1])

    ymin = min(e.s21_db.min() for e in breakdown.values()) - 2.0
    ymax = max(e.s21_db.max() for e in breakdown.values()) + 1.0
    span = ymax - ymin
    ylab = _spread(ends, 0.045 * span)
    for i, name in enumerate(names):
        ax.text(fg[-1] * 1.012, ylab[i], name, color=INK2, fontsize=9.5,
                va="center", ha="left")

    ax.axhline(-3.0, color=MUTED, lw=0.9, ls=(0, (4, 3)))
    ax.text(1.0, -3.0, "−3 dB", color=MUTED, fontsize=9, va="bottom")
    ax.set_title("Bandwidth-limiter breakdown (isolated effects)"
                 + (f" — {label}" if label else ""))
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("EO S21 (dB, electrical)")
    ax.set_xlim(0, fg[-1] * 1.30)
    ax.set_ylim(ymin, ymax)
    return fig


def fig_rf_sparams(line, drive, z_ref=50.0, label=None):
    """|S21| and |S11| of the bare electrode two-port in a z_ref reference."""
    s11, s21 = model.line_sparams(line, drive.length, z_ref)
    fg = line.f / 1e9
    fig, ax = plt.subplots(figsize=(7.8, 4.7), constrained_layout=True)
    ax.plot(fg, 20 * np.log10(np.maximum(np.abs(s21), 1e-4)),
            color=PALETTE[0], label="|S21|")
    ax.plot(fg, 20 * np.log10(np.maximum(np.abs(s11), 1e-4)),
            color=PALETTE[1], label="|S11|")
    ax.legend(loc="lower left")
    ax.set_title(f"Electrode two-port S-parameters ({z_ref:.0f} Ω reference)"
                 + (f" — {label}" if label else ""))
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_xlim(0, fg[-1])
    ax.set_ylim(-40, 2)
    return fig


def fig_voltage_profile(line, drive, freqs_ghz=(2.0, 10.0, 25.0, 50.0),
                        label=None):
    """|V(z)|/Vs along the electrode at several frequencies (standing waves)."""
    fig, ax = plt.subplots(figsize=(7.8, 4.7), constrained_layout=True)
    for i, f0 in enumerate(freqs_ghz[:len(SEQ_BLUES)]):
        z, v, fact = model.voltage_profile(line, drive, f0 * 1e9)
        label = f"{fact / 1e9:.0f} GHz"
        ax.plot(z * 1e3, np.abs(v), color=SEQ_BLUES[i], label=label)
        ax.text(z[-1] * 1e3 * 1.012, np.abs(v[-1]), label, color=INK2,
                fontsize=9.5, va="center", ha="left")
    ax.legend(loc="lower left", title="drive frequency", title_fontsize=9.5)
    ax.set_title("Electrode voltage standing-wave profile"
                 + (f" — {label}" if label else ""))
    ax.set_xlabel("Position z (mm)")
    ax.set_ylabel("|V(z)| / V_source")
    ax.set_xlim(0, drive.length * 1e3 * 1.14)
    return fig
