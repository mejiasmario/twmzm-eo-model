"""
End-to-end demo of the TW-MZM EO response model on PLACEHOLDER data.

Runs the analytical engine's self-tests, then evaluates two example
devices -- a single-ended GSG modulator and a GSGSG series-push-pull
modulator driven differentially -- each in both HFSS interpretation modes
(loaded-line data vs metal-only data plus loading synthesis), prints the
result summaries, and writes the diagnostic figures to ./figures with
_gsg / _gsgsg suffixes.

Usage:
    python run_example.py           # save figures only
    python run_example.py --show    # also open interactive windows

To use your own data, replace inputs.example_inputs() with your CHARGE /
MODE / HFSS exports (see README.md).
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--show", action="store_true",
                        help="open interactive figure windows")
    parser.add_argument("--outdir", default=None,
                        help="output directory for figures (default ./figures)")
    args = parser.parse_args(argv)

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")

    import model
    import plotting
    from inputs import example_inputs, example_inputs_differential

    # ---- sanity: closed-form limiting cases -------------------------------
    print("== self-tests ==")
    failed = False
    for name, ok, detail in model.self_test():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}  ({detail})")
        failed |= not ok
    if failed:
        sys.exit("self-tests failed -- aborting")

    # ---- evaluate both example devices in both HFSS modes ------------------
    plotting.apply_style()
    outdir = (pathlib.Path(args.outdir) if args.outdir
              else pathlib.Path(__file__).parent / "figures")
    examples = [
        ("gsg", "GSG single-ended", example_inputs),
        ("gsgsg", "GSGSG push-pull, differential", example_inputs_differential),
    ]
    figures = {}
    for tag, label, factory in examples:
        print(f"\n== example: {label} ==")
        # mode B: HFSS gave the metal-only electrode; model adds the loading
        rf_metal, junction, optical, drive = factory(
            hfss_includes_junction=False)
        res_b = model.eo_response(rf_metal, junction, optical, drive)
        # mode A: HFSS data already include the doped-silicon loading
        rf_loaded, *_ = factory(hfss_includes_junction=True)
        res_a = model.eo_response(rf_loaded, junction, optical, drive)

        dev = float(np.max(np.abs(res_a.s21_db - res_b.s21_db)))
        print(f"mode A (loaded HFSS) vs mode B (metal-only + synthesis): "
              f"max deviation {dev:.2e} dB\n")
        print(model.summary(res_a, junction, optical, drive))

        figures[f"line_parameters_{tag}"] = plotting.fig_line_parameters(
            res_a.line, rf_metal, optical, label=label)
        figures[f"eo_s21_{tag}"] = plotting.fig_eo_s21(res_a, label=label)
        figures[f"bandwidth_breakdown_{tag}"] = plotting.fig_breakdown(
            model.response_breakdown(rf_loaded, junction, optical, drive),
            label=label)
        figures[f"rf_sparams_{tag}"] = plotting.fig_rf_sparams(
            res_a.line, drive, z_ref=float(abs(drive.z_src)), label=label)
        figures[f"voltage_profile_{tag}"] = plotting.fig_voltage_profile(
            res_a.line, drive, label=label)

    print()
    for name, fig in figures.items():
        print("wrote", plotting.save_figure(fig, name, outdir))

    if args.show:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == "__main__":
    main()
