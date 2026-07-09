# TW-MZM EO response model

Python model that predicts the small-signal EO S21 and 3 dB bandwidth of a
silicon traveling-wave Mach-Zehnder modulator by combining parameters from
physical simulators. It must NOT attempt to replace those simulators: the
user (Mario, silicon-photonics engineer) extracts junction R, C with
Lumerical CHARGE, group index n_g with Lumerical MODE, and line parameters
f, n_RF(f), Z0(f), alpha(f) with Ansys HFSS.

## Commands

- `python run_example.py` — runs the 10 self-tests, both example devices in
  both HFSS modes, prints summaries, writes `figures\*_gsg.png` /
  `*_gsgsg.png`. This is the verification gate: run it after any change to
  the engine; all self-tests must pass.
- No package/venv/test infra beyond this; requires numpy, scipy, matplotlib
  (Anaconda Python 3.11 on this machine).

## Architecture

- `inputs.py` — dataclasses (`RFLine`, `Junction`, `Optical`, `Drive`),
  unit converters, placeholder generators `example_inputs()` (GSG) and
  `example_inputs_differential()` (GSGSG series push-pull).
- `model.py` — analytical traveling-wave engine (`eo_response`), loading
  synthesis (`load_line`), diagnostics (`response_breakdown`,
  `line_sparams`, `voltage_profile`, `combine_dual_drive`), `self_test()`.
- `plotting.py` — figure routines (validated colorblind-safe palette).
- `README.md` — physics, equations, assumptions, references (canonical
  m(f) reference: Gopalakrishnan et al., JLT 12(10), 1994).
- `docs/` — dated session notes for the human; read them for history and
  decisions before proposing changes.

## Critical conventions (do not violate silently)

- Junction units: R in **Ohm*m** (resistance x length; 10 Ohm*mm = 1e-2),
  C in **F/m**. R*C is the junction time constant. R is NEVER "Ohm per m".
- alpha is stored in Np/m internally; constructors take
  `alpha_unit in {dB/cm, dB/mm, dB/m, Np/m}`.
- EO S21 = 20*log10|m(f)/m(0)| in ELECTRICAL dB, normalized to DC;
  `f3db` = first crossing of -3 dB (DC-referenced, not peak-referenced);
  `f6db` (-6 dB electrical) = the 3 dB OPTICAL bandwidth.
- Phasors e^{+jwt}; forward wave e^{-gamma z}.
- Two HFSS modes via `RFLine.includes_junction`: True = loaded-line data
  used directly; False = model adds junction admittance (`load_line`).
  The RC divider H_RC applies in BOTH modes (it converts electrode voltage
  to depletion-cap voltage; not double counting).
- Topologies: single junction -> per-junction `Junction(R_j, C_j)`;
  series push-pull (GSG or differential GSGSG) ->
  `Junction.series_push_pull(R_j, C_j)` (effective 2R_j, C_j/2 —
  MANDATORY in mode B, immaterial in mode A); differential drive uses
  odd-mode n_RF/alpha, differential Z0 and differential Z_S/Z_L (e.g.
  100 Ohm); two independent electrodes -> `combine_dual_drive`.

## Scope boundaries (user's explicit choices)

- In scope: EO S21, bandwidths, RF diagnostics, mismatch/termination
  studies (e.g. under-termination peaking).
- Out of scope: Vpi/VpiL, plasma dispersion, carrier dynamics.
- SPICE/compact models: DEFERRED, not rejected (user reconsidered at the
  end of the first session). Build a netlist/Verilog-A export only when a
  concrete driver co-simulation, transient/eye, or EDA-handoff need
  appears; validate its small-signal EO response against this model (the
  golden reference) before trusting transients.
- Placeholder data in `inputs.example_inputs*()` is FICTIONAL (plausible
  values, invented curves) and is meant to be replaced by the user's real
  CHARGE/MODE/HFSS exports.

## When extending the engine

Add a closed-form limiting-case check to `model.self_test()` for any new
physics, mirroring the existing style (matched/lossless/velocity-matched
flatness, analytic loss and walk-off limits, mode A == mode B equality).
