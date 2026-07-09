# Session notes — 2026-07-08 (initial build)

Human-readable record of the session in which this model was designed and
built. Condensed Q&A summaries first; the full detailed explanation of
**why a low load impedance increases bandwidth** is preserved at the end,
as requested.

---

## 1. What was decided and built

**Goal.** Predict EO S21 and 3 dB bandwidth of a silicon TW-MZM by
combining CHARGE (junction R, C), MODE (n_g), and HFSS (f, n_RF, Z0, α)
outputs — explicitly *not* replacing those simulators. Out of scope by
choice: Vπ/VπL, plasma dispersion, carrier dynamics, SPICE/compact models.

**Approach comparison verdict.** Analytical TW theory, numerically solved
telegrapher equations, and ABCD cascades are the *same physics* at
different generality; the lumped equivalent circuit is a coarse
discretization of it. For a z-uniform electrode described by γ(f), Z0(f) —
exactly what HFSS provides — the closed-form analytical solution
(including RF loss, velocity walk-off, source/load reflections, junction
RC divider) is exact within the quasi-TEM + small-signal assumptions, so
it was chosen. ABCD only pays off for tapered/segmented electrodes or
lumped pad parasitics (possible later extension).

**Key conventions.** R in Ω·m, C in F/m (R·C = junction time constant —
"Ω/m" was a units slip in the original request, resolved); α accepted in
dB/cm, dB/mm, dB/m, Np/m; EO S21 = 20·log10|m(f)/m(0)| in electrical dB;
f₃dB is the first −3 dB crossing referenced to DC; −6 dB electrical =
3 dB optical also reported; e^{+jωt} phasors.

**Both HFSS modes.** `RFLine.includes_junction=True`: HFSS stack-up
contained the doped silicon → γ, Z0 used directly, only the RC divider
H_RC = 1/(1+jωRC) applied on top (not double counting — it converts
electrode voltage to depletion-capacitance voltage). `False`: metal-only
HFSS → model adds Y_j = jωC/(1+jωRC) per unit length and computes the
loaded line (`model.load_line`). Both modes agree to machine precision on
the shared placeholder device (self-tested).

**Verification.** 10 self-tests against closed-form limits (flatness when
matched/lossless/velocity-matched, analytic loss limit, |sinc| walk-off
limit, DC normalization under mismatch, loading trends, mode A ≡ mode B,
passivity, push-pull invariants, dual-drive identity). Placeholder
results: GSG 18.3 GHz, GSGSG differential 17.4 GHz (−3 dB electrical),
both loss-dominated.

## 2. Q&A highlights (condensed)

**"Is a SPICE-style equivalent circuit one of the compared approaches?"**
"Equivalent circuit" means two things. (a) A distributed circuit *solved in
the frequency domain with transmission-line math* — that IS this model.
(b) A lumped RLC ladder run in an actual circuit simulator — approach #4,
rejected for this goal: it needs ~λ/10 segments, frequency-independent
elements (so tabulated α(f), n_RF(f) must be re-fitted with skin-effect
networks), and artificial circuitry for the optical walk-off integral. It
only earns its cost for driver co-simulation, transient/eye/large-signal
work, or PDK handoff. A well-converged ladder and this model agree for a
uniform electrode; this model would be the golden reference if a netlist
export is ever built.

**"Where are the equations and dummy data from?"** Equations re-derived
in-session from standard theory and verified against limiting cases; the
canonical published sources: Gopalakrishnan et al., JLT 12(10) 1994 (full
m(f) with loss, walk-off, reflections — the exact structure implemented);
Alferness 1982 (matched-case factor); Ghione 2009 (textbook treatment);
Yu & Bogaerts JLT 30(11) 2012 (junction loading of the line, silicon);
Patel et al. Opt. Express 23(11) 2015 (silicon TW-MZM model + experiment);
Chrostowski & Hochberg 2015; Pozar (line/network formulas). Citations from
memory — verify page/eq numbers before quoting. The placeholder data is
**invented** (no paper, no real device): per-junction 8 Ω·mm, 240 fF/mm,
n_g = 3.85, CPW-like metal-only curves — each value chosen inside typical
published ranges so the derived quantities (loaded Z0 ≈ 37 Ω, n_RF ≈ 4,
~1–2 dB/mm at 20–30 GHz, ~18 GHz BW at 2.5 mm) land where real silicon
devices do.

**"Does it handle GSG and GSGSG?"** Yes — the engine models one quasi-TEM
mode on one line, so topology enters only through the inputs. Single
junction: per-junction `Junction(R_j, C_j)`. Series push-pull (GSG or
GSGSG): `Junction.series_push_pull(R_j, C_j)` → effective (2R_j, C_j/2);
same RC pole, half the loading — mandatory in mode B, immaterial in mode A.
Differential GSGSG: odd-mode n_RF/α, *differential* Z0 (= 2·Z_odd), and
differential Z_S/Z_L (e.g. 100 Ω), assuming symmetric layout and purely
differential excitation (no common mode / mode conversion — check in
HFSS). Two independent electrodes: `model.combine_dual_drive(res1, res2)`.

## 3. Why a lower load impedance increases bandwidth (kept in detail)

*(Context: to try it, change `z_load` — `inputs.py:253` for the GSG
example (e.g. 50 → 30 Ω), `inputs.py:292` for the differential example
(its loaded Z_diff ≈ 69 Ω, so "low" there means e.g. 45 Ω), or
`inputs.py:218` for the default; in your own scripts just pass
`Drive(length=..., z_src=50.0, z_load=30.0)`. No engine change is needed —
the physics enters via Γ_L at `model.py:127` and the DC normalization at
`model.py:136`.)*

### The one-sentence version

The termination only controls the response near DC; the mid-band response
is set by the source and the line, not the load. Under-terminating *lowers
the DC anchor* that the whole normalized curve is referenced to, so the
(unchanged) mid-band response now sits *above* the reference — a peak —
and the curve takes longer to fall 3 dB below it.

### Three frequency regimes

**Near DC — the termination rules.** With γL → 0 the line is electrically
transparent: every point of the junction sees the plain resistive divider
m(0) = Z_L/(Z_S+Z_L). For the GSG device: 50 Ω gives 0.500·V_s; 30 Ω gives
0.375·V_s. Lower load = *less* modulation voltage at DC, not more.

**Mid-band — the termination disappears.** At higher frequencies the
optical wave effectively interacts only with the forward-traveling wave,
whose launch amplitude is V₊/V_s → Z0/(Z0+Z_S) — for Z0 ≈ 37 Ω against the
50 Ω source, about 0.425·V_s, *independent of Z_L*. The reflected wave
that carries the termination's influence back is killed by two separate
mechanisms, both visible in the formula: its round trip is attenuated by
e^(−2αL), and its contribution to the accumulated optical phase decoheres
because it counter-propagates — the backward overlap factor F_ctr shrinks
like 1/((β+β_o)L) even in a lossless line. (In the code: the
`gl * exp(-x_co) * F(x_ct)` term of m.)

**The transition** happens over the first few GHz, on the scale where the
round-trip phase 2βL reaches ~π — about c/(8·n_RF·L) ≈ 4 GHz for the
2.5 mm, n_RF ≈ 4.1 line — with residual standing-wave ripple of period
c/(2·n_RF·L) ≈ 15 GHz, damped by loss.

Putting the two levels side by side, normalized to DC as the plot does:

- **Z_L = 50 Ω:** anchor 0.500, mid-band 0.425 → mismatch contributes a
  **−1.4 dB shelf**. Reflections *cost* response before loss and RC even
  start.
- **Z_L = 30 Ω:** anchor 0.375, mid-band 0.425 → a **+1.1 dB shelf**. Same
  absolute mid-band voltage, but referenced to a depressed anchor it reads
  as peaking.

The sign of Γ_L = (Z_L−Z0)/(Z_L+Z0) does the work: with Z_L < Z0 the echo
returns *inverted*, cancels line voltage most effectively near DC (where
it returns unattenuated and in phase), and as frequency rises the
cancellation dies, so the response *recovers* upward instead of sagging.

### Equivalent port-side picture

The input impedance of a terminated line is
Z_in = Z0·(Z_L + Z0·tanh γL)/(Z0 + Z_L·tanh γL). At DC, Z_in = Z_L = 30 Ω —
a poor divider against the 50 Ω source. As the line becomes electrically
long and lossy, tanh γL → 1 and Z_in → Z0 = 37 Ω regardless of the
termination: the lossy line "self-matches", and the input voltage division
improves with frequency. Improving division + fixed normalization = rising
normalized response.

### Why that buys bandwidth

The intrinsic roll-offs (junction RC, RF loss, walk-off) are unchanged —
the knob doesn't touch them. But f₃dB is defined relative to the
DC-referenced 0 dB line, so the +1.1 dB shelf is pre-tilt: the roll-offs
must burn through ~4.1 dB (from +1.1 to −3) instead of ~1.6 dB (from −1.4
to −3). The full model's slope around 18 GHz is roughly 0.18 dB/GHz
(−3 dB at 18.3 GHz, −6 dB at 35.1 GHz), so ~2.5 dB of extra headroom
should push f₃dB from ~18 GHz into the high-20s (back-of-envelope; run it).
The `bandwidth_breakdown_gsg.png` "Mismatch only" curve shows the shelf
flipping from below 0 dB to above it — the whole story in one trace; the
standing-wave change is visible in `voltage_profile_gsg.png`.

### What it costs — the part the normalized plot hides

This is equalization, not amplification. Every dB of peaking is paid for
by lowering the absolute low-frequency response: 20·log10(0.375/0.500) ≈
−2.5 dB of junction voltage at DC — a worse effective low-frequency Vπ and
a link-budget penalty. Pushed further it also gives: a visibly non-flat
response (overshoot/ISI in eye diagrams), worse electrical S11 at the
modulator input at low frequencies (the driver looks into ~30 Ω), and more
standing-wave ripple before loss damps it. It cannot rescue a device
limited by junction RC or RF loss — it only spends the mismatch budget
more cleverly. The same physics is why low-output-impedance CMOS drivers
(source-side under-matching, Γ_S < 0) produce similar peaking; some
designs use both ends together.

### Definition caveat

The bandwidth gain exists under the DC-referenced −3 dB definition (what
`model.find_crossing` implements). If a paper defines bandwidth relative
to the response *peak*, the numbers are not comparable — check which
reference the authors used. (A peak-referenced option in `find_crossing`
would be a small change if ever needed.)

Quick sweep snippet:

```python
for zl in (50, 40, 30, 25):
    res = model.eo_response(rf, junction, optical,
                            Drive(2.5e-3, 50.0, zl))
    # overlay res.s21_db; res.f3db per value
```

## 4. Open threads / next steps

- Replace `inputs.example_inputs*()` placeholders with real CHARGE / MODE
  / HFSS exports (the entire point of the model).
- Run the z_load under-termination experiment (30 Ω GSG; ~45 Ω GSGSG).
- Optional: peak-referenced bandwidth variant of `find_crossing`.
- Optional: CSV loaders matched to the actual HFSS export format.
- Optional: bias/length sweep helpers.
- Optional: ABCD layer (pad/wirebond parasitics, tapered or segmented
  electrodes).
- SPICE export: deferred, not rejected (position revised at session end).
  Original exclusion was scope discipline — for small-signal EO S21 a
  ladder adds discretization + RLGC-fitting error and converges to what
  the closed form already gives. Triggers to build it: driver
  co-simulation (real output impedance interacts with the same peaking
  physics as under-termination), transient/eye/large-signal C(V) studies,
  or EDA/PDK handoff. Path: generate the ladder netlist from the same
  CHARGE/MODE/HFSS data and validate its AC EO response against this
  model before trusting transients (~a day of scoped work).
