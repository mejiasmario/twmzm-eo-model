# TW-MZM electro-optic response model

Predicts the small-signal **EO S21** and **3 dB bandwidth** of a silicon
traveling-wave Mach-Zehnder modulator by combining parameters extracted from
physical simulators. It does **not** replace them:

| Tool            | Provides                                   | Model class (`inputs.py`) |
|-----------------|--------------------------------------------|---------------------------|
| Lumerical CHARGE| junction `R` [Ω·m], `C` [F/m] at bias      | `Junction`                |
| Lumerical MODE  | optical group index `n_g`                  | `Optical`                 |
| Ansys HFSS      | `f`, `n_RF(f)`, `Z0(f)`, `α(f)`            | `RFLine`                  |
| design choice   | length `L`, source/termination impedances  | `Drive`                   |

## Method

Analytical traveling-wave theory: the closed-form solution of the telegrapher
equations for a uniform electrode, including RF loss, optical/RF velocity
mismatch, source/load reflections (standing waves), and the junction RC
voltage divider.

```
γ(f)   = α_Np(f) + j·2πf·n_RF(f)/c                     (HFSS)
β_o(f) = 2πf·n_g/c                                     (MODE)
H_RC   = 1/(1 + j·2πf·R·C)                             (CHARGE)

Γ_L = (Z_L − Z0)/(Z_L + Z0),   Γ_S = (Z_S − Z0)/(Z_S + Z0)
V₊/V_s = Z0 / [(Z0 + Z_S)(1 − Γ_S·Γ_L·e^(−2γL))]

F(x)  = (1 − e^(−x))/x
x_co  = (γ − jβ_o)·L            co-propagating term
x_ct  = (γ + jβ_o)·L            counter-propagating (reflected) term

m(f)  = H_RC · (V₊/V_s) · [F(x_co) + Γ_L·e^(−x_co)·F(x_ct)]
EO S21(f) = 20·log10 |m(f)/m(0)|,   m(0) = Z_L/(Z_S + Z_L)
```

### Two HFSS interpretation modes (`RFLine.includes_junction`)

* **`True` (loaded line):** the HFSS stack-up contained the doped silicon, so
  `γ`, `Z0` already describe the loaded electrode. The model applies only the
  junction RC divider `H_RC` (this converts electrode voltage to depletion-
  capacitance voltage and is *not* double counting the loading).
* **`False` (metal-only):** the model first adds the distributed junction
  admittance `Y_j = jωC/(1 + jωRC)` to the line (`model.load_line`), computes
  the loaded `γ`, `Z0`, then proceeds identically.

### Electrode topologies (GSG, GSGSG, dual drive)

The engine models **one quasi-TEM mode on one line with one distributed R–C
load**, so topology enters only through the values you feed it. The fixed
×2 / ×½ topology factors move Vπ (out of scope) and cancel in the
normalized EO S21.

| Topology | Line data (HFSS) | `Junction` | `Drive` Z_S/Z_L |
|---|---|---|---|
| GSG, one junction per cell | driven CPW mode | `Junction(R_j, C_j)` per-junction | single-ended, 50 Ω typ. |
| Series push-pull (GSG or GSGSG), incl. differential drive | driven mode; for differential: odd mode with **differential** Z0 (= 2·Z_odd, as mixed-mode extraction reports) | `Junction.series_push_pull(R_j, C_j)` → effective (2R_j, C_j/2) | differential system, 100 Ω typ. |
| Two independent electrodes (dual drive) | per arm | per arm | per arm → `model.combine_dual_drive(res1, res2)` |

Mode A (`includes_junction=True`) is insensitive to the per-junction vs
effective distinction — only τ = R·C enters, and the series combination
leaves it unchanged. In **mode B the effective values are mandatory**:
per-junction values would double the loading admittance. Differential
modeling additionally assumes a symmetric layout and purely differential
excitation; common-mode excitation and mode conversion are not modeled —
check those in HFSS.

## Quick start

```bash
python run_example.py          # self-tests, summary, figures into ./figures
python run_example.py --show   # additionally open the figure windows
```

Plug in your own data:

```python
import numpy as np
import model
from inputs import RFLine, Junction, Optical, Drive

rf = RFLine.from_hfss(f_hz, n_rf, z0_ohm, alpha, alpha_unit="dB/cm",
                      includes_junction=True)      # or False (metal-only)
junction = Junction(r_ohm_m=8.0e-3,                # 8 Ω·mm
                    c_f_m=2.4e-10)                 # 240 fF/mm
optical = Optical(n_g=3.85)
drive = Drive(length=2.5e-3, z_src=50.0, z_load=50.0)

res = model.eo_response(rf, junction, optical, drive)
print(res.f3db / 1e9, "GHz")                       # -3 dB electrical bandwidth
```

`res` also carries `f6db` (= 3 dB **optical** bandwidth), the complex
response `m(f)`, `s21_db`, and the loaded `line` actually used. Diagnostics:
`model.response_breakdown` (which limiter dominates), `model.line_sparams`
(electrode S11/S21), `model.voltage_profile` (standing waves).

For a differential GSGSG series push-pull device, the same call with
differential-mode quantities:

```python
junction = Junction.series_push_pull(r_j_ohm_m=8.0e-3,  # per-junction values
                                     c_j_f_m=2.4e-10)
drive = Drive(length=2.5e-3, z_src=100.0, z_load=100.0)  # differential system
rf = RFLine.from_hfss(f_hz, n_rf_odd, z0_differential, alpha_odd, "dB/cm",
                      includes_junction=False)
res = model.eo_response(rf, junction, optical, drive)
```

For dual independent electrodes: `model.combine_dual_drive(res_arm1, res_arm2)`.

## Units

| Quantity | Unit    | Example                              |
|----------|---------|--------------------------------------|
| `R`      | Ω·m     | 10 Ω·mm → `1.0e-2` (R×C = time const.)|
| `C`      | F/m     | 250 fF/mm → `2.5e-10`                |
| `α`      | any of `dB/cm`, `dB/mm`, `dB/m`, `Np/m` via `alpha_unit` |
| `f`, `L` | Hz, m   |                                      |

`Junction.from_totals(R_total, C_total, length)` converts total values from a
CHARGE slice of known length.

## Assumptions

1. Single quasi-TEM RF mode; HFSS values are valid de-embedded uniform-line
   parameters over the band (periodic loading is fine if pitch ≪ λ_RF).
2. Small-signal, fixed bias; `R`, `C` are the values at that bias.
3. Phase efficiency dφ/dV_junction is frequency-flat (no carrier dynamics);
   it cancels in the normalized response. Vπ/chirp are out of scope.
4. Modulation is proportional to the voltage across the depletion
   capacitance, hence `H_RC`.
5. Optical envelope travels at `c/n_g`; `n_g` dispersion neglected in-band;
   optical loss drops out of the normalization.
6. Source/termination are lumped impedances at the electrode ends; pad or
   wire-bond parasitics are not included.
7. Conventions: `e^{+jωt}` phasors; EO S21 in electrical dB; the −3 dB
   electrical point (= 1.5 dB optical) is "the" EO bandwidth; −6 dB
   electrical = 3 dB optical is reported alongside.

## Files

| File             | Role                                            |
|------------------|-------------------------------------------------|
| `inputs.py`      | dataclasses, unit converters, placeholder data  |
| `model.py`       | analytical engine, diagnostics, self-tests      |
| `plotting.py`    | figure routines                                 |
| `run_example.py` | end-to-end demo (self-tests → summary → plots)  |

Two placeholder devices are provided, both fictitious but plausible
(per-junction 8 Ω·mm, 240 fF/mm, n_g = 3.85, L = 2.5 mm):
`inputs.example_inputs()` — single-ended GSG, 50 Ω system — and
`inputs.example_inputs_differential()` — GSGSG series push-pull driven
differentially, 100 Ω differential system, junction via
`Junction.series_push_pull`. For each, both HFSS modes describe the *same*
device, so their EO responses must agree — checked at runtime; figures are
written with `_gsg` / `_gsgsg` suffixes.

## References

Standard treatments this implementation follows: G. Ghione,
*Semiconductor Devices for High-Speed Optoelectronics* (traveling-wave
modulator theory); H. Yu & W. Bogaerts, *JLT* **30**, 2012 (equivalent
circuit of loaded TW electrodes); D. Patel *et al.*, *Opt. Express* **23**,
2015 (silicon TW-MZM model + experiment); L. Chrostowski & M. Hochberg,
*Silicon Photonics Design*, ch. 6.
