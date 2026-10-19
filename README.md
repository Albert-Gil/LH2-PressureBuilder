# LH₂ Pressure Builder

**A MATLAB reduced-order model of natural circulation through a liquid-hydrogen pressure-builder vaporizer.**

[![MATLAB](https://img.shields.io/badge/MATLAB-tested%20R2026a-blue?logo=mathworks)](https://www.mathworks.com/products/matlab.html)
[![CoolProp](https://img.shields.io/badge/CoolProp-7.2.0-green)](https://coolprop.org/)
[![License](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

Liquid leaves the bottom of a horizontal cylindrical LH₂ tank, absorbs heat in a horizontal heater, and returns through an unheated vertical pipe into the ullage. The model predicts circulating flow and outlet state by balancing liquid supply head against return-column head, friction, fittings, and acceleration losses.

**Scope:** quasi-steady vaporizer-loop model at imposed tank pressure and liquid temperature. It does **not** calculate tank pressure rise, ullage heating, condensation, or a closed-tank pressurization transient.

<p align="center">
<img src="results/figures/mdot_map_H_Q.png" width="48%" alt="Circulation over height and absorbed heater power">
<img src="results/figures/SOC_scenario_mdot_Tout.png" width="48%" alt="Flow and outlet temperature during prescribed drawdown">
</p>

Blank map regions are cases without an accepted solution (for example above the 250 K outlet cap). They are **not** zero flow.

![Annotated pressure-builder model schematic](docs/model_schematic.svg)

---

## Reduced-order physics

At fixed loop pressure \(P\) and liquid inlet temperature \(T_L\):

| Balance | Expression |
|--------|------------|
| Energy | \(h_{out} = h_L + \dot{Q}/\dot{m}\) |
| Supply head | \(\rho_L(P,T_L)\, g\, H\) |
| Return head | \(\rho_{out}(P,h_{out})\, g\, H_{elev}\) |
| Losses | Effective \(K_L\), \(K_g\), Darcy friction in supply/return, homogeneous acceleration |
| Closure | Solve supply − return − losses = 0 for \(\dot{m}\) |

The horizontal heater adds enthalpy but **no vertical head**. The return is treated as **adiabatic** at constant pressure, so density and viscosity equal the heater outlet. Two-phase properties use a homogeneous no-slip mix, continuous at \(x = 0\) and \(x = 1\). Friction blends laminar and Haaland turbulence over Re = 2300–4000.

\(K_L\) and \(K_g\) lump fittings and unknown heater-passage resistance. Plate geometry is not resolved. Failed solves report `ok = false` and `mdot = NaN` (never a fabricated zero). Outlet temperature is capped by `T_out_max` (default **250 K**).

### Default design

| Parameter | Default |
|-----------|---------|
| Pipe ID `D` | 35.92 mm |
| Loss coefficients | \(K_L = 6.565\), \(K_g = 5.252\) |
| Inlet | 8.273 bar, 25 K |
| Reference | H = 1.54 m, \(\dot{Q}\) = 21.6 kW |
| Sweep | H = 0.01–2 m × \(\dot{Q}\) = 0.36–36 kW |
| Scenarios | 28.8 → 28.8…10.8 kW over 90 min |
| Power decay | \(\dot{Q}_0\) = 42 kW, τ = 15 min |

---

## Run with MATLAB

### Prerequisites

1. Link MATLAB to a Python that has CoolProp (once):

```matlab
pyenv('Version', 'C:\Path\To\python.exe');
```

2. Install Python packages:

```bash
pip install -r requirements.txt
```

### Run the model

Open MATLAB in the **project root** and run:

```matlab
Run_Vaporizer
```

[`Run_Vaporizer.m`](Run_Vaporizer.m) is the single entry point. It adds `model/`, `geometry/`, `studies/`, and `output/` to the path, builds one parameter struct `p`, and runs the studies you enable.

### Where to change dimensions, parameters, and powers

| What | Where |
|------|--------|
| Defaults (geometry, \(K\), pressures, powers, grids, `T_out_max`) | [`model/vap_default_params.m`](model/vap_default_params.m) |
| One-off overrides for this run | After `p = vap_default_params();` in [`Run_Vaporizer.m`](Run_Vaporizer.m) |
| Which studies run | `RUN.*` switches at the top of `Run_Vaporizer.m` |
| MATLAB figures on/off, output folder | `makePlots`, `outputDir` in `Run_Vaporizer.m` |

Example overrides in `Run_Vaporizer.m`:

```matlab
p = vap_default_params();
p.D = 0.04;                              % pipe ID [m]; area recomputed
p.K_g = 5.0;                             % return loss
p.Qdot_list = linspace(1,20,40)*1e3;     % sweep powers [W]
p.scenario_Q0 = 20e3;                    % scenario start [W]
p.scenario_Q_end = [20,15,10]*1e3;
p.T_out_max = 250;                       % outlet cap [K]
```

### Where to find results

| Output | Location |
|--------|----------|
| Datasets | `results/*.mat` (`sweep`, `power_decay`, `SOC_scenario`, `constH_scenario`) |
| Figures | `results/figures/*.png` (when `makePlots = true`) |
| Workspace | struct `results` after `Run_Vaporizer` |

---

## Run with Python

Python does **not** solve the thermosiphon. CoolProp is used from MATLAB for fluid properties. The Python tools regenerate figures (and the schematic) from saved results — useful for headless or batch plotting.

```bash
pip install -r requirements.txt

# Regenerate study figures from results/*.mat
python python/plot_vaporizer_results.py

# Another results folder
python python/plot_vaporizer_results.py --directory path/to/results

# Regenerate the annotated schematic
python python/generate_model_schematic.py
```

You can also set `RUN.python_plots = true` in `Run_Vaporizer.m` so MATLAB calls the figure plotter after the studies.

---

## Layout

```
Run_Vaporizer.m     entry point (switches + optional overrides)
model/              solver, CoolProp wrapper, defaults
geometry/           tank height ↔ volume
studies/            sweep, scenarios, power decay, calibration
output/             MATLAB plotting and result paths
python/             batch figures and schematic
results/            .mat datasets and figures/
docs/               model schematic (SVG/PNG)
tests/              regression checks
```

---

## Limitations

- Adiabatic return; significant pipe heat leak needs a heat-leak model.
- Heater passage resistance is lumped into effective \(K\) losses.
- Homogeneous no-slip two-phase flow; slip and stability are not modeled.
- Quasi-steady only; startup and wall storage are not resolved.
- Numerical convergence is consistency of the equations, not experimental validation.

## Author and license

**Albert Gil**, UC Irvine — [LinkedIn](https://www.linkedin.com/in/albertgilesmendia/).

**GNU GPL v3 or later** — see [`LICENSE`](LICENSE). Open a GitHub issue for questions or bugs.
