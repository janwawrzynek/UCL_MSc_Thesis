# Sterile Neutrino Decays with an Exotic Axion-Like Extension

**MSc Physics Research Thesis — University College London, 2025**
*Supervisors: Professor Frank Deppisch & Professor Robert Thorne*

> *"The work significantly exceeds the requirements for an MSc-level project."* — Project Supervisor
> 
> **Final Grade: Distinction (75.35%)**

---

## Overview

This repository contains the full Python simulation engine developed for my MSc thesis at UCL's High Energy Physics Group. The engine models the decay properties and cosmological abundances of **Heavy Neutral Leptons (HNLs)** — exotic Majorana neutrinos that are compelling candidates for resolving three major Beyond Standard Model (BSM) puzzles:

- **Neutrino oscillations** — via the Type-I Seesaw Mechanism
- **Baryon asymmetry of the universe** — through CP-violating leptogenesis
- **Dark matter** — as a viable cold/warm dark matter candidate

The central physics extension explored here is the model proposed by [Deppisch et al. (2025)](https://doi.org/10.1088/1475-7516/2025/02/054), which introduces an **Axion-Like Particle (ALP)** as a dominant HNL decay channel: **N → aν**. This exotic channel allows HNLs to decay fast enough to evade constraints from Big Bang Nucleosynthesis (BBN), opening up previously excluded regions of the HNL parameter space.

---

## Physics Background

### The Problem: BBN Constraints on HNLs

Without a new exotic decay channel, HNLs with lifetimes τ > 1s would decay dominantly into hadronic channels. The resulting mesons would perturb the proton–neutron conversion ratio in the primordial plasma, altering <sup>4</sup>He abundances and spoiling BBN predictions. This places a strong lower bound on allowed HNL masses.

### The Solution: ALP-Mediated Decays

By introducing a very light pseudo-scalar ALP that couples to the HNL, a new dominant two-body decay channel opens up:

$$\sum_{\alpha=e,\mu,\tau} \Gamma_{N_i \to a\nu_\alpha} = \frac{|U_{\ell\alpha N_i}|^2 m^3_{N_i}}{4\pi f_a^2} \sqrt{1 + \left(\frac{m_a}{m_{N_i}}\right)^2} \left[1 - \left(\frac{m_a}{m_{N_i}}\right)^2\right]^{3/2}$$

This causes HNLs to decay faster, circumventing the BBN bound for lower HNL masses — provided the ALP itself is long-lived (τ_a > 10⁴ s), which is satisfied for m_a ≲ 1 keV where the dominant decay is the suppressed a → νν channel.

---

## Key Results

### HNL Iso-Lifetime Contours
Computed iso-lifetime curves in the (m_N, |U_ℓN|²) plane for all three sterile neutrino flavours (N₁, N₂, N₃) at τ = 0.1s, 1.0s, and 10.0s. A potentially viable region of parameter space was identified for the N₂ HNL at low masses (~0.1 GeV), for exclusive mixing with the μ lepton flavour.

### ALP vs SM Branching Ratios
For active-sterile mixing |U_ℓN|² = 10⁻¹⁰:
| ALP Decay Constant f_a | ALP-Dominant Mass Range |
|---|---|
| 10⁵ GeV | m_N < 5 GeV |
| 10⁶ GeV | m_N ≲ 0.6 GeV (equal at ~0.6 GeV) |
| 10⁷ GeV | m_N ≲ 0.2 GeV |

### Thermally Averaged Decay Rates
Computed thermally averaged decay width densities γ(z) for both SM and ALP channels across all three HNL flavours, providing a foundation for full BBN analysis via coupled Boltzmann equations.

---

## Engine Architecture

The codebase is designed with **modularity and extensibility** as core priorities. New exotic particles and decay channels can be incorporated by subclassing the existing particle hierarchy, with no changes required to the simulation or plotting layers.

```
Model Parameters ──┐
                   ├──► Model.py ──► Particles.py ──► Plotting.py ──► Simulation.py ──► CLI
Rel. DOF (CSV) ────┘
```

### File Structure

```
.
├── OOP/
│   ├── ModelParameters.py       # Physical constants, mixing matrices, ALP & HNL parameters
│   ├── RelDegreesOfFreedom.csv  # g*(T) data from lattice QCD (Table 8 in thesis)
│   ├── RelDegreesOfFreedom.py   # Interpolated g*(T) function in terms of z = m_N/T
│   ├── Model.py                 # Central model class: QCD running, cosmology, kinematics
│   ├── Particles.py             # Full OOP particle hierarchy with all decay width implementations
│   ├── Plotting.py              # High-level plotting routines for all physics scenarios
│   └── Simulation.py           # CLI execution layer; task dispatch and output management
└── README.md
```

---

## Particle Class Hierarchy

The particle inheritance system (defined in `Particles.py`) reflects the physical relationships between particles and enforces type safety — decay width functions cannot accept physically invalid particle classes as inputs.

```
Particle
├── Quark          (up, down, charm, strange, top, bottom)
├── Pseudoscalar   └── ALP
├── Lepton
│   ├── ChargedLepton  (Electron, Muon, Tau)
│   └── Neutrino
│       ├── LightNeutrino   (νe, νμ, ντ)  — mass set via seesaw: m_ν = |U|²·m_N
│       └── SterileNeutrino (N1, N2, N3)  — full decay width implementation
└── Hadron
    └── Meson
        ├── PseudoscalarMeson
        │   ├── ChargedPseudoscalarMeson  (π±, K±, D±, Ds±, B±, Bc)
        │   └── NeutralPseudoscalarMeson  (π⁰, η, η′, ηc)
        └── VectorMeson
            ├── ChargedVectorMeson        (ρ±, D*±, Ds*±)
            └── NeutralVectorMeson        (ρ⁰, ω, ϕ, J/ψ)
```

---

## Decay Channels Implemented

### Standard Model Channels
| Channel | Mediator | Mass Regime |
|---|---|---|
| N → ℓ⁻α νβ ℓ⁺β (α ≠ β) | Charged Current (W±) | All |
| N → να ℓ⁻β ℓ⁺β | NC + CC interference | All |
| N → να νβ ν̄β | Neutral Current (Z) | All |
| N → ℓ⁻α h⁺_P (π±, K±, D±, ...) | CC, exclusive hadronic | m_N < 2 GeV |
| N → να h⁰_P (π⁰, η, η′, ηc) | NC, exclusive hadronic | m_N < 2 GeV |
| N → ℓ⁻α h⁺_V (ρ±, D*±, ...) | CC, exclusive hadronic | m_N < 2 GeV |
| N → να h⁰_V (ρ⁰, ω, ϕ, J/ψ) | NC, exclusive hadronic | m_N < 2 GeV |
| N → ℓ⁻α ui d̄j + Δ_QCD | CC, inclusive quark | m_N > 2 GeV |
| N → να qq̄ + Δ_QCD | NC, inclusive quark | m_N > 2 GeV |

### Exotic ALP Channel (Deppisch et al. 2025)
| Channel | Description |
|---|---|
| N → a να | Dominant exotic decay; opens BBN-safe parameter space |
| a → νν | ALP decay (tree-level dominant for m_a ≲ 1 keV) |

The hadronic regime switch at 2 GeV uses a one-loop running strong coupling with quark flavour thresholds at m_charm = 1.3 GeV and m_bottom = 4.5 GeV:

$$\alpha_s(m_N^2) = \frac{1}{b_0 \ln(m_N^2 / \Lambda^2)}, \quad b_0 = \frac{33 - 2n_f}{12\pi}$$

---

## Installation

**Requirements:** Python 3.10+

```bash
git clone https://github.com/janwawrzynek/UCL_MSc_Thesis.git
cd UCL_MSc_Thesis
pip install numpy scipy pandas matplotlib
```

---

## Usage

All simulations are run from the command line. Set the `TASK_TO_RUN` variable in `Simulation.py` to the desired task, then execute:

```bash
python -m OOP.Simulation
```

### Available Tasks

| Task String | Description |
|---|---|
| `SINGLE_CALCULATION` | Partial width, total width, and branching ratio for a single user-specified channel |
| `PLOT_BRANCHING_RATIOS_HADRONS` | Exclusive hadronic branching ratios vs m_N (< 1 GeV) |
| `PLOT_BRANCHING_RATIOS_QUARKS` | Inclusive quark-level + ALP branching ratios vs m_N |
| `PLOT_QCD_CORRECTION_VS_MASS` | Δ_QCD correction factor as a function of m_N |
| `PLOT_REL_DEGREES_OF_FREEDOM_VS_Z` | Relativistic degrees of freedom g*(z) |
| `PLOT_THERMALLY_AVERAGED_WIDTHS` | Thermally averaged decay rate densities γ(z) |
| `LIFETIME_CONTOUR_PLOT` | HNL iso-lifetime contours in the (m_N, \|U_ℓN\|²) plane |
| `PLOT_COMBINED_BRANCHING_RATIO_SM_ALP` | Combined SM vs ALP branching ratios |
| `PLOT_FERMI_THEORY_ERROR` | Relative error of Fermi contact approximation vs full propagator |

### Example: Computing a Single Decay Channel

In `Simulation.py`, configure the incoming and outgoing particles, then set:
```python
TASK_TO_RUN = "SINGLE_CALCULATION"
```
The engine will print the partial decay width (GeV), total decay width (GeV), lifetime (s), and branching ratio for the specified channel.

### Example: Scanning HNL Parameter Space

To reproduce the iso-lifetime contour plots from the thesis:
```python
TASK_TO_RUN = "LIFETIME_CONTOUR_PLOT"
```
Adjust `mN_min`, `mN_max`, and target lifetimes in the Plotting class as needed.

---

## Configuring Physical Parameters

All physical inputs are centralised in `ModelParameters.py`. To scan over a parameter, simply assign new values before constructing the `Model` instance — the simulation layer never mutates global state, so all runs are fully reproducible.

```python
from OOP.ModelParameters import ModelParameters
from OOP.Model import Model

params = ModelParameters()

# Customise scenario
params.mN = 0.5          # HNL mass in GeV
params.ma = 1e-6         # ALP mass in GeV
params.fa = 1e6          # ALP decay constant in GeV
params.U_matrix[1][1] = 1e-10  # |U_μN2|²

model = Model(params, csv_path="OOP/RelDegreesOfFreedom.csv")
```

**Active-sterile mixing matrix** (rows = HNL flavours N1, N2, N3; columns = active flavours e, μ, τ):
```python
params.U_matrix = [
    [|UeN1|², |UμN1|², |UτN1|²],
    [|UeN2|², |UμN2|², |UτN2|²],
    [|UeN3|², |UμN3|², |UτN3|²]
]
```

---

## Validation & Benchmarks

All SM branching ratios were benchmarked against published results from [Bondarenko & Boyarsky (2018)](https://doi.org/10.1007/JHEP11(2018)032) prior to adding the ALP extension. Key checks:

- ✅ N → νν̄ν invisible channel dominates at low m_N
- ✅ Pion channel activates sharply at m_N = m_π± = 139.6 MeV and becomes dominant
- ✅ η, ρ, K channels activate at their respective mass thresholds
- ✅ Hadronic decays dominate over leptonic above ~1 GeV
- ✅ Fermi contact approximation error scales as (m_N/m_W)⁴ — negligible below m_W

---

## Cosmological Extension

The framework provides thermally averaged decay rate densities:

$$\gamma_{X \to Y} = n_X^{eq}(z) \frac{K_1(z)}{K_2(z)} \Gamma_X, \quad z = \frac{m_N}{T}$$

This provides the foundation for solving the full set of coupled Boltzmann equations for HNL and ALP comoving yields, which would allow direct constraints on (m_N, |U_ℓN|²) from BBN observables. The four coupled equations are derived in Section 7 of the thesis.

---

## Citation

If you use this code in your research, please cite the thesis and the key references it builds upon:

```bibtex
@mastersthesis{wawrzynek2025hnl,
  author  = {Jan Wawrzynek},
  title   = {Sterile Neutrino Decays with an Exotic Axion-Like Extension},
  school  = {University College London},
  year    = {2025},
  note    = {MSc Physics, Distinction}
}

@article{deppisch2025relaxing,
  author  = {Deppisch, Frank F. and Gonzalo},
  title   = {Relaxing limits from Big Bang Nucleosynthesis on Heavy Neutral Leptons with Axion-Like Particles},
  journal = {JCAP},
  volume  = {2025},
  number  = {02},
  pages   = {054},
  year    = {2025}
}

@article{bondarenko2018phenomenology,
  author  = {Bondarenko, K. and Boyarsky, A. and others},
  title   = {Phenomenology of GeV-scale Heavy Neutral Leptons},
  journal = {JHEP},
  volume  = {2018},
  number  = {11},
  year    = {2018}
}
```

---

## Acknowledgements

I would like to thank Professor Frank Deppisch and Professor Robert Thorne for their supervision and guidance throughout this project. This work builds directly on [Deppisch et al. (2025)](https://doi.org/10.1088/1475-7516/2025/02/054) and [Bondarenko & Boyarsky (2018)](https://doi.org/10.1007/JHEP11(2018)032).

---

## License

MIT License — see `LICENSE` for details.
