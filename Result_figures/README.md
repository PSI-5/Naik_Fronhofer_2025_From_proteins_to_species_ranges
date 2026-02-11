# Thermal Adaptation Models - Analysis Results

This repository contains analysis scripts and results for thermal adaptation models during range expansion. The scripts reproduces results for Naik et al. 2025 from simulation data.


## Repository Structure
├── Enzyme_catalysis_model/ # Enzyme catalysis model data
│ ├── eco/ # Ecological simulations (without adaptation)
│ └── evo/ # Evolutionary simulations (with adaptation)
├── Gauss_model/ # Gaussian model data
│ ├── eco/
│ └── evo/
├── Protein_denat_model/ # Protein denaturation model data
│ ├── eco/
│ └── evo/
├── Result_figures/ # Analysis scripts and generated figures
│ ├── plot_fig_2.py # Figure 2: Range expansion dynamics
│ ├── plot_fig_3.py # Figure 3: Trait evolution and niches
│ ├── plot_fig_4.py # Figure 4: Mutation effect distributions
│ └── Paper_figs/ # Generated figure outputs
├── install_dependencies.sh # Dependency installation script
└── README.md # This file


## Available Data

Each model directory contains two types of simulation results:

### 1. **Ecological Simulations** (`eco/`)
- Simulations WITHOUT thermal adaptation
- Fixed traits throughout simulation
- Baseline for comparison with evolutionary runs

### 2. **Evolutionary Simulations** (`evo/`)
- Simulations WITH thermal adaptation
- Traits can evolve via mutation and selection
- Shows adaptation to temperature gradient

### Data Files in Each Model Directory:

**Aggregated Statistics:**
- `pf.csv` - Range front positions over time (with confidence intervals)
- `growth_pf.csv` - Range front Population growth rate statistics
- `mut_eff_pf.csv` - Range front Beneficial mutation effect statistics  
- `niche_pf.csv` - Range front Thermal niche characteristics
- Model-specific trait files (e.g., `DHm_pf.csv`, `dGr_pf.csv`, `sigma_pf.csv`)

**Note**: Plot scripts do not plot all the data provided in the data files. 

**Individual Simulation Runs:**
- `pf_0.csv` to `pf_79.csv` - 80 replicate simulations

> **Note**: The individual mutation data files used by `plot_fig_4.py` are not included in this repository due to their large size. Motivated users can generate these files by running the simulation codes available in the parent directory (outside `Result_figures/`). See the main README in the parent directory for simulation instructions.

---

## Model Descriptions

### 1. Protein_denat_model

A **protein denaturation–based reaction rate model**. Growth and performance are governed by temperature-driven protein stability and denaturation dynamics.

### 2. Enzyme_catalysis_model

An **enzyme catalysis–based temperature-dependent reaction rate model**, where performance depends on enzymatic reaction kinetics and thermal sensitivity.

### 3. Gauss_model

A **phenomenological Gaussian thermal performance curve** model, commonly used as a simplified representation of temperature-dependent fitness.

---
# Script Details

### `plot_fig_2.py`

**Purpose**: Generates Figure 2 of the manuscript comparing range expansion dynamics across models with and without thermal adaptation.

**Key Outputs:**
- Panel A: Range front positions without adaptation (ecological)
- Panel B: Range front positions with adaptation (evolutionary)
- Panel C: Initial expansion speeds without adaptation
- Panel D: Initial expansion speeds with adaptation

**Key Analyses**:
- Range front progression over time (2000-4000 time steps)
- Linear regression to calculate expansion speeds (time < 2200 steps)
- Comparison of hot vs cold front dynamics
- Ecological (no adaptation) vs evolutionary (with adaptation) scenarios

**Usage:**
```bash
python plot_fig_2.py
```

### `plot_fig_3.py`

**Purpose**: Generates Figure 3 of the manuscript analysing trait evolution and thermal niche adaptation.

**Key Outputs** (Figure 3):
- **Panels A-C**: Equilibrium thermal niches for each model

**Key Analyses**:
- Thermal niches at equilibrium (final time step) at the cold and hot front

**Usage:**
```bash
python plot_fig_3.py
```

### `plot_fig_4.py`

**Purpose**: Generates Figure 3 of the manuscript visualising spatial distribution of mutation effects across the temperature gradient.

**Outputs** (Figure 4):
- **Panels A-C**: Mutation effects by spatial position for each model

**Key Analyses**:
- Mutation effects colored by spatial position (temperature)
- Comparison of mutation effect distributions across models

**Usage:**
```bash
python plot_fig_4.py
```

**Important**: This script requires individual mutation data files (`individual_*.csv`) which are not included in this repository due to size constraints. To generate these files, run the simulation codes available in the parent directory following the instructions in the main README.

## Citation

If you use this code or data in your research, please cite:

**Naik, S. & Fronhofer, E.** (2024). From proteins to species ranges: a framework for understanding thermal adaptation during range expansions. *Proceedings of the Royal Society B*. DOI: 10.1098/rspb
