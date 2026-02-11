# Temperature-Dependent Performance Curve Models

This repository contains three alternative **mechanistic models** for studying **temperature-dependent performance curves (TPCs)** and **population dynamics under climate change**.
Each model is implemented in **C++** for computational efficiency and executed through a **Python–MPI wrapper**, enabling large-scale parallel simulation runs.

---

## Repository Structure

```
new_submit_code/
├── Enzyme_catalysis_model/
│   ├── run_ec_tpc.cpp
│   ├── run_ec_tpc_cold_drate.cpp
│   └── ec_eg/
│       ├── parameters.in
│       └── run_pll.py
│
├── Gauss_model/
│   ├── run_gss_tpc.cpp
│   ├── run_gss_tpc_cold_drate.cpp
│   └── gss_eg/
│       ├── parameters.in
│       └── run_pll.py
│
├── Protein_denat_model/
│   ├── run_pd_tpc.cpp
│   ├── run_pd_tpc_cold_drate.cpp
│   ├── prot_make_trait_pf_pll.py
│   └── prot_eg/
│       ├── parameters.in
│       └── run_pll.py
```

Each model directory contains:

* **C++ source files** implementing the temperature-dependent growth model
* An **example directory (`*_eg`)** with a parameter file and Python–MPI launcher script

---

## Model Descriptions

### 1. Protein_denat_model

A **protein denaturation–based reaction rate model**. Growth and performance are governed by temperature-driven protein stability and denaturation dynamics.

### 2. Enzyme_catalysis_model

An **enzyme catalysis–based temperature-dependent reaction rate model**, where performance depends on enzymatic reaction kinetics and thermal sensitivity.

### 3. Gauss_model

A **phenomenological Gaussian thermal performance curve** model, commonly used as a simplified representation of temperature-dependent fitness.

---

## Model Variants

Each model provides two alternative C++ implementations:

* `*_tpc.cpp`
  Baseline temperature-dependent performance curve formulation.

* `*_tpc_cold_drate.cpp`
  Alternative formulation incorporating **cold-biased mortality (death rate)**.

These variants allow direct comparison between standard TPCs and scenarios with asymmetric thermal stress.

---

## Installing Dependencies

From the command line, run:

```bash
./install_dependencies.sh
```

This script installs all required dependencies for compiling the C++ code and running the Python–MPI interface.

---

## Running the Examples

Each model directory contains a ready-to-run example inside its `*_eg` folder. Each model runs for 5000 time steps unless it reaches the end faster.

### Step-by-step Instructions

1. **Navigate to an example directory** (e.g., Gaussian model):

```bash
cd Gauss_model/gss_eg
```

2. **Edit simulation parameters** in `parameters.in`:

```bash
vi parameters.in
```

3. **Run the simulation using MPI**:

```bash
mpiexec -n <N_RANKS> python3 run_pll.py <N_RUNS> <CPP_MODEL>
```

Where:

* `N_RANKS` – Number of CPU cores to use
* `N_RUNS` – Number of replicate simulations
* `CPP_MODEL` – Path to the C++ model file to compile and run

4. **Example command** (2 cores, 2 replicates, Gaussian cold-bias model):

```bash
mpiexec -n 2 python3 run_pll.py 2 ../run_gauss_tpc_cold_drate.cpp
```

---

## Output Files for all the models

Each simulation run generates multiple CSV files for downstream analysis. All the results are written out after a burn-in of 2000 time steps in the centre 10-patches. 

### Gauss Model

### 1. `result_<RUN_NO>_0.csv`

Time-series summary statistics **per patch and per time step**, including:

* `patch` – Patch index
* `time` – Time step
* `cons` – Population size
* `dispersal_events` – Number of dispersal events
* `T0_mut_count` – Mutations affecting optimal temperature
* `sigma_mut_count` – Mutations affecting curve width
* `average_dp` – Mean dispersal trait
* `q25_dp`, `q75_dp` – 25th and 75th percentiles of dispersal trait
* `average_T0` – Mean optimal temperature trait
* `q25_T0`, `q75_T0` – Quartiles of optimal temperature
* `average_sigma` – Mean curve width
* `q25_sigma`, `q75_sigma` – Quartiles of curve width
* `average_r0` – Mean maximum growth rate
* `q25_r0`, `q75_r0` – Quartiles of maximum growth rate
* `average_mut_eff` – Mean mutation effect size
* `q25_mut_eff`, `q75_mut_eff` – Quartiles of mutation effect size
* `sum_beta` – Intraspecific competition coefficient

---

### 2. `individual_trait_<TIME_STEP>_<RUN_NO>_0.csv`
Individual-level trait data at a given time step:

* individual 
* T0 - optimal temperature trait
* sigma - curve width
* disp_rate - dispersal trait
* patch_no - occupied patch
* mut_eff - if mutation occured, mutation effect


### Enzyme Catalysis Model

### 1. `result_<RUN_NO>_0.csv`

Time-series summary statistics **per patch and per time step**, including:

* `patch` – Patch index
* `time` – Time step
* `cons` – Population size
* `dispersal_events` – Number of dispersal events
* `DHm_mut_count` – Mutations affecting Activation enthalpy
* `DSm_mut_count` – Mutations affecting  Activation entropy
* `DCp_mut_count` – Mutations affecting Heat capacity change
* `average_dp` – Mean dispersal trait
* `q25_dp`, `q75_dp` – 25th and 75th percentiles of dispersal trait
* `average_DHm` – Mean Activation enthalpy
* `q25_DHm`, `q75_DHm` – Quartiles of Activation enthalpy
* `average_DSm` – Mean  Activation entropy
* `q25_DSm`, `q75_DSm` – Quartiles of  Activation entropy
* `average_DCp` – Mean Heat capacity change
* `q25_DCp`, `q75_DCp` – Quartiles of Heat capacity change
* `average_r0` – Mean maximum growth rate
* `q25_r0`, `q75_r0` – Quartiles of maximum growth rate
* `average_mut_eff` – Mean mutation effect size
* `q25_mut_eff`, `q75_mut_eff` – Quartiles of mutation effect size
* `sum_beta` – Intraspecific competition coefficient

---

### 2. `individual_trait_<TIME_STEP>_<RUN_NO>_0.csv`
Individual-level trait data at a given time step:

* individual 
* DHm - Activation enthalpy
* DSm - Activation entropy 
* DCp - Heat capacity change
* disp_rate - dispersal trait
* patch_no - occupied patch
* mut_eff - if mutation occured, mutation effect
---


### Protein Denaturation Model

### 1. `result_<RUN_NO>_0.csv`

Time-series summary statistics **per patch and per time step**, including:

* `patch` – Patch index
* `time` – Time step
* `cons` – Population size
* `dispersal_events` – Number of dispersal events
* `dGr_mut_count` – Mutations affecting Activation free energy
* `dSr_mut_count` – Mutations affecting Activation entropy
* `average_dp` – Mean dispersal trait
* `q25_dp`, `q75_dp` – 25th and 75th percentiles of dispersal trait
* `average_dGr` – Mean Activation free energy trait
* `q25_dGr`, `q75_dGr` – Quartiles of Activation free energy
* `average_dSr` – Mean Activation entropy
* `q25_dSr`, `q75_dSr` – Quartiles of Activation entropy
* `average_r0` – Mean maximum growth rate
* `q25_r0`, `q75_r0` – Quartiles of maximum growth rate
* `average_mut_eff` – Mean mutation effect size
* `q25_mut_eff`, `q75_mut_eff` – Quartiles of mutation effect size
* `sum_beta` – Intraspecific competition coefficient

---

### 2. `individual_trait_<TIME_STEP>_<RUN_NO>_0.csv`
Individual-level trait data at a given time step:

* individual 
* dGr - Activation free energy trait
* dSr - Activation entropy
* disp_rate - dispersal trait
* patch_no - occupied patch
* mut_eff - if mutation occured, mutation effect

## Citation

If you use this code or data in your research, please cite:

**Naik, S. & Fronhofer, E.** (2024). From proteins to species ranges: a framework for understanding thermal adaptation during range expansions. *Proceedings of the Royal Society B*. DOI: 10.1098/rspb


## Contact
Emanuel A. Fronhofer
Institut Des Sciences De L’evolution De Montpellier, France
Email: emanuel.fronhofer@umontpellier.fr

