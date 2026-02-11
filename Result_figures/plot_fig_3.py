import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
import seaborn as sns
import statistics
import os
import re
import matplotlib.ticker as ticker

# Set global plotting parameters
fsize = 28
plt.rcParams.update({'font.size': 28})
plt.rcParams['lines.markersize'] = 6
plt.rcParams['lines.linewidth'] = 2.5

# Define a function to extract numbers and convert them to float
def extract_floats(x):
    """Extract all floating-point numbers from a string and convert to float list."""
    numbers = re.findall(r"[-+]?(?:\d*\.*\d+e[-+]*\d+)|[-+]?(?:\d*\.*\d+)", x)
    return list(map(float, numbers)) if numbers else []

# List of columns in niche_pf that contain string representations of lists
columns_to_convert = [
    "hot_niche", "cold_niche", "core_niche",
    "q_25_hot_br_niche", "q_25_cold_br_niche", "q_25_core_br_niche",
    "q_75_hot_br_niche", "q_75_cold_br_niche", "q_75_core_br_niche"
]

def trait_pf_plot(filename, disp_stat, full_list, time_steps, temp_range, 
                  axs_grw_vr, axs_niche, axs_trt1, axs_trt2, colors, model_label):
    """
    Plot trait evolution data for a given thermal adaptation model.
    
    Parameters:
    -----------
    filename : str
        Model directory name (e.g., 'Enzyme_catalysis_model')
    disp_stat : str
        Simulation type ('eco' or 'evo')
    full_list : list
        List of files in the simulation directory
    time_steps : range
        Time steps to analyze
    temp_range : numpy array
        Temperature range for niche plots
    axs_grw_vr : matplotlib.axes.Axes (2-element list)
        Axes for growth rate plots (between-replicate and within-replicate variation)
    axs_niche : matplotlib.axes.Axes (3-element list)
        Axes for niche plots (one per model)
    axs_trt1 : matplotlib.axes.Axes (2-element list)
        Axes for trait 1 evolution (e.g., DHm, dGr, or sigma)
    axs_trt2 : matplotlib.axes.Axes (2-element list)
        Axes for trait 2 evolution (e.g., DCp, dSr, or T0)
    colors : numpy array
        Color code for the model
    model_label : numpy array
        Label for the model (for legends)
    """
    
    ind = 0  # Index counter (unused in current code)
    para_file = filename + '/' + disp_stat + '/' + '/parameters_eavar.in'
    
    # Load processed fitness data files
    growth_pf = pd.read_csv(filename + '/' + disp_stat + '/' + 'growth_pf.csv', sep=" ")
    mut_eff_pf = pd.read_csv(filename + '/' + disp_stat + '/' + 'mut_eff_pf.csv', sep=" ")
    niche_pf = pd.read_csv(filename + '/' + disp_stat + '/' + 'niche_pf.csv', sep=" ")
    
    # Convert string representations of lists to actual lists of floats
    niche_pf[columns_to_convert] = niche_pf[columns_to_convert].map(extract_floats)
    
    # ========================================================================
    # MUTATION EFFECT PLOTS (Between-replicate and within-replicate variation)
    # ========================================================================
    axs_mut_eff[0].plot(mut_eff_pf["time_step"], mut_eff_pf["hot_mut_eff"], 
                        color=colors[0], label=model_label[0])
    axs_mut_eff[0].fill_between(mut_eff_pf["time_step"],
                                mut_eff_pf["q_25_hot_br_mut_eff"],
                                mut_eff_pf["q_75_hot_br_mut_eff"],
                                color=colors[0], alpha=0.4)
    axs_mut_eff[1].plot(mut_eff_pf["time_step"], mut_eff_pf["hot_wr_mut_eff"],
                        color=colors[0])
    axs_mut_eff[1].fill_between(mut_eff_pf["time_step"],
                                mut_eff_pf["q_25_hot_wr_mut_eff"],
                                mut_eff_pf["q_75_hot_wr_mut_eff"],
                                alpha=0.4, color=colors[0])
    
    # Cold front mutation effects (with cross-hatch pattern)
    axs_mut_eff[0].plot(mut_eff_pf["time_step"], mut_eff_pf["cold_mut_eff"],
                        color=colors[0], label=model_label[0])
    axs_mut_eff[0].fill_between(mut_eff_pf["time_step"],
                                mut_eff_pf["q_25_cold_br_mut_eff"],
                                mut_eff_pf["q_75_cold_br_mut_eff"],
                                hatch="X", edgecolor='white', facecolor=colors[0], alpha=0.4)
    axs_mut_eff[1].plot(mut_eff_pf["time_step"], mut_eff_pf["cold_wr_mut_eff"],
                        color=colors[0])
    axs_mut_eff[1].fill_between(mut_eff_pf["time_step"],
                                mut_eff_pf["q_25_cold_wr_mut_eff"],
                                mut_eff_pf["q_75_cold_wr_mut_eff"],
                                hatch="X", edgecolor='white', alpha=0.4, facecolor=colors[0])
    
    axs_mut_eff[0].set_xlabel('Time step')
    axs_mut_eff[0].set_ylabel("Beneficial mutation effect")
    axs_mut_eff[1].set_xlabel('Time step')
    axs_mut_eff[1].set_ylabel("Beneficial mutation effect" + ' within replicate variation')
    
    # ========================================================================
    # GROWTH RATE PLOTS
    # ========================================================================
    axs_grw_vr[0].plot(growth_pf["time_step"], growth_pf["hot_growth"],
                       color=colors[0], label=model_label[0])
    axs_grw_vr[0].fill_between(growth_pf["time_step"],
                               growth_pf["q_25_hot_br_growth"],
                               growth_pf["q_75_hot_br_growth"],
                               color=colors[0], alpha=0.4)
    axs_grw_vr[1].plot(growth_pf["time_step"], growth_pf["hot_wr_growth"],
                       color=colors[0])
    axs_grw_vr[1].fill_between(growth_pf["time_step"],
                               growth_pf["q_25_hot_wr_growth"],
                               growth_pf["q_75_hot_wr_growth"],
                               alpha=0.4, color=colors[0])
    
    axs_grw_vr[0].plot(growth_pf["time_step"], growth_pf["cold_growth"],
                       color=colors[0])
    axs_grw_vr[0].fill_between(growth_pf["time_step"],
                               growth_pf["q_25_cold_br_growth"],
                               growth_pf["q_75_cold_br_growth"],
                               color=colors[0], hatch="X", edgecolor='white', alpha=0.4)
    axs_grw_vr[1].plot(growth_pf["time_step"], growth_pf["cold_wr_growth"],
                       color=colors[0], linestyle='dotted')
    axs_grw_vr[1].fill_between(growth_pf["time_step"],
                               growth_pf["q_25_cold_wr_growth"],
                               growth_pf["q_75_cold_wr_growth"],
                               alpha=0.4, color=colors[0], hatch="X", edgecolor='white')
    
    axs_grw_vr[0].set_xlabel('Time step')
    axs_grw_vr[0].set_ylabel('Population growth rate')
    axs_grw_vr[1].set_xlabel('Time step')
    axs_grw_vr[1].set_ylabel('Population variation in growth rate')
    
    # ========================================================================
    # OPTIMUM TEMPERATURE PLOTS
    # ========================================================================
    axs_topt[0].plot(niche_pf["time_step"], niche_pf["hot_Topt"],
                     color=colors[0], label=model_label[0])
    axs_topt[0].plot(niche_pf["time_step"], niche_pf["cold_Topt"],
                     color=colors[0], label=model_label[0])
    axs_topt[0].fill_between(niche_pf["time_step"],
                             niche_pf["q_25_hot_br_Topt"],
                             niche_pf["q_75_hot_br_Topt"],
                             color=colors[0], label=model_label[0], alpha=0.4)
    axs_topt[0].fill_between(niche_pf["time_step"],
                             niche_pf["q_25_cold_br_Topt"],
                             niche_pf["q_75_cold_br_Topt"],
                             color=colors[0], hatch="X", edgecolor='white', alpha=0.4)
    
    # Plot environmental temperatures (hot and cold patch temperatures)
    axs_topt[0].plot(niche_pf["time_step"], 273 + 0.5 * niche_pf["hot_patch_no"],
                     color='r', label=model_label[0])
    axs_topt[0].plot(niche_pf["time_step"], 273 + 0.5 * niche_pf["cold_patch_no"],
                     color='b', label=model_label[0])
    
    axs_topt[1].plot(niche_pf["time_step"], niche_pf["cold_wr_Topt"],
                     color=colors[0], linestyle='dotted')
    axs_topt[1].fill_between(niche_pf["time_step"],
                             niche_pf['q_25_cold_wr_Topt'],
                             niche_pf["q_75_cold_wr_Topt"],
                             alpha=0.4, color=colors[0], hatch="X", edgecolor='white')
    axs_topt[1].plot(niche_pf["time_step"], niche_pf["hot_wr_Topt"],
                     color=colors[0], linestyle='dotted')
    axs_topt[1].fill_between(niche_pf["time_step"],
                             niche_pf["q_25_hot_wr_Topt"],
                             niche_pf["q_75_hot_wr_Topt"],
                             alpha=0.4, color=colors[0], hatch="X", edgecolor='white')
    
    axs_topt[0].set_xlabel('Time step')
    axs_topt[0].set_ylabel('Optimum temperature')
    axs_topt[1].set_xlabel('Time step')
    axs_topt[1].set_ylabel('Population variation in optimum temperature')
    
    # ========================================================================
    # MODEL-SPECIFIC TRAIT PLOTS
    # ========================================================================
    
    if filename == 'Enzyme_catalysis_model':
        # Load enzyme catalysis model specific trait data
        DHm_pf = pd.read_csv(filename + '/' + disp_stat + '/' + 'DHm_pf.csv', sep=" ")
        DSm_pf = pd.read_csv(filename + '/' + disp_stat + '/' + 'DSm_pf.csv', sep=" ")
        DCp_pf = pd.read_csv(filename + '/' + disp_stat + '/' + 'DCp_pf.csv', sep=" ")
        
        # Plot enthalpy change (ΔH‡) evolution
        axs_trt1[0].plot(DHm_pf["time_step"], DHm_pf["hot_DHm"],
                         color=colors[0], label=model_label[0])
        axs_trt1[0].fill_between(DHm_pf["time_step"],
                                 DHm_pf["q_25_hot_br_DHm"],
                                 DHm_pf["q_75_hot_br_DHm"],
                                 color=colors[0], alpha=0.4)
        axs_trt1[1].plot(DHm_pf["time_step"], DHm_pf["hot_wr_DHm"],
                         color=colors[0])
        axs_trt1[1].fill_between(DHm_pf["time_step"],
                                 DHm_pf["q_25_hot_wr_DHm"],
                                 DHm_pf["q_75_hot_wr_DHm"],
                                 alpha=0.4, color=colors[0])
        
        axs_trt1[0].plot(DHm_pf["time_step"], DHm_pf["cold_DHm"],
                         color=colors[0], label=model_label[0])
        axs_trt1[0].fill_between(DHm_pf["time_step"],
                                 DHm_pf["q_25_cold_br_DHm"],
                                 DHm_pf["q_75_cold_br_DHm"],
                                 hatch="X", edgecolor='white', facecolor=colors[0], alpha=0.4)
        axs_trt1[1].plot(DHm_pf["time_step"], DHm_pf["cold_wr_DHm"],
                         color=colors[0])
        axs_trt1[1].fill_between(DHm_pf["time_step"],
                                 DHm_pf["q_25_cold_wr_DHm"],
                                 DHm_pf["q_75_cold_wr_DHm"],
                                 hatch="X", edgecolor='white', alpha=0.4, facecolor=colors[0])
        
        axs_trt1[0].set_xlabel('Time step')
        axs_trt1[0].set_ylabel("$dΔH^‡_{T{0}}$")  # Enthalpy change
        axs_trt1[1].set_xlabel('Time step')
        axs_trt1[1].set_ylabel("$dΔH^‡_{T{0}}$" + ' within replicate variation')
        
        # Plot heat capacity change (ΔCp) evolution
        axs_trt2[0].plot(DCp_pf["time_step"], DCp_pf["hot_DCp"],
                         color=colors[0], label=model_label[0])
        axs_trt2[0].fill_between(DCp_pf["time_step"],
                                 DCp_pf["q_25_hot_br_DCp"],
                                 DCp_pf["q_75_hot_br_DCp"],
                                 color=colors[0], alpha=0.4)
        axs_trt2[1].plot(DCp_pf["time_step"], DCp_pf["hot_wr_DCp"],
                         color=colors[0])
        axs_trt2[1].fill_between(DCp_pf["time_step"],
                                 DCp_pf["q_25_hot_wr_DCp"],
                                 DCp_pf["q_75_hot_wr_DCp"],
                                 alpha=0.4, color=colors[0])
        
        axs_trt2[0].plot(DCp_pf["time_step"], DCp_pf["cold_DCp"],
                         color=colors[0], label=model_label[0])
        axs_trt2[0].fill_between(DCp_pf["time_step"],
                                 DCp_pf["q_25_cold_br_DCp"],
                                 DCp_pf["q_75_cold_br_DCp"],
                                 hatch="X", edgecolor='white', facecolor=colors[0], alpha=0.4)
        axs_trt2[1].plot(DCp_pf["time_step"], DCp_pf["cold_wr_DCp"],
                         color=colors[0])
        axs_trt2[1].fill_between(DCp_pf["time_step"],
                                 DCp_pf["q_25_cold_wr_DCp"],
                                 DCp_pf["q_75_cold_wr_DCp"],
                                 hatch="X", edgecolor='white', alpha=0.4, facecolor=colors[0])
        
        axs_trt2[0].set_xlabel('Time step')
        axs_trt2[0].set_ylabel("$dΔC_{p}$")  # Heat capacity change
        axs_trt2[1].set_xlabel('Time step', fontsize=fsize)
        axs_trt2[1].set_ylabel("$dΔC_{p}$" + ' within replicate variation', fontsize=fsize)
        
        # Plot thermal niches at equilibrium (final time step)
        axs_niche[0].plot(temp_range, niche_pf["hot_niche"].iloc[-1],
                          color='r', label="Hot front")
        axs_niche[0].fill_between(temp_range,
                                  niche_pf["q_25_hot_br_niche"].iloc[-1],
                                  niche_pf["q_75_hot_br_niche"].iloc[-1],
                                  color='r', alpha=0.4)
        axs_niche[0].plot(temp_range, niche_pf["cold_niche"].iloc[-1],
                          color='b', label='Cold front')
        axs_niche[0].fill_between(temp_range,
                                  niche_pf["q_25_cold_br_niche"].iloc[-1],
                                  niche_pf["q_75_cold_br_niche"].iloc[-1],
                                  color='b', alpha=0.4)
        axs_niche[0].set_ylim([0.0, 5.1])
        axs_niche[0].set_xlim([273, 345])
        axs_niche[0].set_xlabel('Temperature (K)')
        axs_niche[0].set_ylabel('Population growth rate')
        axs_niche[0].tick_params(axis='x', pad=15)
    
    elif filename == 'Protein_denat_model':
        # Load protein denaturation model specific trait data
        dGr_pf = pd.read_csv(filename + '/' + disp_stat + '/' + 'dGr_pf.csv', sep=" ")
        dSr_pf = pd.read_csv(filename + '/' + disp_stat + '/' + 'dSr_pf.csv', sep=" ")
        
        # Plot Gibbs free energy change (ΔG‡) evolution
        axs_trt1[0].plot(dGr_pf["time_step"], dGr_pf["hot_dGr"],
                         color=colors[0], label=model_label[0])
        axs_trt1[0].fill_between(dGr_pf["time_step"],
                                 dGr_pf["q_25_hot_br_dGr"],
                                 dGr_pf["q_75_hot_br_dGr"],
                                 color=colors[0], alpha=0.4)
        axs_trt1[1].plot(dGr_pf["time_step"], dGr_pf["hot_wr_dGr"],
                         color=colors[0])
        axs_trt1[1].fill_between(dGr_pf["time_step"],
                                 dGr_pf["q_25_hot_wr_dGr"],
                                 dGr_pf["q_75_hot_wr_dGr"],
                                 alpha=0.4, color=colors[0])
        
        axs_trt1[0].plot(dGr_pf["time_step"], dGr_pf["cold_dGr"],
                         color=colors[0])
        axs_trt1[0].fill_between(dGr_pf["time_step"],
                                 dGr_pf["q_25_cold_br_dGr"],
                                 dGr_pf["q_75_cold_br_dGr"],
                                 hatch="X", edgecolor='white', facecolor=colors[0], alpha=0.4)
        axs_trt1[1].plot(dGr_pf["time_step"], dGr_pf["cold_wr_dGr"],
                         color=colors[0])
        axs_trt1[1].fill_between(dGr_pf["time_step"],
                                 dGr_pf["q_25_cold_wr_dGr"],
                                 dGr_pf["q_75_cold_wr_dGr"],
                                 hatch="X", edgecolor='white', alpha=0.4, facecolor=colors[0])
        
        axs_trt1[0].set_xlabel('Time step', fontsize=fsize)
        axs_trt1[0].set_ylabel("$dΔG^‡_{r}$", fontsize=fsize)  # Gibbs free energy
        axs_trt1[1].set_xlabel('Time step', fontsize=fsize)
        axs_trt1[1].set_ylabel("$dΔG^‡_{r}$" + ' within replicate variation', fontsize=fsize)
        
        # Plot entropy change (ΔS‡) evolution
        axs_trt2[0].plot(dSr_pf["time_step"], dSr_pf["hot_dSr"],
                         color=colors[0])
        axs_trt2[0].fill_between(dSr_pf["time_step"],
                                 dSr_pf["q_25_hot_br_dSr"],
                                 dSr_pf["q_75_hot_br_dSr"],
                                 color=colors[0], alpha=0.4)
        axs_trt2[1].plot(dSr_pf["time_step"], dSr_pf["hot_wr_dSr"],
                         color=colors[0])
        axs_trt2[1].fill_between(dSr_pf["time_step"],
                                 dSr_pf["q_25_hot_wr_dSr"],
                                 dSr_pf["q_75_hot_wr_dSr"],
                                 alpha=0.4, color=colors[0])
        
        axs_trt2[0].plot(dSr_pf["time_step"], dSr_pf["cold_dSr"],
                         color=colors[0], label=model_label[0])
        axs_trt2[0].fill_between(dSr_pf["time_step"],
                                 dSr_pf["q_25_cold_br_dSr"],
                                 dSr_pf["q_75_cold_br_dSr"],
                                 hatch="X", edgecolor='white', facecolor=colors[0], alpha=0.4)
        axs_trt2[1].plot(dSr_pf["time_step"], dSr_pf["cold_wr_dSr"],
                         color=colors[0])
        axs_trt2[1].fill_between(dSr_pf["time_step"],
                                 dSr_pf["q_25_cold_wr_dSr"],
                                 dSr_pf["q_75_cold_wr_dSr"],
                                 hatch="X", edgecolor='white', alpha=0.4, facecolor=colors[0])
        
        axs_trt2[0].set_xlabel('Time step', fontsize=fsize)
        axs_trt2[0].set_ylabel("$dΔS^‡_{r}$", fontsize=fsize)  # Entropy change
        axs_trt2[1].set_xlabel('Time step', fontsize=fsize)
        axs_trt2[1].set_ylabel("$dΔS^‡_{r}$" + ' within replicate variation', fontsize=fsize)
        
        # Plot thermal niches at equilibrium
        axs_niche[1].plot(temp_range, niche_pf["hot_niche"].iloc[-1], color='r')
        axs_niche[1].fill_between(temp_range,
                                  niche_pf["q_25_hot_br_niche"].iloc[-1],
                                  niche_pf["q_75_hot_br_niche"].iloc[-1],
                                  color='r', alpha=0.4)
        axs_niche[1].plot(temp_range, niche_pf["cold_niche"].iloc[-1], color='b')
        axs_niche[1].fill_between(temp_range,
                                  niche_pf["q_25_cold_br_niche"].iloc[-1],
                                  niche_pf["q_75_cold_br_niche"].iloc[-1],
                                  color='b', alpha=0.4)
        axs_niche[1].set_ylim([0.0, 0.8])
        axs_niche[1].set_xlim([273, 345])
        axs_niche[1].set_xlabel('Temperature (K)')
        axs_niche[1].tick_params(axis='x', pad=15)
    
    elif filename == 'Gauss_model':
        # Load Gaussian model specific trait data
        sigma_pf = pd.read_csv(filename + '/' + disp_stat + '/' + 'sigma_pf.csv', sep=" ")
        T0_pf = pd.read_csv(filename + '/' + disp_stat + '/' + 'T0_pf.csv', sep=" ")
        
        # Plot niche width (σ) evolution
        axs_trt1[0].plot(sigma_pf["time_step"], sigma_pf["hot_sigma"],
                         color=colors[0], label=model_label[0])
        axs_trt1[0].fill_between(sigma_pf["time_step"],
                                 sigma_pf["q_25_hot_br_sigma"],
                                 sigma_pf["q_75_hot_br_sigma"],
                                 color=colors[0], alpha=0.4)
        axs_trt1[1].plot(sigma_pf["time_step"], sigma_pf["hot_wr_sigma"],
                         color=colors[0])
        axs_trt1[1].fill_between(sigma_pf["time_step"],
                                 sigma_pf["q_25_hot_wr_sigma"],
                                 sigma_pf["q_75_hot_wr_sigma"],
                                 alpha=0.4, color=colors[0])
        
        axs_trt1[0].plot(sigma_pf["time_step"], sigma_pf["cold_sigma"],
                         color=colors[0])
        axs_trt1[0].fill_between(sigma_pf["time_step"],
                                 sigma_pf["q_25_cold_br_sigma"],
                                 sigma_pf["q_75_cold_br_sigma"],
                                 hatch="X", edgecolor='white', facecolor=colors[0], alpha=0.4)
        axs_trt1[1].plot(sigma_pf["time_step"], sigma_pf["cold_wr_sigma"],
                         color=colors[0])
        axs_trt1[1].fill_between(sigma_pf["time_step"],
                                 sigma_pf["q_25_cold_wr_sigma"],
                                 sigma_pf["q_75_cold_wr_sigma"],
                                 hatch="X", edgecolor='white', alpha=0.4, facecolor=colors[0])
        
        axs_trt1[0].set_xlabel('Time step', fontsize=fsize)
        axs_trt1[0].set_ylabel(r"$\sigma$", fontsize=fsize)  # Niche width
        axs_trt1[1].set_xlabel('Time step', fontsize=fsize)
        axs_trt1[1].set_ylabel(r"$\sigma$" + ' within replicate variation', fontsize=fsize)
        
        # Plot optimum temperature (T₀) evolution
        axs_trt2[0].plot(T0_pf["time_step"], T0_pf["hot_T0"],
                         color=colors[0], label=model_label[0])
        axs_trt2[0].fill_between(T0_pf["time_step"],
                                 T0_pf["q_25_hot_br_T0"],
                                 T0_pf["q_75_hot_br_T0"],
                                 color=colors[0], alpha=0.4)
        axs_trt2[1].plot(T0_pf["time_step"], T0_pf["hot_wr_T0"],
                         color=colors[0])
        axs_trt2[1].fill_between(T0_pf["time_step"],
                                 T0_pf["q_25_hot_wr_T0"],
                                 T0_pf["q_75_hot_wr_T0"],
                                 alpha=0.4, color=colors[0])
        
        axs_trt2[0].plot(T0_pf["time_step"], T0_pf["cold_T0"],
                         color=colors[0])
        axs_trt2[0].fill_between(T0_pf["time_step"],
                                 T0_pf["q_25_cold_br_T0"],
                                 T0_pf["q_75_cold_br_T0"],
                                 hatch="X", edgecolor='white', facecolor=colors[0], alpha=0.4)
        axs_trt2[1].plot(T0_pf["time_step"], T0_pf["cold_wr_T0"],
                         color=colors[0])
        axs_trt2[1].fill_between(T0_pf["time_step"],
                                 T0_pf["q_25_cold_wr_T0"],
                                 T0_pf["q_75_cold_wr_T0"],
                                 hatch="X", edgecolor='white', alpha=0.4, facecolor=colors[0])
        
        axs_trt2[0].set_xlabel('Time step', fontsize=fsize)
        axs_trt2[0].set_ylabel(r"$T_0$", fontsize=fsize)  # Optimum temperature
        axs_trt2[1].set_xlabel('Time step', fontsize=fsize)
        axs_trt2[1].set_ylabel(r"$T_0$" + ' within replicate variation', fontsize=fsize)
        
        # Plot thermal niches at equilibrium
        axs_niche[2].plot(temp_range, niche_pf["hot_niche"].iloc[-1], color='r')
        axs_niche[2].fill_between(temp_range,
                                  niche_pf["q_25_hot_br_niche"].iloc[-1],
                                  niche_pf["q_75_hot_br_niche"].iloc[-1],
                                  color='r', alpha=0.4)
        axs_niche[2].plot(temp_range, niche_pf["cold_niche"].iloc[-1], color='b')
        axs_niche[2].fill_between(temp_range,
                                  niche_pf["q_25_cold_br_niche"].iloc[-1],
                                  niche_pf["q_75_cold_br_niche"].iloc[-1],
                                  color='b', alpha=0.4)
        axs_niche[2].set_ylim([0.0, 0.8])
        axs_niche[2].set_xlim([273, 345])
        axs_niche[2].set_xlabel('Temperature (K)')
        axs_niche[2].tick_params(axis='x', pad=15)

# ============================================================================
# FIGURE SETUP - GROWTH RATE PLOTS
# ============================================================================
fig_grw_vr, axs_grw_vr = plt.subplots(1, 2)
axs_grw_vr[0].text(0.95, 0.95, 'A', transform=axs_grw_vr[0].transAxes,
                   verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_grw_vr[1].text(0.95, 0.95, 'B', transform=axs_grw_vr[1].transAxes,
                   verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

# ============================================================================
# FIGURE SETUP - NICHE PLOTS (3 panels, one for each model)
# ============================================================================
fig_niche, axs_niche = plt.subplots(1, 3)
axs_niche[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
titles = ['A', 'B', 'C']

for j in range(3):
    axs_niche[j].text(0.95, 0.95, titles[j], transform=axs_niche[j].transAxes,
                      verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

# Adjust layout and set titles
plt.subplots_adjust(wspace=0.5, hspace=0.3)
axs_niche[0].set_title('Enzyme catalysis based model', loc='center', fontsize=28)
axs_niche[1].set_title('Protein denaturation based model', loc='center', fontsize=28)
axs_niche[2].set_title('Gaussian model', loc='center', fontsize=28)

# ============================================================================
# FIGURE SETUP - TRAIT PLOTS (Recreated for each model)
# ============================================================================

# Enzyme catalysis model plots
fig_trt1, axs_trt1 = plt.subplots(1, 2)
axs_trt1[0].text(0.95, 0.95, 'A', transform=axs_trt1[0].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_trt1[1].text(0.95, 0.95, 'B', transform=axs_trt1[1].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

fig_trt2, axs_trt2 = plt.subplots(1, 2)
axs_trt2[0].text(0.95, 0.95, 'A', transform=axs_trt2[0].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_trt2[1].text(0.95, 0.95, 'B', transform=axs_trt2[1].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

fig_mut_eff, axs_mut_eff = plt.subplots(1, 2)
axs_mut_eff[0].text(0.95, 0.95, 'A', transform=axs_mut_eff[0].transAxes,
                    verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_mut_eff[1].text(0.95, 0.95, 'B', transform=axs_mut_eff[1].transAxes,
                    verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

fig_topt, axs_topt = plt.subplots(1, 2)
axs_topt[0].text(0.95, 0.95, 'A', transform=axs_topt[0].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_topt[1].text(0.95, 0.95, 'B', transform=axs_topt[1].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

# ============================================================================
# ENZYME CATALYSIS MODEL PLOTTING
# ============================================================================
filename = 'Enzyme_catalysis_model'
disp_stat = 'evo'
full_list = os.listdir(filename + '/' + disp_stat + '/')
time_steps = range(2000, 2800, 2)
temp_range = np.arange(273, 400, 0.5)
colors = np.array(['coral'])
model_label = np.array(['Enzyme catalysis based model'])
trait_pf_plot(filename, disp_stat, full_list, time_steps, temp_range,
              axs_grw_vr, axs_niche, axs_trt1, axs_trt2, colors, model_label)

# ============================================================================
# PROTEIN DENATURATION MODEL PLOTTING
# ============================================================================
# Recreate figures for protein denaturation model
fig_trt1, axs_trt1 = plt.subplots(1, 2)
axs_trt1[0].text(0.95, 0.95, 'A', transform=axs_trt1[0].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_trt1[1].text(0.95, 0.95, 'B', transform=axs_trt1[1].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

fig_trt2, axs_trt2 = plt.subplots(1, 2)
axs_trt2[0].text(0.95, 0.95, 'A', transform=axs_trt2[0].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_trt2[1].text(0.95, 0.95, 'B', transform=axs_trt2[1].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

fig_mut_eff, axs_mut_eff = plt.subplots(1, 2)
axs_mut_eff[0].text(0.95, 0.95, 'A', transform=axs_mut_eff[0].transAxes,
                    verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_mut_eff[1].text(0.95, 0.95, 'B', transform=axs_mut_eff[1].transAxes,
                    verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

fig_topt, axs_topt = plt.subplots(1, 2)
axs_topt[0].text(0.95, 0.95, 'A', transform=axs_topt[0].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_topt[1].text(0.95, 0.95, 'B', transform=axs_topt[1].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

filename = 'Protein_denat_model'
disp_stat = 'evo'
full_list = os.listdir(filename + '/' + disp_stat + '/')
time_steps = range(2000, 2800, 2)
temp_range = np.arange(273, 345, 0.5)
colors = np.array(['goldenrod'])
model_label = np.array(['Protein denaturation based model'])
trait_pf_plot(filename, disp_stat, full_list, time_steps, temp_range,
              axs_grw_vr, axs_niche, axs_trt1, axs_trt2, colors, model_label)

# ============================================================================
# GAUSSIAN MODEL PLOTTING
# ============================================================================
# Recreate figures for Gaussian model
fig_trt1, axs_trt1 = plt.subplots(1, 2)
axs_trt1[0].text(0.95, 0.95, 'A', transform=axs_trt1[0].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_trt1[1].text(0.95, 0.95, 'B', transform=axs_trt1[1].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

fig_trt2, axs_trt2 = plt.subplots(1, 2)
axs_trt2[0].text(0.95, 0.95, 'A', transform=axs_trt2[0].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_trt2[1].text(0.95, 0.95, 'B', transform=axs_trt2[1].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

fig_mut_eff, axs_mut_eff = plt.subplots(1, 2)
axs_mut_eff[0].text(0.95, 0.95, 'A', transform=axs_mut_eff[0].transAxes,
                    verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_mut_eff[1].text(0.95, 0.95, 'B', transform=axs_mut_eff[1].transAxes,
                    verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

fig_topt, axs_topt = plt.subplots(1, 2)
axs_topt[0].text(0.95, 0.95, 'A', transform=axs_topt[0].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)
axs_topt[1].text(0.95, 0.95, 'B', transform=axs_topt[1].transAxes,
                 verticalalignment='top', horizontalalignment='right', fontweight='bold', fontsize=30)

filename = 'Gauss_model'
disp_stat = 'evo'
full_list = os.listdir(filename + '/' + disp_stat + '/')
time_steps = range(2000, 2800, 2)
temp_range = np.arange(273, 345, 0.5)
colors = np.array(['silver'])
model_label = np.array(['Gaussian model'])
trait_pf_plot(filename, disp_stat, full_list, time_steps, temp_range,
              axs_grw_vr, axs_niche, axs_trt1, axs_trt2, colors, model_label)

# ============================================================================
# FINAL FORMATTING AND SAVING
# ============================================================================
# Add legends to appropriate plots
axs_grw_vr[0].legend(loc=2, frameon=False)
axs_niche[0].legend(loc=2, frameon=False)

# Adjust figure size and save niche plot
figure = fig_niche
figure.set_size_inches(30, 10)
fig_niche.tight_layout(pad=0.75)
fig_niche.savefig('Paper_figs/Fig_3_pf_niche.png', dpi=500)
