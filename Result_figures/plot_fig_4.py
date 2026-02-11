import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
import seaborn as sns
import statistics
import os
import re
from alive_progress import alive_bar
import matplotlib.ticker as mticker

# Set global plotting parameters
fsize = 21
plt.rcParams.update({'font.size': 21})
plt.rcParams['lines.markersize'] = 6
plt.rcParams['lines.linewidth'] = 0.5

def ind_plot(filename, disp_stat, full_list, save_fol):
    """
    Plot individual mutation effects across the spatial gradient.
    
    Parameters:
    -----------
    filename : str
        Model directory name
    disp_stat : str
        Simulation type ('evo' for evolutionary simulations)
    full_list : list
        List of files in the simulation directory
    save_fol : str
        Folder to save supplementary figures
    """
    
    fig, ax = plt.subplots(1, 2)
    
    if filename == 'Enzyme_catalysis_model':
        # Process 4 individual replicates (replicates 4-7)
        for rep in range(4, 8, 1):
            # Filter files for this specific replicate
            filt_list_rep = [x for x in full_list if ("_" + str(rep)) in x]
            filt_list_rep = [x for x in filt_list_rep if ("individual") in x]
            
            # Process individual mutation data files with progress bar
            with alive_bar(len(filt_list_rep), bar='fish') as bar:
                for i in filt_list_rep:
                    # Read individual mutation data
                    data = pd.read_csv(filename + '/' + disp_stat + '/' + i, 
                                       skiprows=[0], delimiter=" ", header=None)
                    data.columns = data.iloc[0]
                    data = data[1:]
                    data = data.astype(float)
                    
                    # Filter data: exclude zero mutation effects and unrealistic trait values
                    data = data[(data["mut_eff"] != 0)]
                    data = data[(data["DCp"] < -10**(-8))]
                    data = data[(data["DHm"] < 0.5 - 10**(-8))]
                    data = data[(data["DSm"] < -10**(-8))]
                    
                    # Create color gradient based on patch number (spatial position)
                    colors = cm.coolwarm(np.linspace(0, 1, lan_size))
                    
                    # Plot mutation effects by spatial position
                    hist_ax[0].scatter(data["patch_no"], data["mut_eff"], 
                                       color=colors[data["patch_no"].astype(int)])
                    
                    # Plot mutation effects vs. ΔCp changes (excluding zero changes)
                    ax[0].scatter(data[data["DDCp"] != 0]["DDCp"],
                                 data[data["DDCp"] != 0]["mut_eff"],
                                 color=colors[data[data["DDCp"] != 0]["patch_no"].astype(int)])
                    
                    # Plot mutation effects vs. ΔHm changes (excluding zero changes)
                    ax[1].scatter(data[data["DDHm"] != 0]["DDHm"],
                                 data[data["DDHm"] != 0]["mut_eff"],
                                 color=colors[data[data["DDHm"] != 0]["patch_no"].astype(int)])
                    bar()
        
        # Label panels and set axis labels
        ax[0].text(0.95, 0.95, 'A', transform=ax[0].transAxes,
                   verticalalignment='top', horizontalalignment='right', fontweight='bold')
        ax[1].text(0.95, 0.95, 'B', transform=ax[1].transAxes,
                   verticalalignment='top', horizontalalignment='right', fontweight='bold')
        
        ax[1].set_xlabel("$dΔH^‡_{T{0}}$")  # Enthalpy change
        ax[1].set_ylabel("Mutation effect")
        ax[0].set_xlabel("$dΔC_{p}$")  # Heat capacity change
        ax[0].set_ylabel("Mutation effect")
        ax[1].set_ylim(-0.1, 0.1)
        
        # Configure histogram plot
        hist_ax[0].set_xlabel("Patch Number")
        hist_ax[0].set_ylabel("Fitness effect")
        fig.set_size_inches(20, 10)
        fig.tight_layout(pad=1)
    
    elif filename == 'Protein_denat_model':
        # Process protein denaturation model data
        for rep in range(4, 8, 1):
            filt_list_rep = [x for x in full_list if ("_" + str(rep)) in x]
            filt_list_rep = [x for x in filt_list_rep if ("individual") in x]
            
            with alive_bar(len(filt_list_rep), bar='fish') as bar:
                for i in filt_list_rep:
                    data = pd.read_csv(filename + '/' + disp_stat + '/' + i,
                                       skiprows=[0], delimiter=" ", header=None)
                    data.columns = data.iloc[0]
                    data = data[1:]
                    data = data.astype(float)
                    data = data[(data["mut_eff"] != 0)]
                    
                    colors = cm.coolwarm(np.linspace(0, 1, lan_size))
                    
                    # Plot mutation effects by spatial position
                    hist_ax[1].scatter(data["patch_no"], data["mut_eff"],
                                       color=colors[data["patch_no"].astype(int)])
                    
                    # Plot mutation effects vs. ΔGr changes (when ΔSr = 0)
                    ax[0].scatter(data.loc[(data["DdGr"] != 0) & (data["DdSr"] == 0)]["DdGr"],
                                 data.loc[(data["DdGr"] != 0) & (data["DdSr"] == 0)]["mut_eff"],
                                 color=colors[data[(data["DdGr"] != 0) & (data["DdSr"] == 0)]["patch_no"].astype(int)])
                    
                    # Plot mutation effects vs. ΔSr changes (when ΔGr = 0)
                    ax[1].scatter(data.loc[(data["DdGr"] == 0) & (data["DdSr"] != 0)]["DdSr"],
                                 data.loc[(data["DdGr"] == 0) & (data["DdSr"] != 0)]["mut_eff"],
                                 color=colors[data.loc[(data["DdGr"] == 0) & (data["DdSr"] != 0)]["patch_no"].astype(int)])
                    bar()
        
        # Label panels and set axis labels
        ax[0].text(0.95, 0.95, 'A', transform=ax[0].transAxes,
                   verticalalignment='top', horizontalalignment='right', fontweight='bold')
        ax[1].text(0.95, 0.95, 'B', transform=ax[1].transAxes,
                   verticalalignment='top', horizontalalignment='right', fontweight='bold')
        
        ax[1].set_xlabel("$dΔS^‡_{r}$")  # Entropy change
        ax[1].set_ylabel("Mutation effect")
        ax[1].set_ylim(-0.6, 0.6)
        ax[0].set_xlabel("$dΔG^‡_{r}$")  # Gibbs free energy change
        ax[0].set_ylabel("Mutation effect")
        ax[0].set_ylim(-0.6, 0.6)
        
        # Format x-axis ticks for better readability
        ax[1].xaxis.set_major_locator(mticker.MaxNLocator(nbins=3))
        hist_ax[1].set_xlabel("Patch Number")
        fig.set_size_inches(20, 10)
        fig.tight_layout(pad=1)
    
    elif filename == 'Gauss_model':
        # Process Gaussian model data
        for rep in range(4, 8, 1):
            filt_list_rep = [x for x in full_list if ("_" + str(rep)) in x]
            filt_list_rep = [x for x in filt_list_rep if ("individual") in x]
            
            with alive_bar(len(filt_list_rep), bar='fish') as bar:
                for i in filt_list_rep:
                    data = pd.read_csv(filename + '/' + disp_stat + '/' + i,
                                       skiprows=[0], delimiter=" ", header=None)
                    data.columns = data.iloc[0]
                    data = data[1:]
                    data = data.astype(float)
                    data = data[(data["mut_eff"] != 0)]
                    
                    colors = cm.coolwarm(np.linspace(0, 1, lan_size))
                    
                    # Plot mutation effects vs. T0 (optimum temperature)
                    ax[0].scatter(data["T0"], data["mut_eff"],
                                 color=colors[data["patch_no"].astype(int)])
                    
                    # Plot mutation effects vs. σ (niche width)
                    ax[1].scatter(data["sigma"], data["mut_eff"],
                                 color=colors[data["patch_no"].astype(int)])
                    
                    # Plot mutation effects by spatial position
                    hist_ax[2].scatter(data["patch_no"], data["mut_eff"],
                                       color=colors[data["patch_no"].astype(int)])
                    bar()
        
        # Label panels and set axis labels
        ax[0].text(0.95, 0.95, 'A', transform=ax[0].transAxes,
                   verticalalignment='top', horizontalalignment='right', fontweight='bold')
        ax[1].text(0.95, 0.95, 'B', transform=ax[1].transAxes,
                   verticalalignment='top', horizontalalignment='right', fontweight='bold')
        
        ax[1].set_xlabel(r"$\sigma$")  # Niche width
        ax[1].set_ylabel("Mutation effect")
        ax[0].set_xlabel(r"$T_0$")  # Optimum temperature
        ax[0].set_ylabel("Mutation effect")
        hist_ax[2].set_xlabel("Patch Number")
        fig.set_size_inches(20, 10)
        fig.tight_layout(pad=1)
    
    # Save the histogram figure
    hist.tight_layout(pad=0.75)
    hist.savefig(save_fol+'/Fig_4_mut_eff_spatial_dis.png', dpi=500)

# ============================================================================
# SETUP HISTOGRAM FIGURE (3 panels for the 3 models)
# ============================================================================
hist, hist_ax = plt.subplots(1, 3)

# Label panels A, B, C
hist_ax[0].text(0.95, 0.95, 'A', transform=hist_ax[0].transAxes,
                verticalalignment='top', horizontalalignment='right', fontweight='bold')
hist_ax[1].text(0.95, 0.95, 'B', transform=hist_ax[1].transAxes,
                verticalalignment='top', horizontalalignment='right', fontweight='bold')
hist_ax[2].text(0.95, 0.95, 'C', transform=hist_ax[2].transAxes,
                verticalalignment='top', horizontalalignment='right', fontweight='bold')

# Set titles for each panel
hist_ax[0].set_title('Enzyme catalysis based model', loc='center')
hist_ax[1].set_title('Protein denaturation based model', loc='center')
hist_ax[2].set_title('Gaussian model', loc='center')

# Set consistent axis limits for comparison
hist_ax[0].set_ylim([-0.8, 0.6])
hist_ax[1].set_ylim([-0.8, 0.6])
hist_ax[2].set_ylim([-0.8, 0.6])

# Hide y-axis labels for middle and right panels for cleaner layout
hist_ax[1].set_yticklabels([])
hist_ax[2].set_yticklabels([])

# Set x-axis limits to match landscape size
hist_ax[0].set_xlim([0, 144])
hist_ax[1].set_xlim([0, 144])
hist_ax[2].set_xlim([0, 144])

# Adjust figure size
hist.set_size_inches(20, 15)

# Landscape size (number of patches in the spatial gradient)
lan_size = 144

# Folder for saving supplementary figures
save_fol = 'Paper_figs'

# ============================================================================
# PROCESS EACH MODEL
# ============================================================================

# Protein denaturation model
filename = 'Protein_denat_model'
disp_stat = "evo"
full_list = os.listdir(filename + '/' + disp_stat + '/')
ind_plot(filename, disp_stat, full_list, save_fol)
print('Protein denaturation model done')

# Enzyme catalysis model
filename = 'Enzyme_catalysis_model'
disp_stat = 'evo'
full_list = os.listdir(filename + '/' + disp_stat + '/')
ind_plot(filename, disp_stat, full_list, save_fol)
print('Enzyme catalysis model done')

# Gaussian model
filename = 'Gauss_model'
disp_stat = 'evo'
full_list = os.listdir(filename + '/' + disp_stat + '/')
ind_plot(filename, disp_stat, full_list, save_fol)
print('Gauss model done')
