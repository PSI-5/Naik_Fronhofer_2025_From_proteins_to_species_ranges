import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Set global plotting parameters
plt.rcParams.update({'font.size': 27})
plt.rcParams['lines.markersize'] = 12
plt.rcParams['lines.linewidth'] = 0.5

# Define the linear function for curve fitting
def func(x, m, c):
    """Linear function: y = m*x + c"""
    return m * x + c

def plot(filename, axs, axs_speed, model_label, axs_pred, axs_pred_iqr):
    """
    Plot range front evolution data for a given model.
    
    Parameters:
    -----------
    filename : str
        Directory containing the model data
    axs : matplotlib.axes.Axes
        Axes for plotting range front positions
    axs_speed : matplotlib.axes.Axes
        Axes for plotting expansion speeds
    model_label : list
        Model name/label for legend
    axs_pred : matplotlib.axes.Axes
        Axes for trajectories of individual simulation runs
    axs_pred_iqr : matplotlib.axes.Axes
        Axes for interquartile range (IQR) of range fronts
    """
    
    world_size = 144  # Size of the simulated world/grid
    print('Starting ' + filename + '/' + disp_stat)
    
    # Load aggregated statistics if available
    if os.path.isfile(filename + '/' + disp_stat + '/pf.csv'):
        df = pd.read_csv(filename + '/' + disp_stat + '/pf.csv', sep='\t')
        noz_time = df['time']
        evo_range_hot = df['range_hot']
        evo_low_hot = df['range_hot_low']
        evo_high_hot = df['range_hot_high']
        evo_range_cold = df['range_cold']
        evo_low_cold = df['range_cold_low']
        evo_high_cold = df['range_cold_high']
    
    hot_pf = []  # Store expansion speeds for hot front
    cold_pf = []  # Store expansion speeds for cold front
    
    # Process individual simulation runs
    for i in range(80):
        if os.path.isfile(filename + '/' + disp_stat + f'/pf_{i}.csv'):
            df_speed = pd.read_csv(filename + '/' + disp_stat + f'/pf_{i}.csv', sep='\t')
            noz_time_speed = df_speed['time']
            evo_range_hot_speed = df_speed['range_hot']
            evo_range_cold_speed = df_speed['range_cold']
            
            # Prepare data for linear regression (early time steps only)
            A = np.vstack([
                noz_time_speed[noz_time_speed < 2200],
                np.ones(len(noz_time_speed[noz_time_speed < 2200]))
            ]).T
            
            try:
                # Fit linear models to estimate expansion speeds
                m_hot, _ = curve_fit(func, A[:, 0], world_size - evo_range_hot_speed[noz_time_speed < 2200])
                m_cold, _ = curve_fit(func, A[:, 0], evo_range_cold_speed[noz_time_speed < 2200])
                
                hot_pf.append(m_hot[0])
                cold_pf.append(m_cold[0])
            except:
                continue
            
            # Plot individual trajectories
            axs_pred.scatter(noz_time_speed, evo_range_hot_speed, color=colors[0], alpha=0.3)
            axs_pred.scatter(noz_time_speed, evo_range_cold_speed, color=colors[0], alpha=0.3)
            axs_pred.set_xlabel('Time step')
            axs_pred.set_ylabel('Range front position')
    
    # Plot expansion speeds with IQR
    axs_speed.plot(['Cold'+'\n'+'front', 'Hot'+'\n'+'front'], 
                   [abs(np.median(cold_pf)), abs(np.median(hot_pf))], 
                   marker='o', color=colors[0], label=model_label[0])
    axs_speed.fill_between(['Cold'+'\n'+'front', 'Hot'+'\n'+'front'],
                          [abs(np.quantile(cold_pf, 0.25)), abs(np.quantile(hot_pf, 0.25))],
                          [abs(np.quantile(cold_pf, 0.75)), abs(np.quantile(hot_pf, 0.75))],
                          color=colors[0], alpha=0.4)
    axs_speed.set_ylim([0, 0.2])
    axs_speed.tick_params(axis='x', pad=10)
    
    # Plot range front evolution with confidence intervals
    axs.plot(noz_time, evo_range_hot, color=colors[0], marker='o', label=model_label[0])
    axs.fill_between(noz_time, evo_high_hot, evo_low_hot, color=colors[0], alpha=0.4)
    axs.plot(noz_time, evo_range_cold, color=colors[0], marker='o', markerfacecolor='none')
    axs.fill_between(noz_time, evo_high_cold, evo_low_cold, color=colors[0], alpha=0.4)
    axs.set_xlim([2000, 4000])
    axs.tick_params(axis='x', pad=10)
    axs.set_ylim([0, world_size])
    axs.set_xlabel('Time step')
    
    # Plot IQR of range fronts
    axs_pred_iqr.plot(noz_time, evo_high_hot - evo_low_hot, color=colors[0], marker='o', label=model_label[0])
    axs_pred_iqr.plot(noz_time, evo_high_cold - evo_low_cold, color=colors[0], marker='o', markerfacecolor='none')
    axs_pred_iqr.set_xlabel('Time step')
    axs_pred_iqr.set_ylabel('Range front IQR')
    axs_pred_iqr.set_xlim([2000, 4000])
    
    axs_pred.set_xlim([2000, 4000])

# ============================================================================
# MAIN PLOTTING SECTION
# ============================================================================

# Create figure for individual trajectories (A and B panels)
fig_pred, axs_pred = plt.subplots(1, 2)
axs_pred[0].text(0.95, 0.95, 'A', transform=axs_pred[0].transAxes, fontsize=30,
                 verticalalignment='top', horizontalalignment='left', fontweight='bold')
axs_pred[1].text(0.95, 0.95, 'B', transform=axs_pred[1].transAxes, fontsize=30,
                 verticalalignment='top', horizontalalignment='left', fontweight='bold')

# Create figure for evolution plots (2x2 grid)
fig_evo, axs_evo = plt.subplots(2, 2)
axs_evo[0][0].set_title('Without Thermal Adaptation', fontsize=25, loc='center')
axs_evo[0][1].set_title('With Thermal Adaptation', fontsize=25, loc='center')

# Label panels A-D
axs_evo[0][0].text(0.95, 0.95, 'A', transform=axs_evo[0][0].transAxes, fontsize=30,
                   verticalalignment='top', horizontalalignment='left', fontweight='bold')
axs_evo[0][1].text(0.95, 0.95, 'B', transform=axs_evo[0][1].transAxes, fontsize=30,
                   verticalalignment='top', horizontalalignment='left', fontweight='bold')
axs_evo[1][0].text(0.95, 0.95, 'C', transform=axs_evo[1][0].transAxes, fontsize=30,
                   verticalalignment='top', horizontalalignment='left', fontweight='bold')
axs_evo[1][1].text(0.95, 0.95, 'D', transform=axs_evo[1][1].transAxes, fontsize=30,
                   verticalalignment='top', horizontalalignment='left', fontweight='bold')

# ============================================================================
# ECOLOGICAL SIMULATIONS (Without Thermal Adaptation)
# ============================================================================
disp_stat = 'eco'  # Ecological simulations (no evolution)

# Plot Enzyme catalysis model
colors = np.array(['coral'])
labels = ['Enzyme catalysis based model']
plot('Enzyme_catalysis_model', axs_evo[0][0], axs_evo[1][0], labels, axs_pred[0], axs_pred[1])

# Plot Protein denaturation model
colors = np.array(['goldenrod'])
labels = ['Protein denaturation based model']
plot('Protein_denat_model', axs_evo[0][0], axs_evo[1][0], labels, axs_pred[0], axs_pred[1])

# Plot Gaussian model
colors = np.array(['silver'])
labels = ['Gauss model']
plot('Gauss_model', axs_evo[0][0], axs_evo[1][0], labels, axs_pred[0], axs_pred[1])

# ============================================================================
# EVOLUTIONARY SIMULATIONS (With Thermal Adaptation)
# ============================================================================
disp_stat = 'evo'  # Evolutionary simulations (with adaptation)

# Plot Enzyme catalysis model
colors = np.array(['coral'])
labels = ['Enzyme catalysis based model']
plot('Enzyme_catalysis_model', axs_evo[0][1], axs_evo[1][1], labels, axs_pred[0], axs_pred[1])

# Plot Protein denaturation model
colors = np.array(['goldenrod'])
labels = ['Protein denaturation based model']
plot('Protein_denat_model', axs_evo[0][1], axs_evo[1][1], labels, axs_pred[0], axs_pred[1])

# Plot Gaussian model
colors = np.array(['silver'])
labels = ['Gauss model']
plot('Gauss_model', axs_evo[0][1], axs_evo[1][1], labels, axs_pred[0], axs_pred[1])

# ============================================================================
# FINAL FIGURE FORMATTING
# ============================================================================
# Add labels and adjust layout
axs_evo[0][0].set_ylabel('Range front position')
axs_evo[1][0].set_ylabel('Initial expansion speed')
axs_evo[0][0].legend(loc=2, frameon=False, fontsize=22)

# Hide y-axis labels for right column panels
axs_evo[1][1].yaxis.set_ticklabels([])
axs_evo[0][1].yaxis.set_ticklabels([])

# Adjust figure size and layout
figure = fig_evo
figure.set_size_inches(20, 20)
fig_evo.tight_layout(pad=1)

# Save figure
fig_evo.savefig('Paper_figs/Fig_2_evo_eco_pf.png', dpi=500)
