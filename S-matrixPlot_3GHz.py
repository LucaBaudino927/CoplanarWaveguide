from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
from pathlib import Path
import numpy as np
#import scipy.special as special
#from scipy.special import ellipk as K
import mpmath as mp
from mpmath import ellipk as K
from typing import Sequence, Tuple, List

# --- CONFIGURATION ---
BASE_PATH_WAVE = Path('postpro/Rame/Wave/CPW')
BASE_PATH_LUMPED = Path('postpro/Rame/Lumped/CPW')
BASE_PATH_102OHM = Path('postpro/Rame/Lumped/CPW/Impedenza_102ohm')
BASE_PATH_102OHM_DIVERSE_CONDUCIBILITA = Path('postpro/Rame/Lumped/CPW/Impedenza_102ohm_mesh_22um_diverse_conducibilita')
BASE_PATH_TEST = Path('postpro/Rame/Lumped/CPW/Test_diversi_parametri_e_geometrie')
BASE_PATH_DIVERSE_PORTE = BASE_PATH_TEST / 'NO_PEC/Diverse_porte'
BASE_PATH_NUOVI_TEST_PEC = Path('postpro/Rame/CPW/NuoviTestPEC')
BASE_PATH_REPORT = Path('postpro/Rame/CPW/Report')

IMPEDANCES = [100, 102, 102.1, 102.5]
CONDUCTIVITIES = [33112582, 57471264, 59600000, 62500000, 72500000] #S/m
PORT_WIDTHS = [10, 50, 75, 100, 150, 200, 490]  #um
COLORS = ['green', 'red', 'blue', 'orange', 'cyan', 'black', 'magenta', 'purple', 'brown', 'pink', 'gray']
MARKERSIZE = 5
SIERRA_PATH = Path('postpro/Rame/Lumped/CPW/SierraData_122ohm')

CONDUCTIVITY_COPPER = 59600000  # S/m
CONDUCTIVITY_ALUMINIUM = 32894736.84 #S/m = 1/3.04x10-8  #33112582.78 #S/m 
SKIN_EFF_FRACTION = 0.5  # Fraction increase in conductor losses to consider skin effect relevant
Z0 = 105  # Characteristic impedance in Ohm
CPW_LENGTH = 100e-6
SIMULATION_CONDUCTIVITY = CONDUCTIVITY_ALUMINIUM

# --- FUNCTION DEFINITIONS ---


def read_csv_data(path, read_ports=None):
    """
    Reads S-parameter data from a CSV file.
    `read_ports` is a list of indices you want to read (e.g., [0, 1, 3, 5])
    """
    data = []
    with open(path, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            data.append([float(row[i]) for i in read_ports])
    return list(zip(*data))  # Returns list of columns


def plot_S11(xs, ys, y_err, skin_effect_thr, chi2vsSierra,                                                  #data
            labels, colors, dof,                                                                            #descriptions
            show_skin_eff_thr=True, show_ratio_plot=True, show_Sierra=True, show_error_bars= False,
            show_frequency_marker=False):                                                                   #options

    SierraData = []
    myDatasets = []
    if not show_error_bars: y_err = 0
    if show_ratio_plot and show_Sierra:
        fig, (ax, axRatio) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(15, 10))
        axRatio.set_xlabel('Frequency [GHz]', fontsize=15)
        axRatio.set_ylabel('Ratio', fontsize=15)
        axRatio.grid(True)
    else:
        fig, ax = plt.subplots(figsize=(15, 7))
        ax.set_xlabel('Frequency [GHz]', fontsize=15)

    for x, y, label, color in zip(xs, ys, labels, colors):
        if 'Sierra' in label:
            if show_Sierra: ax.errorbar(x, y, yerr=0, label=label, color='purple', linestyle='solid', marker='x', markersize=MARKERSIZE)
            SierraData = y
        else:
            #ax.plot(x, y, label=label, color=color, linestyle='solid', marker='o')
            ax.errorbar(x, y, xerr=0, yerr=y_err, color=color, linestyle='solid', marker='o', markersize=MARKERSIZE, label=label)
            myDatasets.append(y)
            if show_frequency_marker:
                for i, freq in enumerate(x):
                    if freq == 0.5:
                        ax.errorbar(x[i], y[i], xerr=0, yerr=0, color='red', marker='*', markersize=3*MARKERSIZE)
    ax.set_title('|S11|', fontsize=18)
    ax.grid(True)
    ax.set_ylabel('|S11| [dB]', fontsize=15)
    
    #ratio vs Sierra
    if show_ratio_plot and show_Sierra:
        loss_err_3_percent = np.array(np.multiply(SierraData, 0.03))
        ##plot "error bars" around sierra
        ax.errorbar(x, np.add(y,loss_err_3_percent), yerr=0, label='Sierra 3% Dispersion', color='purple', linestyle='--', marker='', markersize=MARKERSIZE)
        ax.errorbar(x, np.add(y,-loss_err_3_percent), yerr=0, color='purple', linestyle='--', marker='', markersize=MARKERSIZE)
        #Nel ratio plot stampo una sola banda di dispersione attorno a 1 tanto le barre definite come 0.03*ratio saranno tutte molto vicine a 0.97 e 1.03
        axRatio.hlines(1+0.03, min(x), max(x), color='purple', linestyle='--', label='Insertion Loss 3% Dispersion')
        axRatio.hlines(1-0.03, min(x), max(x), color='purple', linestyle='--')
        ##
        for i, myData in enumerate(myDatasets):
            ratio = np.array(myData) / np.array(SierraData)
            ratio_err = np.abs(ratio) * (y_err / np.abs(y))
            axRatio.errorbar(x, ratio, yerr=ratio_err, label=f'{labels[i]}/[Sierra 122Ohm]', color=COLORS[i], linestyle='solid', marker='o', markersize=MARKERSIZE)
        axRatio.axhline(1, color='black', linestyle='--')
        axRatio.set_ylim(0.5, 1.1)

    if skin_effect_thr > 0 and show_skin_eff_thr:
        ax.axvline(skin_effect_thr, color='green', linestyle='--', label=f'Skin Effect Threshold = {skin_effect_thr:.2f} GHz')
        if show_ratio_plot and show_Sierra: axRatio.axvline(skin_effect_thr, color='green', linestyle='--', label=f'Skin Effect Threshold = {skin_effect_thr:.2f} GHz')

    if show_Sierra:
        at = AnchoredText(describe_helper(labels, dof, chi2vsSierra=chi2vsSierra, chi2vsIPC=[]),
                        loc='upper left', prop=dict(size=15), frameon=True,
                        bbox_to_anchor=(0.96, 0.4), bbox_transform=ax.transAxes
                        )
        at.patch.set_boxstyle(boxstyle="round,pad=0.,rounding_size=0.2")
        at.patch.set_edgecolor((0, 0, 0, 0.2))
        ax.add_artist(at)
    
    handles, legend_labels = ax.get_legend_handles_labels()
    handles = [h[0] if isinstance(h, mpl.container.ErrorbarContainer) else h for h in handles]
    if show_Sierra:
        ax.legend(handles, legend_labels, loc='upper left', bbox_to_anchor=(0.96, 1), prop={'size': 15})
    else:
        ax.legend(handles, legend_labels, prop={'size': 15})
    if show_ratio_plot and show_Sierra:
        handles, legend_labels = axRatio.get_legend_handles_labels()
        handles = [h[0] if isinstance(h, mpl.container.ErrorbarContainer) else h for h in handles]
        axRatio.legend(handles, legend_labels, loc='upper left', bbox_to_anchor=(0.96, 1), prop={'size': 15})

    plt.xticks(rotation=25)
    plt.tight_layout()
    #plt.xscale('log')


def plot_S21(xs, ys, y_err, dielectric_losses, conductor_losses, skin_effect_thr, chi2vsSierra, chi2vsIPC,                  #data
            labels, colors, dof,                                                                                            #descriptions
            show_skin_eff_thr=True, show_loss=True, show_ratio_plot=True, show_Sierra=True, show_error_bars=False,
            show_frequency_marker=False):                                                                                   #options

    SierraData = []
    myDatasets = []
    if not show_error_bars: y_err = 0
    if show_ratio_plot:
        fig, (ax, axRatio) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(15, 10))
        axRatio.set_xlabel('Frequency [GHz]', fontsize=15)
        axRatio.set_ylabel('Ratio', fontsize=15)
        axRatio.grid(True)
    else:
        fig, ax = plt.subplots(figsize=(15, 7))
        ax.set_xlabel('Frequency [GHz]', fontsize=15)

    for x, y, label, color in zip(xs, ys, labels, colors):
        if 'Sierra' in label:
            if show_Sierra: ax.plot(x, y, label=label, color='purple', linestyle='solid', marker='x', markersize=MARKERSIZE)
            SierraData = y
        else:
            ax.errorbar(x, y, xerr=0, yerr=y_err, color=color, linestyle='solid', marker='o', markersize=MARKERSIZE, label=label)
            myDatasets.append(y)
            if show_frequency_marker:
                for i, freq in enumerate(x):
                    if freq == 0.5:
                        ax.errorbar(x[i], y[i], xerr=0, yerr=0, color='red', marker='*', markersize=3*MARKERSIZE)
    ax.set_title('|S21|', fontsize=18)
    ax.grid(True)
    ax.set_ylabel('|S21| [dB]', fontsize=15)
    
    #ratio vs Sierra
    if show_ratio_plot:
        if show_Sierra:
            for (i, myData), x in zip(enumerate(myDatasets), xs):
                ratio = np.array(myData) / np.array(SierraData)
                ratio_err = np.abs(ratio) * (y_err / np.abs(y))
                axRatio.errorbar(x, ratio, yerr=ratio_err, label=f'{labels[i]}/[Sierra 122Ohm]', color=COLORS[i+4], linestyle='solid', marker='o', markersize=MARKERSIZE)
        axRatio.axhline(1, color='black', linestyle='--')
        axRatio.set_ylim(0.8, 1.3)

    if skin_effect_thr > 0 and show_skin_eff_thr:
        ax.axvline(skin_effect_thr, color='green', linestyle='--', label=f'Skin Effect Threshold = {skin_effect_thr:.2f} GHz')
        if show_ratio_plot: axRatio.axvline(skin_effect_thr, color='green', linestyle='--', label=f'Skin Effect Threshold = {skin_effect_thr:.2f} GHz')

    if show_ratio_plot:
        axRatio.hlines(1+0.03, min(x), max(x), color='orange', linestyle='--', label='Insertion Loss 3% Dispersion')
        axRatio.hlines(1-0.03, min(x), max(x), color='orange', linestyle='--')

    if show_loss:
        first = True
        for (i, myData), x in zip(enumerate(myDatasets), xs):
            ratioIPC = np.divide(np.array(myData), -np.array(np.add(dielectric_losses, conductor_losses)))
            ratio_err = y_err / np.array(np.add(dielectric_losses, conductor_losses))
            if show_ratio_plot: 
                axRatio.errorbar(x, ratioIPC, yerr=ratio_err, label=f'{labels[i]}/IPC', color=colors[i], linestyle='solid', marker='^', markersize=MARKERSIZE)
                if show_frequency_marker:
                    for i, freq in enumerate(x):
                        if freq == 0.5:
                            axRatio.errorbar(x[i], ratioIPC[i], xerr=0, yerr=0, color='red', marker='*', markersize=3*MARKERSIZE)

        loss_err_3_percent = np.array(np.multiply(np.add(dielectric_losses, conductor_losses), 0.03))
        for i in range(len(xs[0]) - 1):
            ax.hlines(-dielectric_losses[i], x[i], x[i+1], color='blue', linestyle='--', label='IPC Dielectric Loss')
            ax.fill_betweenx([-dielectric_losses[i], 0], x[i], x[i+1], color='blue', alpha=0.1)
            ax.hlines(-dielectric_losses[i]-conductor_losses[i], x[i], x[i+1], color='red', linestyle='--', label='IPC Insertion Loss')
            ax.fill_betweenx([-dielectric_losses[i]-conductor_losses[i], -dielectric_losses[i]], x[i], x[i+1], color='red', alpha=0.1)
            ax.hlines(-dielectric_losses[i]-conductor_losses[i]-loss_err_3_percent[i], x[i], x[i+1], color='orange', linestyle='--', label='Insertion Loss 3% Dispersion')
            ax.hlines(-dielectric_losses[i]-conductor_losses[i]+loss_err_3_percent[i], x[i], x[i+1], color='orange', linestyle='--')
            '''
            ax.hlines(-dielectric_loss_thr[i]-conductor_loss_thr[i]-(calculate_conductor_losses(1e-9)*(SKIN_EFF_FRACTION)), x[i], x[i+1], color='yellow', linestyle='--', label=f'Skin effect losses < {SKIN_EFF_FRACTION*100}%')
            ax.fill_betweenx([-dielectric_loss_thr[i]-conductor_loss_thr[i]-(calculate_conductor_losses(1e-9)*(SKIN_EFF_FRACTION)), 
                            -dielectric_loss_thr[i]-conductor_loss_thr[i]], 
                            x[i], x[i+1], color='yellow', alpha=0.1)
            '''

            #set legend only one time
            if first:
                first = False
                at = AnchoredText(describe_helper(labels, dof, chi2vsSierra=chi2vsSierra, chi2vsIPC=chi2vsIPC),
                                  loc='upper left', prop=dict(size=15), frameon=True,
                                  bbox_to_anchor=(0.96, 0.55), bbox_transform=ax.transAxes
                                  )
                at.patch.set_boxstyle(boxstyle="round,pad=0.,rounding_size=0.2")
                at.patch.set_edgecolor((0, 0, 0, 0.2))
                ax.add_artist(at)
                handles, legend_labels = ax.get_legend_handles_labels()
                handles = [h[0] if isinstance(h, mpl.container.ErrorbarContainer) else h for h in handles]
                ax.legend(handles, legend_labels, loc='upper left', bbox_to_anchor=(0.96, 1), fontsize='large', prop={'size': 15})
                if show_ratio_plot:
                    handles, legend_labels = axRatio.get_legend_handles_labels()
                    handles = [h[0] if isinstance(h, mpl.container.ErrorbarContainer) else h for h in handles]
                    axRatio.legend(handles, legend_labels, loc='upper left', bbox_to_anchor=(0.96, 1.1), fontsize='large', prop={'size': 15})

    if not show_loss:
        at = AnchoredText(describe_helper(labels, dof, chi2vsSierra=chi2vsSierra, chi2vsIPC=chi2vsIPC),
                          loc='center left', prop=dict(size=10), frameon=True,
                          bbox_to_anchor=(0.96, 0.5), bbox_transform=ax.transAxes
                          )
        at.patch.set_boxstyle(boxstyle="round,pad=0.,rounding_size=0.2")
        at.patch.set_edgecolor((0, 0, 0, 0.2))
        ax.add_artist(at)
        ax.legend(loc='upper left', bbox_to_anchor=(0.96, 1), fontsize='large', prop={'size': 15})
        if show_ratio_plot:
            axRatio.legend(loc='upper left', bbox_to_anchor=(0.96, 1), fontsize='large', prop={'size': 15})
 
    plt.xticks(rotation=25)
    plt.tight_layout()
    #plt.xscale('log')
    #plt.subplots_adjust(right=0.9)


def describe_helper(labels, dof, chi2vsSierra=None, chi2vsIPC=None):
    value = ""
    if chi2vsSierra:
        value += f"$\\chi^2$ vs Sierra with {dof} dof\n"
        for label, chi2 in zip(labels, chi2vsSierra):
            value += f"{label}  $\\chi^2$ = {chi2:.2f}\n"
    if chi2vsIPC:
        value += f"$\\chi^2$ vs IPC with {dof} dof\n"
        for label, chi2 in zip(labels, chi2vsIPC):
            value += f"{label}  $\\chi^2$ = {chi2:.2f}\n"
    #remove last \n
    value = value[:-2]
    return value


def calculate_skin_effect_threshold(x = SKIN_EFF_FRACTION, conductivity=SIMULATION_CONDUCTIVITY, trace_length=CPW_LENGTH, trace_width=90e-6, trace_thickness=20e-6, Z0=Z0):
    
    """
    Metodo 1:
    Calculate the skin effect threshold frequency for a given conductivity and trace properties.
    We consider the skin effect relevant when the resistence become x higher than ohmic resistance. Default x=0.1 means 10% increase.
    """
    k = np.sqrt(4 * 1e-7 * np.pi * np.pi)
    perimeter = 2 * (trace_thickness + trace_width)
    R = trace_length / (conductivity * trace_thickness * trace_width)
    frequency = conductivity * np.pow(R * (x + 1), 2) * np.pow(perimeter, 2) / (np.pow(trace_length, 2) * np.pow(k, 2))

    #return frequency/1e9  # Convert to GHz

    """
    Metodo 2:
    Calculate the skin effect threshold frequency for a given conductivity and trace properties.
    We consider the skin effect relevant when the conductor loss become x higher than conductor loss at low frequency. Default x=0.1 means 10% increase.
    """
    #formula vecchia
    frequency = np.pow(x + 1, 2) / (conductivity * np.pi * 4e-7 * np.pi *np.pow(trace_thickness, 2))

    #formula nuova
    #frequency = np.pow(x + 1, 2) / (conductivity * np.pi * 4e-7 * np.pi *np.pow(trace_thickness, 2)) * np.pow(trace_length / trace_width, 2)

    """
    Metodo 3:
    Calculate the skin effect threshold frequency for a given conductivity and trace properties.
    We consider the skin effect relevant when the conductor loss become a factor N times higher than conductor loss at low frequency.
    """
    N = 14
    conductor_loss_low_freq = calculate_conductor_losses(1e-9, conductivity, trace_width, trace_thickness, trace_length, 100e-6, Z0)  # dB
    frequency = np.pow((N*conductor_loss_low_freq)/trace_length,2)*np.pow(2*Z0*conductivity*trace_width/8.686,2)/(np.pi*4e-7*np.pi*conductivity)
    mu0 = 4e-7 * np.pi
    delta = np.sqrt(1 / (np.pi * frequency * mu0 * conductivity))  # Skin depth
    print(f"Skin effect threshold frequency: {frequency/1e9} GHz, skin depth delta = {delta*1e6} um at N = {N}")

    return float(frequency)/1e9  # GHz


def calculate_dielectric_losses(frequency=3, trace_length=CPW_LENGTH, trace_width=90e-6, gap_width=100e-6, trace_thickness=20e-6, dielectric_constant=3.3, loss_tangent=0.0013):
    """
    Fraction of power loss in the dielectric for a CPW line.
    It calculates Ad in dB/m and return the value in dB.
    """
    def KoverKprime(k, kprime):
        s = lambda x: np.log(2 * (np.sqrt(1 + x) + np.pow(4 * x, 1 / 4)) / (np.sqrt(1 + x) - np.pow(4 * x, 1 / 4)))
        if k >= 1.0 / np.sqrt(2):
            return s(k) / (2 * np.pi)
        else:
            return 2 * np.pi / s(kprime)

    k = trace_width / (trace_width + 2 * gap_width)
    k1 = np.sinh(np.pi * trace_width / (4 * trace_thickness)) / np.sinh(np.pi * (trace_width + 2 * gap_width) / (4 * trace_thickness))

    kprime = np.sqrt(1 - np.pow(k, 2))
    k1prime = np.sqrt(1 - np.pow(k1, 2))

    koverkprime = KoverKprime(k, kprime)
    k1overk1prime = KoverKprime(k1, k1prime)

    e_eff = 1 + ((dielectric_constant - 1) / 2) * k1overk1prime / koverkprime

    # Dielectric loss in dB/m
    #with e_eff
    Ad = (8.686 * np.pi) * np.sqrt(e_eff) * loss_tangent * frequency*1e9 / 2.99792458e8
    #without e_eff
    Ad = (8.686 * np.pi) * np.sqrt(dielectric_constant) * loss_tangent * frequency*1e9 / 2.99792458e8
    print(f"Frequency: {frequency:.2f} GHz, Effective Dielectric Constant: {e_eff:.4f}, Dielectric Loss: {1e-2 * Ad} dB/cm,  Dielectric Loss: {Ad * trace_length} dB")

    return float(Ad) * (trace_length)  # Convert to dB for the given trace length in meters


def calculate_conductor_losses(frequency=3, conductivity=SIMULATION_CONDUCTIVITY, trace_width=90e-6, trace_thickness=20e-6, trace_length=CPW_LENGTH, gap_width=100e-6, Z0=Z0):
    """
    Fraction of power loss in the conductor for a CPW line.
    It calculates Ac in dB/m and return the value in dB.
    """
    mu0 = 4e-7 * np.pi
    delta = np.sqrt(1 / (np.pi * frequency*1e9 * mu0 * conductivity))  # Skin depth

    Rs = (1 + 1j) / (conductivity * delta) # Surface resistance in Ohm
    Rstrip = Rs.real / (trace_width) # Ohm/m
    # Conductor loss in dB/m
    if delta > trace_thickness:
        Ac = (8.686) / (2 * Z0 * conductivity * trace_thickness * trace_width)
    else:
        Ac = (8.686 * Rstrip) / (2 * Z0)

    
    '''
    Test dal paper: "A new analytical, cad-oriented model for the ohmic and radiation loesses of asymmetric coplanar waveguides in
    hybrid and monolithic mic's", G. Ghione, C.U. Naldi
    '''
    if delta > trace_thickness: Rs = (1 + 1j) / (conductivity * trace_thickness)
    a = trace_width / 2
    b = a + gap_width
    ks = a / b
    ksprime = np.sqrt(1 - ks*ks)
    k2 = np.sinh(np.pi * a / (2 * 25e-6)) / np.sinh(np.pi * b / (2 * 25e-6))
    k2prime = np.sqrt(1 - k2*k2)
    #per il calcolo di e_eff provare la formula del paper An Equivalent Circuit Model for the Coplanar Waveguide Step Discontinuity
    #Colin Sinclair and Stephen J Nightingale
    e_eff = 1 + ((3.3 - 1) / 2) * (K(ksprime)/K(ks)) * (K(k2)/K(k2prime)) 
    #anche questa è un'approssimazione per e_eff che non va bene, chiedere ad Alessandro se ha accesso al paper
    #G. Ghione and C. Naldi, “Analytical formulas for coplanar lines in hybrid and monolithic MIC’s,” Electron. Lett., vol. 20, pp. 179–181, 1984


    '''
    k = np.sqrt(2*a*2*b/((a+b)*(a+b)))
    kprime = np.sqrt(1 - k*k)
    k1 = np.sqrt( (2*np.sinh(np.pi*a/(2*25e-6))*2*np.sinh(np.pi*b/(2*25e-6))) / ((np.sinh(np.pi*a/(2*25e-6))+np.sinh(np.pi*b/(2*25e-6)))*(np.sinh(np.pi*a/(2*25e-6))+np.sinh(np.pi*b/(2*25e-6)))) )
    k1prime = np.sqrt(1 - k1*k1)
    e_eff = 1 + ((3.3 - 1) / 2) * (K(k1prime)/K(k1)) * (K(k)/K(kprime))
    '''

    Ac = 8.686 * Rs.real * np.sqrt(3.3) / (480 * np.pi * K(ks) * K(ksprime) * (1 - ks*ks)) * (1/a * (np.pi + np.log(8 * np.pi * a * (1 - ks) / (trace_thickness * (1 + ks)))) + 
                                                                                                1/b * (np.pi + np.log(8 * np.pi * b * (1 - ks) / (trace_thickness * (1 + ks)))))
    
    print(f"Frequency: {frequency:.2f} GHz, Skin Depth: {delta*1e6:.4f} um, e_eff: {e_eff}, Rs: {Rs.real:.4f} Ohm, Conductor Loss: {1e-2 * Ac} dB/cm, Conductor Loss: {Ac * trace_length} dB")

    return float(Ac) * (trace_length)  # Convert to dB for the given trace length in meters


def calculate_free_space_radiation_losses(frequency=3, dielectric_constant=3.3, trace_width=90e-6, gap_width=90e-6, trace_length=CPW_LENGTH):
    """
    Fraction of power loss due to free-space radiation.
    It calculates Afsr in dB/m and return the value in dB.
    """

    '''
    Test dal paper: "Terahertz attenuation and dispersion characteristics of coplanar transmission lines"
    M.Y. Frankel;S. Gupta;J.A. Valdmanis;G.A. Mouron
    '''
    a = trace_width / 2
    b = a + gap_width
    ks = a / b
    ksprime = np.sqrt(1 - ks*ks)
    k2 = np.sinh(np.pi * a / (2 * 25e-6)) / np.sinh(np.pi * b / (2 * 25e-6))
    k2prime = np.sqrt(1 - k2*k2)
    e_eff = 1 + ((3.3 - 1) / 2) * (K(ksprime)/K(ks)) * (K(k2)/K(k2prime))
    k = trace_width / (trace_width + 2 * gap_width)
    kprime = np.sqrt(1 - k*k)
    Afsr = (8.686 * 2 * np.pow(np.pi/2, 5) * np.pow(1 - e_eff/dielectric_constant, 2) * np.pow(trace_width + 2*gap_width, 2) * np.pow(dielectric_constant, 1.5) * np.pow(frequency*1e9, 3)) / (np.sqrt(e_eff/dielectric_constant) * np.pow(2.99792458e8, 3) * K(kprime) * K(k))

    print(f'K(kprime): {float(K(kprime))}, K(k): {float(K(k))}')
    print(f'Frequency: {frequency:.2f} GHz, Free-Space Radiation Loss: {1e-2 * Afsr} dB/cm, Free-Space Radiation Loss: {Afsr * trace_length} dB')
    return float(Afsr) * (trace_length)  # Convert to dB for the given trace length in meters


def calculate_Z0_from_conductor_losses(trace_width=90e-6, trace_thickness=20e-6, gap_width=100e-6, substrate_thickness=25e-6, dielectric_constant=3.3):
    """
    Calculate Z0 from the predicted conductor loss using following paper IPC formula and comparing it with the usual conductor loss equation.
    "A new analytical, cad-oriented model for the ohmic and radiation loesses of asymmetric coplanar waveguides in
    hybrid and monolithic mic's", G. Ghione, C.U. Naldi
    """

    a = trace_width / 2
    b = a + gap_width
    ks = a / b
    ksprime = np.sqrt(1 - ks*ks)
    k2 = np.sinh(np.pi * a / (2 * substrate_thickness)) / np.sinh(np.pi * b / (2 * substrate_thickness))
    k2prime = np.sqrt(1 - k2*k2)
    e_eff = 1 + ((dielectric_constant - 1) / 2) * (K(ksprime)/K(ks)) * (K(k2)/K(k2prime))
    #Ac = 8.686 * Rs.real * np.sqrt(dielectric_constant) / (480 * np.pi * K(ks) * K(ksprime) * (1 - ks*ks)) * (1/a * (np.pi + np.log(8 * np.pi * a * (1 - ks) / (trace_thickness * (1 + ks)))) + 
    #                                                                                            1/b * (np.pi + np.log(8 * np.pi * b * (1 - ks) / (trace_thickness * (1 + ks)))))
    #print(f'CPW impedance v1: {float((30 * np.pi) / np.sqrt(dielectric_constant) * (K(ksprime)/K(ks))):.2f} Ohm')
    Z0 = (1/trace_width) * (240 * np.pi * K(ks) * K(ksprime) * (1 - ks*ks)) / (np.sqrt(dielectric_constant) * (1/a * (np.pi + np.log(8 * np.pi * a * (1 - ks) / (trace_thickness * (1 + ks)))) + 1/b * (np.pi + np.log(8 * np.pi * b * (1 - ks) / (trace_thickness * (1 + ks))))))

    #versione approssimata, infatti non funziona
    #print(f'CPW impedance v1: {float((30 * np.pi) / np.sqrt(3.3) * (K(ksprime)/K(ks))):.2f} Ohm')

    print(f'CPW impedance: {Z0} Ohm')

    return float(Z0)


def plot_Z0_vs_software():
    """
    Plot the characteristic impedance Z0 calculated by different software.
    """
    # Data from different software
    widths = np.array([50, 60, 70, 80, 90, 100])  # in micrometers
    gap = np.array([50, 60, 70, 80, 90, 100])  # in micrometers
    #z0_ansys_gap90um = np.array([186.5, 142.3, 112.5, 97.3, 86.9, 73.2, 64.8, 54.3, 47.7, 42.9])
    z0_sierra_gap100um = np.array([143.86789, 136.85319, 131.23733, 126.58972, 122.64803, 119.24168])
    z0_txline_gap100um = np.array([122.341, 118.392, 114.858, 111.948, 109.388, 107.107])
    z0_sierra_trace90um = np.array([86.59548, 94.99007, 102.67336, 109.78425, 116.41895, 122.64803])
    z0_txline_trace90um = np.array([80.977, 87.8673, 93.9013, 99.5041, 104.642, 109.388])

    fig, ax = plt.subplots()
    #plt.plot(widths, z0_ansys, label='Ansys HFSS', marker='o')
    ax.plot(widths, z0_sierra_gap100um, label='Sierra', marker='s')
    ax.plot(widths, z0_txline_gap100um, label='TXLINE', marker='d')
    ax.errorbar(90, 105, xerr=0, yerr=2, linestyle='None', marker='^', label='Palace')
    ax.set_title('Characteristic Impedance Z0 vs Trace Width @100um Gap Width', fontsize=18)
    ax.grid(True)
    ax.legend(fontsize='x-large')
    plt.xlabel('Trace Width [um]')
    plt.ylabel('Characteristic Impedance Z0 [Ohm]')
    plt.xticks(rotation=25)
    plt.tight_layout()

    fig, ax = plt.subplots()
    #plt.plot(widths, z0_ansys, label='Ansys HFSS', marker='o')
    ax.plot(gap, z0_sierra_trace90um, label='Sierra', marker='s')
    ax.plot(gap, z0_txline_trace90um, label='TXLINE', marker='d')
    ax.errorbar(100, 105, xerr=0, yerr=2, linestyle='None', marker='^', label='Palace')
    ax.set_title('Characteristic Impedance Z0 vs Gap Width @90um Trace Width', fontsize=18)
    ax.grid(True)
    ax.legend(fontsize='x-large')
    plt.xlabel('Gap Width [um]')
    plt.ylabel('Characteristic Impedance Z0 [Ohm]')
    plt.xticks(rotation=25)
    plt.tight_layout()

    plt.show()


def chisq(obs, exp, error):
    return np.sum(np.pow(obs - exp, 2) / np.pow(error, 2))


def estimate_systematic_uncertainty():

    s11_all = []
    s21_all = []
    s11_errors = []
    s21_errors = []
    paths = ['DiverseImpedenze/DiverseImpedenze50MHz/cpw_lumped_22um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
             'DiverseImpedenze/DiverseImpedenze50MHz/cpw_lumped_22um_104ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
             'DiverseImpedenze/DiverseImpedenze50MHz/cpw_lumped_22um_106ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
             'DiverseMesh/DiverseMesh50MHz/cpw_lumped_15um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
             #'DiverseMesh/cpw_lumped_15.0um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
             #'DiverseConducibilita/cpw_lumped_22.5um_105ohm_cond_58000000_porte_75um_NOPEC_Absorbing_4'
            ]
    for i, path in enumerate(paths):
        path = BASE_PATH_REPORT / f'{path}/port-S.csv'
        data = read_csv_data(path, [1, 3])
        s11_all.append(data[0])
        s21_all.append(data[1])

    for index_freq in range(0, len(s11_all[0])):
        s11_fixed_frequency = []
        s21_fixed_frequency = []
        for s11, s21 in zip(s11_all, s21_all):
            s11_fixed_frequency.append(s11[index_freq])
            s21_fixed_frequency.append(s21[index_freq])
        s11_errors.append(abs((max(s11_fixed_frequency) - min(s11_fixed_frequency)) / 2))
        s21_errors.append(abs((max(s21_fixed_frequency) - min(s21_fixed_frequency)) / 2))
        s11_fixed_frequency.clear()
        s21_fixed_frequency.clear()

    return s11_errors, s21_errors

def calculate_Z_and_propagation_constant_from_S_matrix(s11_all, s11_phases_all, s21_all, s21_phases_all, Z0=Z0, trace_length=CPW_LENGTH):
    
    def back_from_db(x_db: float) -> float:
        return 10.0**(x_db / 20.0)

    def deg2rad(x_deg: float) -> float:
        return np.deg2rad(x_deg)

    def polar_to_complex(r: float, phi: float) -> complex:
        return r * np.exp(1j * phi)

    # --- Build S11 complex datasets --------------------------
    s11_complex_all: List[List[complex]] = []
    for magnitudes_db, phases_deg in zip(s11_all, s11_phases_all):
        if not np.isscalar(magnitudes_db):
            s11_complex_list = []
            for m_db, p_deg in zip(magnitudes_db, phases_deg):
                s11_complex_list.append( polar_to_complex(back_from_db(m_db), deg2rad(p_deg)) )
            s11_complex_all.append(s11_complex_list)
        else:
            s11_complex_all.append( polar_to_complex(back_from_db(magnitudes_db), deg2rad(phases_deg)) )
    #print(f's11_complex_all={s11_complex_all}')

    # --- Build S21 complex datasets --------------------------
    s21_complex_all: List[List[complex]] = []
    for magnitudes_db, phases_deg in zip(s21_all, s21_phases_all):
        if not np.isscalar(magnitudes_db):
            s21_complex_list = []
            for m_db, p_deg in zip(magnitudes_db, phases_deg):
                s21_complex_list.append( polar_to_complex(back_from_db(m_db), deg2rad(p_deg)) )
            s21_complex_all.append(s21_complex_list)
        else:
            s21_complex_all.append( polar_to_complex(back_from_db(magnitudes_db), deg2rad(phases_deg)) )
    #print(f's21_complex_all={s21_complex_all}')

    Z: List[List[complex]] = []
    gamma: List[List[complex]] = []
    for s11_dataset, s21_dataset in zip(s11_complex_all, s21_complex_all):
        if not np.isscalar(s11_dataset):
            Z_list = []
            gamma_list = []
            for s11, s21 in zip(s11_dataset, s21_dataset):
                Z_list.append( np.sqrt( Z0*Z0 * (np.pow(1 + s11, 2) - np.pow(s21, 2)) / (np.pow(1 - s11, 2) - np.pow(s21, 2))) )
                gamma_list.append( np.log(( (1 - np.pow(s11, 2) + np.pow(s21, 2)) / (2*s21) ) + np.sqrt(( (np.pow(s11, 2) - np.pow(s21, 2) + 1) - np.pow(2*s11, 2)) / (np.pow(2*s21, 2)))) / trace_length
            )
            Z.append(Z_list)
            gamma.append(gamma_list)
        else:
            Z.append(
                np.sqrt( Z0*Z0 * (np.pow(1 + s11_dataset, 2) - np.pow(s21_dataset, 2)) / (np.pow(1 - s11_dataset, 2) - np.pow(s21_dataset, 2)) )
            )
            gamma.append(
                np.log(( (1 - np.pow(s11_dataset, 2) + np.pow(s21_dataset, 2)) / (2*s21_dataset) ) + np.sqrt(( (np.pow(s11_dataset, 2) - np.pow(s21_dataset, 2) + 1) - np.pow(2*s11_dataset, 2)) / (np.pow(2*s21_dataset, 2)))) / trace_length
            )
    #print(f'Z={Z}')        
    #print(f'propagation constant={gamma}')

    return Z, gamma


def calculate_RGLC_and_plot_Z0(Z, gamma, x_all, labels, colors=COLORS):

    R: List[List[complex]] = []
    G: List[List[complex]] = []
    L: List[List[complex]] = []
    C: List[List[complex]] = []
    Z0_calc: List[List[complex]] = []
    for Z_dataset, gamma_dataset, freq_dataset in zip(Z, gamma, x_all):
        if not np.isscalar(Z_dataset):
            R_list = []
            G_list = []
            L_list = []
            C_list = []
            Z0_list =[]
            for Z_i, gamma_i, freq_i in zip(Z_dataset, gamma_dataset, freq_dataset):
                freq_i = freq_i*2*np.pi
                R_list.append((gamma_i*Z_i).real)
                G_list.append((gamma_i/Z_i).real)
                L_list.append(((gamma_i*Z_i).imag)/(freq_i*1e9))
                C_list.append(((gamma_i/Z_i).imag)/(freq_i*1e9))
                Z0_list.append( np.sqrt((((gamma_i*Z_i).real)+1j*(freq_i*1e9)*(((gamma_i*Z_i).imag)/(freq_i*1e9)))/(((gamma_i/Z_i).real)+1j*(freq_i*1e9)*(((gamma_i/Z_i).imag)/(freq_i*1e9)))).real )
            R.append(R_list)
            G.append(G_list)
            L.append(L_list)
            C.append(C_list)
            Z0_calc.append(Z0_list)
        else:
            freq_dataset = freq_dataset*2*np.pi
            R.append((gamma_dataset*Z_dataset).real)
            G.append((gamma_dataset/Z_dataset).real)
            L.append(((gamma_dataset*Z_dataset).real)/(freq_dataset*1e9))
            C.append(((gamma_dataset/Z_dataset).real)/(freq_dataset*1e9))
            Z0_calc.append( np.sqrt((((gamma_dataset*Z_dataset).real)+1j*(freq_dataset*1e9)*(((gamma_dataset*Z_dataset).imag)/(freq_dataset*1e9)))/(((gamma_dataset/Z_dataset).real)+1j*(freq_dataset*1e9)*(((gamma_dataset/Z_dataset).imag)/(freq_dataset*1e9)))).real )

    fig, ax = plt.subplots()
    if not np.isscalar(Z0_calc[0]):
        for x, y, label, color in zip(x_all, Z0_calc, labels, colors):
            ax.errorbar(x, y, xerr=0, yerr=0, color=color, linestyle='solid', marker='o', markersize=MARKERSIZE, label=label)
            for i, freq in enumerate(x):
                if freq == 0.5:
                    ax.errorbar(x[i], y[i], xerr=0, yerr=0, color='red', marker='*', markersize=3*MARKERSIZE)
    else:
        #ax.errorbar(x_all, Z0_calc, xerr=0, yerr=0, color=colors[0], linestyle='solid', marker='o', markersize=MARKERSIZE, label='Z0')
        #ax.errorbar(x_all, R,       xerr=0, yerr=0, color=colors[1], linestyle='solid', marker='o', markersize=MARKERSIZE, label='R')
        #ax.errorbar(x_all, G,       xerr=0, yerr=0, color=colors[2], linestyle='solid', marker='o', markersize=MARKERSIZE, label='G')
        ax.errorbar(x_all, L,       xerr=0, yerr=0, color=colors[3], linestyle='solid', marker='o', markersize=MARKERSIZE, label='L')
        ax.errorbar(x_all, C,       xerr=0, yerr=0, color=colors[4], linestyle='solid', marker='o', markersize=MARKERSIZE, label='C')

    ax.set_title('Characteristic Impedance Z0 vs Frequency', fontsize=18)
    ax.grid(True)
    ax.legend(fontsize='x-large')

    
    fig1, ax1 = plt.subplots()
    expected_s11 = []
    print(f'len(Z0_calc)={len(Z0_calc)}, len(x_all)={len(x_all)}')
    for x, y in zip(x_all, Z0_calc):
        print(f'20*np.log(np.abs((Z0 - y)/(Z0 + y)))={20*np.log(np.abs((Z0 - y)/(Z0 + y)))}, Z0={Z0}, y={y}, x={x}')
        expected_s11.append( 20*np.log(1/(1 - (np.abs((Z0 - y)/(Z0 + y)))*(np.abs((Z0 - y)/(Z0 + y))) )) )
    for x, y, label, color in zip(x_all, expected_s11, labels, colors):
        ax1.errorbar(x, y, xerr=0, yerr=0, linestyle='solid', marker='o', markersize=MARKERSIZE)
    ax1.set_title('Expected S11 vs Frequency', fontsize=18)
    ax1.grid(True)
    ax1.legend(fontsize='x-large')
    

    plt.xlabel('Frequency [GHz]')
    plt.ylabel('Characteristic Impedance Z0 [Ohm]')
    plt.xticks(rotation=25)
    plt.tight_layout()

    return Z0_calc[0]


def calculate_RGLC_and_plot_Z0_from_the_105Ohm_simulation(show_plot=True):

    #Questo era per dare un'idea dei valori di RLGC ma ora dovrei partire da dei valori fissi per la nostra linea senza andare
    #a leggere dalla simulazione per ricostruirli
    path = BASE_PATH_REPORT / 'DiverseMesh/DiverseMesh50MHz/cpw_lumped_22um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4/port-S.csv'
    x, s11, s11_phase, s21, s21_phase = read_csv_data(path, [0, 1, 2, 3, 4])
    labels = ['[$Z_{S,L}=105\\Omega$]']
    Z, gamma = calculate_Z_and_propagation_constant_from_S_matrix(s11, s11_phase, s21, s21_phase, Z0, trace_length=CPW_LENGTH)
    Z0_calc = calculate_RGLC_and_plot_Z0(Z, gamma, x, labels[0])
    if show_plot: plt.show()
    return Z0_calc



def plot_single_simulation(path, skin_effect_thr, label, sierra_path=SIERRA_PATH):
    x, s11, s11_phase, s21, s21_phase = read_csv_data(path, [0, 1, 2, 3, 4])
    #s11_err, s21_err = estimate_systematic_uncertainty()
    x_sierra, sierra_s11, sierra_s21 = read_csv_data(sierra_path, [0, 1, 3])

    dielectric_losses = []
    conductor_losses = []
    for freq in x:
        dielectric_losses.append(calculate_dielectric_losses(freq))
        conductor_losses.append(calculate_conductor_losses(freq))

    plot_S11([x, x_sierra], [s11, sierra_s11], 0, skin_effect_thr, [],
            [label, '[Sierra 122Ohm]'], ['blue', 'r'], 0,
            False, True, False, False, True)
    
    chi2vsIPC = []
    print("\n===========Chi-square S21 vs IPC with errors on IPC value===========")
    print("Cutting frequencies over skin effect threshold")
    loss_err_2_percent = np.array(np.multiply(np.add(dielectric_losses, conductor_losses), 0.02))
    loss_err_3_percent = np.array(np.multiply(np.add(dielectric_losses, conductor_losses), 0.03))
    loss_err_4_percent = np.array(np.multiply(np.add(dielectric_losses, conductor_losses), 0.04))
    
    freq_index_after_thr = len(s21)
    if skin_effect_thr < 2.9:
        freq_index_after_thr = next(j for j, freq in enumerate(x) if freq > skin_effect_thr)
    myData = s21[:freq_index_after_thr]
    chi2 = chisq(-np.array(np.add(dielectric_losses[:freq_index_after_thr], conductor_losses[:freq_index_after_thr])), np.array(myData[:freq_index_after_thr]), loss_err_2_percent[:freq_index_after_thr])
    print(f"S21 Chi-squared with {len(myData) - 2} dof, assuming dispersion from IPC of 2% for {label} vs IPC: {chi2}")
    chi2 = chisq(-np.array(np.add(dielectric_losses[:freq_index_after_thr], conductor_losses[:freq_index_after_thr])), np.array(myData[:freq_index_after_thr]), loss_err_3_percent[:freq_index_after_thr])
    chi2vsIPC.append(chi2) #per ora salvo solo quello al 3% perchè sembra il risultato migliore
    dof = len(myData) - 2
    print(f"S21 Chi-squared with {len(myData) - 2} dof, assuming dispersion from IPC of 3% for {label} vs IPC: {chi2}")
    chi2 = chisq(-np.array(np.add(dielectric_losses[:freq_index_after_thr], conductor_losses[:freq_index_after_thr])), np.array(myData[:freq_index_after_thr]), loss_err_4_percent[:freq_index_after_thr])
    print(f"S21 Chi-squared with {len(myData) - 2} dof, assuming dispersion from IPC of 4% for {label} vs IPC: {chi2}")

    plot_S21([x, x_sierra], [s21, sierra_s21], 0, dielectric_losses, conductor_losses, skin_effect_thr, [], chi2vsIPC,
            [label, '[Sierra 122Ohm]'], ['blue', 'r'], dof,
            True, True, True, False, False, True)

    #calculate chi squared for S11
    #chi2 = chisq(np.array(s11[:]), -np.array(np.add(dielectric_losses[:], conductor_losses[:])), s21_err[:])
    #print(f"S21 Chi-squared with {len(s11) - 1} dof for {label} vs IPC: {chi2}")

    plt.show()


def plot_multiple_simulations(base_path):
    x_all = []
    s11_all = []
    s11_phases_all = []
    s21_all = []
    s21_phases_all = []
    labels = []
    sierra_path = ''
    skin_effect_thr = -1

    if base_path == BASE_PATH_102OHM:
        sierra_path = SIERRA_PATH / 'port-S_10MHz_3GHz.csv'
        labels = [str(imp) for imp in IMPEDANCES]
        for i, imp in enumerate(IMPEDANCES):
            path = base_path / f'coplanar_waveguide_lumped_22um_{imp}ohm_cond_62500000/port-S.csv'
            data = read_csv_data(path, [1, 3])
            s11_all.append(data[0])
            s21_all.append(data[1])
        # Calculate skin effect threshold
        skin_effect_thr = calculate_skin_effect_threshold(SKIN_EFF_FRACTION, 62500000, 1000e-6, 90e-6, 20e-6)


    if base_path == BASE_PATH_DIVERSE_PORTE:
        sierra_path = SIERRA_PATH / 'port-S_10MHz_3GHz.csv'
        labels = [str(port) for port in PORT_WIDTHS]
        for i, port in enumerate(PORT_WIDTHS):
            path = base_path / f'coplanar_waveguide_lumped_22um_122.64ohm_cond_59600000_porte_{port}um_length_1mm/port-S.csv'
            data = read_csv_data(path, [1, 3])
            s11_all.append(data[0])
            s21_all.append(data[1])
        # Calculate skin effect threshold
        skin_effect_thr = calculate_skin_effect_threshold(SKIN_EFF_FRACTION, SIMULATION_CONDUCTIVITY, 1000e-6, 90e-6, 20e-6)


    if base_path == BASE_PATH_NUOVI_TEST_PEC:
        sierra_path = SIERRA_PATH / 'port-S_10MHz_3GHz.csv'
        labels = ['lumped_NOPEC_Absorbing_4',
                  'lumped_PEC_9-11_Absorbing_4',
                  'lumped_PEC_9_Absorbing_4',
                  'lumped_PEC_9-11_NOAbsorbing',
                  'lumped_NOPEC_Absorbing_4_105ohm',
                  'wave_PEC_9-11_Absorbing_8',
                  #'wave_PEC_9_Absorbing_6-7-8',
                  'wave_PEC_9-11_NOAbsorbing',
                  'wave_PEC_6-7-8-9-11_NOAbsorbing',
                  #'wave_PEC_6-7-9-11_Absorbing_8',
                  'wave_PEC_9-11_Absorbing_8_PMC_12',
                  'wave_PEC_9-11_NOAbsorbing_PMC_12',
                  ]
        paths = ['cpw_lumped_22um_122.64ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 'cpw_lumped_22um_122.64ohm_cond_59600000_porte_75um_PEC_9-11_Absorbing_4',
                 'cpw_lumped_22um_122.64ohm_cond_59600000_porte_75um_PEC_9_Absorbing_4',
                 'cpw_lumped_22um_122.64ohm_cond_59600000_porte_75um_PEC_9-11_NOAbsorbing',
                 'cpw_lumped_22um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 'cpw_wave_22um_porte720x450_cond_59600000_PEC_9-11_Absorbing_8',
                 #'cpw_wave_22um_porte720x450_cond_59600000_PEC_9_Absorbing_6-7-8',
                 'cpw_wave_22um_porte720x450_cond_59600000_PEC_9-11_NOAbsorbing',
                 'cpw_wave_22um_porte720x450_cond_59600000_PEC_6-7-8-9-11_NOAbsorbing',
                 #'cpw_wave_22um_porte720x450_cond_59600000_PEC_6-7-9-11_Absorbing_8',
                 'cpw_wave_22um_porte720x450_cond_59600000_PEC_9-11_Absorbing_8_PMC_12',
                 'cpw_wave_22um_porte720x450_cond_59600000_PEC_9-11_NOAbsorbing_PMC_12',
                 ]
        for i, path in enumerate(paths):
            path = base_path / f'{path}/port-S.csv'
            data = read_csv_data(path, [0, 1, 3])
            x_all.append(data[0])
            s11_all.append(data[1])
            s21_all.append(data[2])
        # Calculate skin effect threshold
        skin_effect_thr = calculate_skin_effect_threshold(SKIN_EFF_FRACTION, SIMULATION_CONDUCTIVITY, 1000e-6, 90e-6, 20e-6)


    if base_path == BASE_PATH_REPORT:
        sierra_path = SIERRA_PATH / 'port-S_50MHz_3GHz.csv'

        # all simulation here are LUMPED
        
        
        #Diverse Impedenze
        # all simulations here are LUMPED, NOPEC, Absorbing_4
        labels = ['[$Z_{S,L}=122\\Omega$]',
                  '[$Z_{S,L}=104\\Omega$]',
                  '[$Z_{S,L}=105\\Omega$]',
                  #'[$Z_{S,L}=106\\Omega$]',
                  #'[$Z_{S,L}=107\\Omega$]',
                  #'[$Z_{S,L}=108\\Omega$]',
                  #'[$Z_{S,L}=109\\Omega$]',
                  #'[$Z_{S,L}=110\\Omega$]'
                  ]
        paths = ['cpw_lumped_22um_122ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 'cpw_lumped_22um_104ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 'cpw_lumped_22um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 #'cpw_lumped_22um_106ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 #'cpw_lumped_22um_107ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 #'cpw_lumped_22um_108ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 #'cpw_lumped_22um_109ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 #'cpw_lumped_22um_110ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4'
                 ]
        for i, path in enumerate(paths):
            path = base_path / f'DiverseImpedenze/DiverseImpedenze50MHz/{path}/port-S.csv'
            data = read_csv_data(path, [0, 1, 2, 3, 4])
            x_all.append(data[0])
            s11_all.append(data[1])
            s11_phases_all.append(data[2])
            s21_all.append(data[3])
            s21_phases_all.append(data[4])
        
        '''
        #Diverse Mesh
        labels = ['[mesh size=15.0$\\mu$m]',
                  '[mesh size=22.5$\\mu$m]',
                  '[mesh size=45.0$\\mu$m]',
                  '[mesh size=67.5$\\mu$m]'
                  ]
        paths = ['cpw_lumped_15um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 'cpw_lumped_22um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 'cpw_lumped_45um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 'cpw_lumped_67um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4'
                ]
        for i, path in enumerate(paths):
            path = base_path / f'DiverseMesh/DiverseMesh50MHz/{path}/port-S.csv'
            data = read_csv_data(path, [0, 1, 2, 3, 4])
            x_all.append(data[0])
            s11_all.append(data[1])
            s11_phases_all.append(data[2])
            s21_all.append(data[3])
            s21_phases_all.append(data[4])
        '''
        '''
        #Diverse Conducibilita
        labels = ['[$\\sigma=58.0MS/m$]',
                  #'[$\\sigma=59.0MS/m$]',
                  '[$\\sigma=59.6MS/m$]',
                  '[$\\sigma=62.5MS/m$]'
                ]
        paths = ['cpw_lumped_22um_105ohm_cond_58000000_porte_75um_NOPEC_Absorbing_4',
                 #'cpw_lumped_22.5um_105ohm_cond_59000000_porte_75um_NOPEC_Absorbing_4',
                 'cpw_lumped_22um_105ohm_cond_59600000_porte_75um_NOPEC_Absorbing_4',
                 'cpw_lumped_22um_105ohm_cond_62500000_porte_75um_NOPEC_Absorbing_4'
                ]
        for i, path in enumerate(paths):
            path = base_path / f'DiverseConducibilita/DiverseConducibilita50MHz/{path}/port-S.csv'
            data = read_csv_data(path, [0, 1, 2, 3, 4])
            x_all.append(data[0])
            s11_all.append(data[1])
            s11_phases_all.append(data[2])
            s21_all.append(data[3])
            s21_phases_all.append(data[4])
        '''

        # Calculate skin effect threshold
        skin_effect_thr = calculate_skin_effect_threshold(SKIN_EFF_FRACTION, SIMULATION_CONDUCTIVITY, 1000e-6, 90e-6, 20e-6)


    # Sierra data
    sierra_x, sierra_s11, sierra_s21 = read_csv_data(sierra_path, [0, 1, 3])

    dielectric_losses = []
    conductor_losses = []
    print("\n===========Calculating dielectric losses===========")
    dielectric_losses.extend(calculate_dielectric_losses(i) for i in x_all[0])
    print("\n\n===========Calculating conductor losses===========")
    conductor_losses.extend(calculate_conductor_losses(i) for i in x_all[0])
    
    s11_err, s21_err = estimate_systematic_uncertainty()
    calculate_RGLC_and_plot_Z0_from_the_105Ohm_simulation(False)

    compare_with_Sierra = True
    chi2vsSierra = []

    ####################################S11####################################
    if compare_with_Sierra:
        print("\n===========Chi-square S11 vs Sierra===========")
        for i, myData in enumerate(s11_all):
            chi2 = chisq(np.array(myData), np.array(sierra_s11), s11_err)
            chi2vsSierra.append(chi2)
            print(f"S11 Chi-squared with {len(myData) - 2} dof for {labels[i]} vs Sierra: {chi2}")

    # Plot S11
    plot_S11(x_all + [sierra_x], s11_all + [sierra_s11], s11_err, skin_effect_thr, chi2vsSierra,
            labels + ['Sierra 122.64 Ohm'], COLORS + ['r'], 
            True, True, compare_with_Sierra, False, False)


    ####################################S21####################################
    chi2vsSierra.clear()
    chi2vsIPC = []
    if compare_with_Sierra:
        print("\n===========Chi-square S21 vs Sierra===========")
        for i, myData in enumerate(s21_all):
            chi2 = chisq(np.array(myData), np.array(sierra_s21), s21_err)
            chi2vsSierra.append(chi2)
            print(f"S21 Chi-squared with {len(myData) - 1} dof for {labels[i]} vs Sierra: {chi2}")

    print("\n===========Chi-square S21 vs IPC===========")
    print("Without applying cuts on frequency range")
    for i, myData in enumerate(s21_all):
        chi2 = chisq(np.array(myData), -np.array(np.add(dielectric_losses, conductor_losses)), s21_err)
        print(f"S21 Chi-squared with {len(myData) - 1} dof for {labels[i]} vs IPC: {chi2}")

    print("Cutting frequencies over skin effect threshold")
    for i, myData in enumerate(s21_all):
        freq_index_after_thr = len(myData)
        if skin_effect_thr < 2.9:
            freq_index_after_thr = next(j for j, freq in enumerate(x_all[i]) if freq > skin_effect_thr)
        myData = myData[:freq_index_after_thr]
        chi2 = chisq(np.array(myData), -np.array(np.add(dielectric_losses[:freq_index_after_thr], conductor_losses[:freq_index_after_thr])), s21_err[:freq_index_after_thr])
        print(f"S21 Chi-squared with {len(myData) - 1} dof for {labels[i]} vs IPC: {chi2}")

    print("\n===========Chi-square S21 vs IPC with errors on IPC value===========")
    print("Cutting frequencies over skin effect threshold")
    loss_err_2_percent = np.array(np.multiply(np.add(dielectric_losses, conductor_losses), 0.02))
    loss_err_3_percent = np.array(np.multiply(np.add(dielectric_losses, conductor_losses), 0.03))
    loss_err_4_percent = np.array(np.multiply(np.add(dielectric_losses, conductor_losses), 0.04))
    for i, myData in enumerate(s21_all):
        freq_index_after_thr = len(myData)
        if skin_effect_thr < 2.9: 
            freq_index_after_thr = next(j for j, freq in enumerate(x_all[i]) if freq > skin_effect_thr)
        myData = myData[:freq_index_after_thr]
        chi2 = chisq(-np.array(np.add(dielectric_losses[:freq_index_after_thr], conductor_losses[:freq_index_after_thr])), np.array(myData), loss_err_2_percent[:freq_index_after_thr])
        print(f"S21 Chi-squared with {len(myData) - 2} dof, assuming dispersion from IPC of 2% for {labels[i]} vs IPC: {chi2}")
    for i, myData in enumerate(s21_all):
        freq_index_after_thr = len(myData)
        if skin_effect_thr < 2.9: 
            freq_index_after_thr = next(j for j, freq in enumerate(x_all[i]) if freq > skin_effect_thr)
        myData = myData[:freq_index_after_thr]
        chi2 = chisq(-np.array(np.add(dielectric_losses[:freq_index_after_thr], conductor_losses[:freq_index_after_thr])), np.array(myData), loss_err_3_percent[:freq_index_after_thr])
        chi2vsIPC.append(chi2) #per ora salvo solo quello al 3% perchè sembra il risultato migliore
        dof = len(myData) - 2
        print(f"S21 Chi-squared with {len(myData) - 2} dof, assuming dispersion from IPC of 3% for {labels[i]} vs IPC: {chi2}")
    for i, myData in enumerate(s21_all):
        freq_index_after_thr = len(myData)
        if skin_effect_thr < 2.9:
            freq_index_after_thr = next(j for j, freq in enumerate(x_all[i]) if freq > skin_effect_thr)
        myData = myData[:freq_index_after_thr]
        chi2 = chisq(-np.array(np.add(dielectric_losses[:freq_index_after_thr], conductor_losses[:freq_index_after_thr])), np.array(myData), loss_err_4_percent[:freq_index_after_thr])
        print(f"S21 Chi-squared with {len(myData) - 2} dof, assuming dispersion from IPC of 4% for {labels[i]} vs IPC: {chi2}")

    plot_S21(x_all + [sierra_x], s21_all + [sierra_s21], s21_err, dielectric_losses, conductor_losses, skin_effect_thr, chi2vsSierra, chi2vsIPC,
            labels + ['Sierra 122.64 Ohm'], COLORS + ['r'], dof,
            True, True, True, compare_with_Sierra, False, False)

    plt.show()


def S_parameters_vs_length():
    x_all = []
    s11_all = []
    s21_all = []
    base_path = Path('postpro/Alluminio/Lumped/CPW/DiverseLunghezze/500MHz')
    labels = ['L=1mm',
              'L=2mm',
              'L=3mm',
              'L=4mm',
              'L=5mm',
              'L=6mm',
              'L=7mm',
              'L=8mm',
              'L=9mm',
              'L=10mm',
              ]
    paths = ['cpw_lumped_22um_105ohm_cond_33112582_port0.075_1000um_NOPEC_Absorbing_4',
             'cpw_lumped_22um_105ohm_cond_33112582_port0.075_2000um_NOPEC_Absorbing_4',
             'cpw_lumped_22um_105ohm_cond_33112582_port0.075_3000um_NOPEC_Absorbing_4',
             'cpw_lumped_22um_105ohm_cond_33112582_port0.075_4000um_NOPEC_Absorbing_4',
             'cpw_lumped_22um_105ohm_cond_33112582_port0.075_5000um_NOPEC_Absorbing_4',
             'cpw_lumped_22um_105ohm_cond_33112582_port0.075_6000um_NOPEC_Absorbing_4',
             'cpw_lumped_22um_105ohm_cond_33112582_port0.075_7000um_NOPEC_Absorbing_4',
             'cpw_lumped_22um_105ohm_cond_33112582_port0.075_8000um_NOPEC_Absorbing_4',
             'cpw_lumped_22um_105ohm_cond_33112582_port0.075_9000um_NOPEC_Absorbing_4',
             'cpw_lumped_22um_105ohm_cond_33112582_port0.075_10000um_NOPEC_Absorbing_4',
             ]
    x_all = [1, 
             2, 
             3,
             4,
             5,
             6,
             7,
             8,
             9,
             10
             ]
    for i, path in enumerate(paths):
        path = base_path / f'{path}/port-S.csv'
        data = read_csv_data(path, [1, 3])
        s11_all.extend(data[0])
        s21_all.extend(data[1])
    # Calculate skin effect threshold
    skin_effect_thr = calculate_skin_effect_threshold(SKIN_EFF_FRACTION, SIMULATION_CONDUCTIVITY, 1000e-6, 90e-6, 20e-6)

    fig1, ax1 = plt.subplots()
    ax1.errorbar(x_all, s11_all, xerr=0, yerr=0, color='red', linestyle='solid', marker='o', markersize=MARKERSIZE)
    ax1.set_title('S11 vs Line Length @0.5GHz', fontsize=18)
    ax1.grid(True)
    #ax1.legend(fontsize='x-large')
    ax1.set_xlabel('Line Length [mm]')
    ax1.set_ylabel('|S11|')

    fig2, ax2 = plt.subplots()
    ax2.errorbar(x_all, s21_all, xerr=0, yerr=0, color='red', linestyle='solid', marker='o', markersize=MARKERSIZE)
    ax2.set_title('S21 vs Line Length @0.5GHz', fontsize=18)
    ax2.grid(True)
    #ax2.legend(fontsize='x-large')
    ax2.set_xlabel('Line Length [mm]')
    ax2.set_ylabel('|S21|')

    '''
    Aggiungere queste simulazioni a frequenza più alta per vedere oscillazioni di S11?
    '''
    '''
    s11_all = []
    s21_all = []
    base_path = Path('postpro/Alluminio/DiverseLunghezze/3GHz')
    for i, path in enumerate(paths):
        path = base_path / f'{path}/port-S.csv'
        data = read_csv_data(path, [1, 3])
        s11_all.append(data[0])
        s21_all.append(data[1])
    # Calculate skin effect threshold
    skin_effect_thr = calculate_skin_effect_threshold(SKIN_EFF_FRACTION, SIMULATION_CONDUCTIVITY, 1000e-6, 90e-6, 20e-6)

    fig3, ax3 = plt.subplots()
    for x, y, label, color in zip(x_all, s11_all, labels, COLORS):
        ax3.errorbar(x, y, xerr=0, yerr=0, color='red', linestyle='solid', marker='o', markersize=MARKERSIZE)
    ax3.set_title('S11 vs Line Length @3GHz', fontsize=18)
    ax3.grid(True)
    #ax3.legend(fontsize='x-large')
    ax3.set_xlabel('Line Length [mm]')
    ax3.set_ylabel('|S11|')

    fig4, ax4 = plt.subplots()
    for x, y, label, color in zip(x_all, s21_all, labels, COLORS):
        ax4.errorbar(x, y, xerr=0, yerr=0, color='red', linestyle='solid', marker='o', markersize=MARKERSIZE)
    ax4.set_title('S21 vs Line Length @3GHz', fontsize=18)
    ax4.grid(True)
    #ax4.legend(fontsize='x-large')
    ax4.set_xlabel('Line Length [mm]')
    ax4.set_ylabel('|S21|')
    '''

    plt.xticks(rotation=25)
    #plt.tight_layout()
    plt.show()


def plot_data_vs_frequency(path='', label='', one_port_only=False):
    
    if path == '':
        return
    
    base_path = Path('postpro/Alluminio/Lumped/CPW/DiverseLunghezze/500MHz')

    full_path = f'{path}/port-S.csv'
    if one_port_only:
        data = read_csv_data(full_path, [0, 1])
    else:
        data = read_csv_data(full_path, [0, 1, 3])
    
    # Calculate skin effect threshold
    skin_effect_thr = calculate_skin_effect_threshold(SKIN_EFF_FRACTION, SIMULATION_CONDUCTIVITY, 1000e-6, 90e-6, 20e-6)

    fig1, ax1 = plt.subplots()
    ax1.errorbar(data[0], data[1], xerr=0, yerr=0, color='red', linestyle='solid', marker='o', markersize=MARKERSIZE, label=label)
    ax1.set_title('S11 vs Frequency', fontsize=18)
    ax1.grid(True)
    #ax1.legend(fontsize='x-large')
    ax1.set_xlabel('Frequency [GHz]')
    ax1.set_ylabel('|S11|')

    if not one_port_only:
        fig2, ax2 = plt.subplots()
        ax2.errorbar(data[0], data[2], xerr=0, yerr=0, color='blue', linestyle='solid', marker='o', markersize=MARKERSIZE, label=label)
        ax2.set_title('S21 vs Frequency', fontsize=18)
        ax2.grid(True)
        #ax2.legend(fontsize='x-large')
        ax2.set_xlabel('Frequency [GHz]')
        ax2.set_ylabel('|S21|')

    plt.xticks(rotation=25)
    #plt.tight_layout()
    plt.show()




# --- MAIN EXECUTION ---
if __name__ == '__main__':

    #plot_multiple_simulations(BASE_PATH_102OHM)
    #plot_multiple_simulations(BASE_PATH_102OHM_DIVERSE_CONDUCIBILITA)
    #plot_multiple_simulations(BASE_PATH_TEST)
    #plot_multiple_simulations(BASE_PATH_DIVERSE_PORTE)
    #plot_multiple_simulations(BASE_PATH_NUOVI_TEST_PEC)

    if SIMULATION_CONDUCTIVITY == CONDUCTIVITY_COPPER and CPW_LENGTH == 1000e-6: plot_multiple_simulations(BASE_PATH_REPORT)

    ################################ TEST ALTRE FUNZIONI ################################

    #calculate_free_space_radiation_losses()
    #plot_Z0_vs_software()
    #calculate_Z0_from_conductor_losses()
    #calculate_free_space_radiation_losses()

    #calculate_RGLC_and_plot_Z0_from_the_105Ohm_simulation()
    
    if SIMULATION_CONDUCTIVITY == CONDUCTIVITY_ALUMINIUM and CPW_LENGTH == 10000e-6:
        plot_single_simulation('postpro/Alluminio/Lumped/CPW/TestBeole_1cm/cpw_lumped_22um_105ohm_cond_33112582_porte_75um_NOPEC_Absorbing_4/port-S.csv',
                            calculate_skin_effect_threshold(SKIN_EFF_FRACTION, SIMULATION_CONDUCTIVITY, CPW_LENGTH, 90e-6, 20e-6),
                            '[FGCPW-Aluminum]', sierra_path=SIERRA_PATH / 'port-S_50MHz_3GHz.csv')

    if SIMULATION_CONDUCTIVITY == CONDUCTIVITY_COPPER and CPW_LENGTH == 1000e-6:
        plot_single_simulation('postpro/Rame/TestPEC_and_small_ports/cpw_lumped_22um_105ohm_cond_59600000_port0.0075_1000um_PEC_9_Absorbing_4/port-S.csv',
                                calculate_skin_effect_threshold(SKIN_EFF_FRACTION, SIMULATION_CONDUCTIVITY, CPW_LENGTH, 90e-6, 20e-6),
                                '[FGCPW-Copper_PEC_9_and_port_0.75%]', sierra_path=SIERRA_PATH / 'port-S_50MHz_3GHz.csv')
        
    if SIMULATION_CONDUCTIVITY == CONDUCTIVITY_ALUMINIUM and CPW_LENGTH == 20e-6:
        plot_single_simulation('postpro/Alluminio/Lumped/CPW/TestRefinedMesh/cpw_lumped_2um_105ohm_cond_32894736_port0.0075_20um_NOPEC_Absorbing_4/port-S.csv',
                                calculate_skin_effect_threshold(SKIN_EFF_FRACTION, SIMULATION_CONDUCTIVITY, CPW_LENGTH, 90e-6, 20e-6),
                                '[FGCPW-Aluminum]', sierra_path=SIERRA_PATH / 'port-S_50MHz_3GHz.csv')
        
    if SIMULATION_CONDUCTIVITY == CONDUCTIVITY_ALUMINIUM and CPW_LENGTH == 100e-6:
        plot_single_simulation('postpro/Alluminio/Lumped/CPW/TestRefinedMesh/cpw_lumped_2um_105ohm_cond_32894736_port0.075_100um_NOPEC_Absorbing_4/port-S.csv',
                                calculate_skin_effect_threshold(SKIN_EFF_FRACTION, SIMULATION_CONDUCTIVITY, CPW_LENGTH, 90e-6, 20e-6),
                                '[FGCPW-Aluminum]', sierra_path=SIERRA_PATH / 'port-S_50MHz_3GHz.csv')
    

    """
    Analisi della CPW short del paper "Study of the Radio Frequency (RF) performance 
    of a Wafer-Level Package (WLP) with Through Silicon Vias (TSVs) for the integration 
    of RF-MEMS and micromachined waveguides in the context of 5G and Internet of Things (IoT) applications: 
    Part 1—validation of the 3D modelling approach"
    (design A)
    """
    calculate_Z0_from_conductor_losses(trace_width=116e-6, trace_thickness=20e-6, gap_width=65e-6, substrate_thickness=525e-6, dielectric_constant=11.7)
    plot_data_vs_frequency('postpro/AltriMateriali/Lumped/ShortCPW/cpw_lumped_22um_40ohm_cond_58800000_port0.075_1350um_NOPEC_Absorbing_4',
                           '[Cu-HRS]', #HRS = High Resistivity Silicon @dielectric constant=11.7 and loss tangent=0.0013 (?)
                           one_port_only=True)
    
    
    
    #S_parameters_vs_length()
    #calculate_skin_effect_threshold(0, CONDUCTIVITY_COPPER, 1000e-6, 90e-6,20e-6, Z0)
    #calculate_skin_effect_threshold(0, CONDUCTIVITY_ALUMINIUM, 1000e-6, 90e-6,20e-6, Z0)

    
    calculate_Z0_from_conductor_losses(trace_width=90e-6*0.7, trace_thickness=20e-6, gap_width=100e-6+0.3*90e-6, substrate_thickness=25e-6, dielectric_constant=3.3)
    
    
    
    

