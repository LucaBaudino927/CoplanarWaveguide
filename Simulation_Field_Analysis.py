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
BASE_PATH_REPORT = Path('postpro/Rame/Lumped/CPW/Report')

COLORS = ['green', 'red', 'blue', 'orange', 'cyan', 'black', 'magenta', 'purple', 'brown', 'pink', 'gray']
MARKERSIZE = 5

# --- FUNCTION DEFINITIONS ---


def read_csv_field_data(path, probe_number=1):
    data = []
    column_indeces = [0, 
                      6*(probe_number - 1) + 1, 
                      6*(probe_number - 1) + 2, 
                      6*(probe_number - 1) + 3, 
                      6*(probe_number - 1) + 4, 
                      6*(probe_number - 1) + 5, 
                      6*(probe_number - 1) + 6]
    with open(path, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            data.append([float(row[i]) for i in column_indeces])
    """
    Ex = [np.sqrt(data[i][1]**2 + data[i][2]**2) for i in range(len(data))]
    print(f'data[0][1]: {data[0][1]}, data[0][2]: {data[0][2]}, Ex[0]: {Ex[0]}')
    Ey = [np.sqrt(data[i][3]**2 + data[i][4]**2) for i in range(len(data))]
    print(f'data[0][3]: {data[0][3]}, data[0][4]: {data[0][4]}, Ey[0]: {Ey[0]}')
    Ez = [np.sqrt(data[i][5]**2 + data[i][6]**2) for i in range(len(data))]
    print(f'data[0][5]: {data[0][5]}, data[0][6]: {data[0][6]}, Ez[0]: {Ez[0]}')
    E = [np.sqrt(Ex[i]**2 + Ey[i]**2 + Ez[i]**2) for i in range(len(data))]
    print(f'E[0]: {E[0]}')
    """
    Ey = [np.sqrt(data[i][3]**2 + data[i][4]**2) for i in range(len(data))]
    Ez = [np.sqrt(data[i][5]**2 + data[i][6]**2) for i in range(len(data))]
    E = [np.sqrt(Ey[i]**2 + Ez[i]**2) for i in range(len(data))]
    return list(zip([d[0] for d in data], E))  # Returns list of columns


def plot_data_vs_frequency(path='', label='', probe_number=251):
    
    if path == '':
        return

    y_probe = list(range(250, 500 +1))  # y positions of probes in um
    z_probe = np.full(probe_number, 6.0)      # z positions of probes in um
    for j in range(len(y_probe)):
        y_probe[j] = y_probe[j] - 378.0 # center the probe positions around zero
    r_probe = np.sqrt(np.array(y_probe)**2 + np.array(z_probe)**2)
    for j in range(len(y_probe)):
        if y_probe[j] < 0:
            r_probe[j] = -r_probe[j]
    #print(f'Probe radial positions (um): {r_probe}')

    E_fixed_freq = []
    full_path = f'{path}/probe-E.csv'
    for i in range(1, probe_number + 1):
        data = read_csv_field_data(full_path, i)
        marker_size = MARKERSIZE
        if len(data) > 100:
            marker_size = 1
        
        if (i) % 10 == 0: #plot every 10 probes
            fig1, ax1 = plt.subplots()
            ax1.errorbar([d[0] for d in data], [d[1] for d in data], xerr=0, yerr=0, color=COLORS[i % len(COLORS)], linestyle='solid', marker='o', markersize=marker_size, label=f'Probe {i} {label}')
            ax1.set_title(f'E field at Probe {i}, y = {y_probe[i-1]}', fontsize=18)
            ax1.grid(True)
            ax1.legend(fontsize='x-large')
            ax1.set_xlabel('Frequency [GHz]')
            ax1.set_ylabel('|E| [V/m]')

        ## frequency range: i = 0 -> 0.05 GHz, i = 296 -> 3.0 GHz
        freq_index = 0
        fixed_frequency = data[freq_index][0] # data is a vector of (frequency, E) tuples, data[0][0] takes the lowest frequency
        E_fixed_freq.append(data[freq_index][1]) # data is a vector of (frequency, E) tuples, data[0][1] takes E at first frequency
        print(f'Probe {i}, y = {y_probe[i-1]} um, f = {fixed_frequency} GHz, |E| = {data[freq_index][1]} V/m')

    fig2, ax2 = plt.subplots()
    ax2.errorbar(r_probe, E_fixed_freq, xerr=0, yerr=0, color=COLORS[i % len(COLORS)], linestyle='solid', marker='o', markersize=marker_size, label=f'Probe {i} {label}')
    ax2.set_title(f'E field vs probe position, f = {fixed_frequency} GHz', fontsize=18)
    ax2.grid(True)
    ax2.legend(fontsize='x-large')
    ax2.set_xlabel('y [um]')
    ax2.set_ylabel('|E| [V/m]')

    width = 50.0 #um
    gap = 23.0 #um
    a = width / 2
    b = a + gap
    k = width / (width + 2 * gap)
    np.seterr(divide = 'ignore')
    ax2.plot(r_probe, 100000000 / (K(k) * np.sqrt(np.abs(r_probe**2 - a**2)) * np.sqrt(np.abs(r_probe**2 - b**2))), 'g--')

    plt.xticks(rotation=25)
    #plt.tight_layout()
    plt.show()


# --- MAIN EXECUTION ---
if __name__ == '__main__':

    path = 'postpro/Gold/Lumped/CalKit/LOAD103R/TWOPORT_cpw_lumped_3.33um_50_50ohm_cond_45000000_port0.075_50um_NOPEC_Absorbing_4'
    #'postpro/Gold/Lumped/CalKit/0509/TWOPORT_cpw_lumped_10um_50_99999999ohm_cond_45000000_port0.075_1800um_NOPEC_Absorbing_4'
    #'postpro/Gold/Lumped/CalKit/0509/TWOPORT_cpw_lumped_8um_50_50ohm_cond_45000000_port0.075_1600um_NOPEC_Absorbing_4'
    #'postpro/Gold/Wave/CalKit/0509/ONEPORT_cpw_wave_10um_removeMetal_condSubstrate_port500x785_1800um_PEC_5_6_8_10_Absorbing_7'
    #'postpro/Gold/Wave/CalKit/0509/ONEPORT_cpw_wave_10um_cond_45000000_port500x785_1800um_PEC_5_6_8_10_Absorbing_7'
    #'postpro/Gold/Wave/CalKit/0509/ONEPORT_cpw_wave_10um_removeMetal_port500x785_1800um_PEC_7_9_NOAbsorbing'
    #'postpro/Gold/Wave/CalKit/0509/ONEPORT_cpw_wave_10um_removeMetal_port500x785_1800um_PEC_4_5_7_9_Absorbing_6'
    #'postpro/Gold/Lumped/CalKit/0509/TWOPORT_cpw_lumped_10um_50_99999999ohm_cond_45000000_port0.075_1800um_NOPEC_Absorbing_4'
    #'postpro/Gold/Lumped/CalKit/0509/TWOPORT_cpw_lumped_10um_50_378ohm_cond_45000000_port0.075_1800um_NOPEC_Absorbing_4'
    label = '[CalKitLOAD103R]'
    #'[CalKit0509-Open]'
    plot_data_vs_frequency(path, label, 251)


    #========================================================================================================================================#
    

