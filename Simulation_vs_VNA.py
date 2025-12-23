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

IMPEDANCES = [100, 102, 102.1, 102.5]
CONDUCTIVITIES = [33112582, 57471264, 59600000, 62500000, 72500000] #S/m
PORT_WIDTHS = [10, 50, 75, 100, 150, 200, 490]  #um
COLORS = ['green', 'red', 'blue', 'orange', 'cyan', 'black', 'magenta', 'purple', 'brown', 'pink', 'gray']
MARKERSIZE = 5
SIERRA_PATH = Path('postpro/Rame/Lumped/CPW/SierraData_122ohm')

CONDUCTIVITY_COPPER = 59600000  # S/m
CONDUCTIVITY_ALUMINIUM = 32894736.84 #S/m = 1/3.04x10-8  #33112582.78 #S/m 
SKIN_EFF_FRACTION = 0.5  # Fraction increase in conductor losses to consider skin effect relevant
Z0_STANDARD_CPW = 105  # Characteristic impedance in Ohm of my 90-100-20 CPW
Z0_CALKIT = 50
CPW_LENGTH = 1000e-6
SIMULATION_CONDUCTIVITY = CONDUCTIVITY_COPPER
SIMULATION_Z0 = Z0_CALKIT

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


def calculate_skin_effect_threshold(conductivity=SIMULATION_CONDUCTIVITY, trace_length=CPW_LENGTH, trace_width=90e-6, trace_thickness=20e-6, Z0=SIMULATION_Z0):
    
    """
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


def calculate_conductor_losses(frequency=3, conductivity=SIMULATION_CONDUCTIVITY, trace_width=90e-6, trace_thickness=20e-6, trace_length=CPW_LENGTH, gap_width=100e-6, Z0=SIMULATION_Z0):
    """
    Fraction of power loss in the conductor for a CPW line.
    It calculates Ac in dB/m and return the value in dB.
    """
    mu0 = 4e-7 * np.pi
    delta = np.sqrt(1 / (np.pi * frequency*1e9 * mu0 * conductivity))  # Skin depth

    Rs = (1 + 1j) / (conductivity * delta) # Surface resistance in Ohm
    if delta > trace_thickness: Rs = (1 + 1j) / (conductivity * trace_thickness)
    Rstrip = Rs.real / (trace_width) # Ohm/m
    # Conductor loss in dB/m
    Ac = (8.686 * Rstrip) / (2 * Z0)

    
    '''
    Test dal paper: "A new analytical, cad-oriented model for the ohmic and radiation loesses of asymmetric coplanar waveguides in
    hybrid and monolithic mic's", G. Ghione, C.U. Naldi
    '''
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

    InverseIntegral = (240 * np.pi * K(ks) * K(ksprime) * (1 - ks*ks)) / (1/a * (np.pi + np.log(8 * np.pi * a * (1 - ks) / (trace_thickness * (1 + ks)))) + 
                                                                                                1/b * (np.pi + np.log(8 * np.pi * b * (1 - ks) / (trace_thickness * (1 + ks)))))
    Z0_calc = InverseIntegral / trace_width / np.sqrt(3.3)
    #print(f'Z0: {Z0_calc}')
    Ac = 8.686 * Rstrip / (2 * Z0_calc)

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


def calculate_Z_and_propagation_constant_from_S_matrix(s11_all, s11_phases_all, s21_all, s21_phases_all, Z0=SIMULATION_Z0, trace_length=CPW_LENGTH):
    
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
        #print(f'20*np.log(np.abs((Z0 - y)/(Z0 + y)))={20*np.log(np.abs((Z0 - y)/(Z0 + y)))}, Z0={Z0}, y={y}, x={x}')
        expected_s11.append( 20*np.log(1/(1 - (np.abs((SIMULATION_Z0 - y)/(SIMULATION_Z0 + y)))*(np.abs((SIMULATION_Z0 - y)/(SIMULATION_Z0 + y))) )) )
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
    Z, gamma = calculate_Z_and_propagation_constant_from_S_matrix(s11, s11_phase, s21, s21_phase, SIMULATION_Z0, trace_length=CPW_LENGTH)
    Z0_calc = calculate_RGLC_and_plot_Z0(Z, gamma, x, labels[0])
    if show_plot: plt.show()
    return Z0_calc


def plot_data_vs_frequency(path='', label='', one_port_only=False):
    
    if path == '':
        return

    full_path = f'{path}/port-S.csv'
    if one_port_only:
        data = read_csv_data(full_path, [0, 1, 2])
    else:
        data = read_csv_data(full_path, [0, 1, 2, 3, 4])

    marker_size = MARKERSIZE
    if len(data[0]) > 100:
        marker_size = 1
    
    # Calculate skin effect threshold
    skin_effect_thr = calculate_skin_effect_threshold(SIMULATION_CONDUCTIVITY, 1000e-6, 90e-6, 20e-6)

    fig1, ax1 = plt.subplots()
    ax1.errorbar(data[0], data[1], xerr=0, yerr=0, color='red', linestyle='solid', marker='o', markersize=marker_size, label=label)
    ax1.set_title('S11 vs Frequency', fontsize=18)
    ax1.grid(True)
    ax1.legend(fontsize='x-large')
    ax1.set_xlabel('Frequency [GHz]')
    ax1.set_ylabel('|S11|')

    if not one_port_only:
        fig2, ax2 = plt.subplots()
        ax2.errorbar(data[0], data[3], xerr=0, yerr=0, color='blue', linestyle='solid', marker='o', markersize=marker_size, label=label)
        ax2.set_title('S21 vs Frequency', fontsize=18)
        ax2.grid(True)
        ax2.legend(fontsize='x-large')
        ax2.set_xlabel('Frequency [GHz]')
        ax2.set_ylabel('|S21| [dB]')

    if one_port_only:
        Z_modulo = []
        Z_phase = []
        for i_freq in enumerate(data[0]):
            i = i_freq[0]
            S11_modulo = 10**(data[1][i]/20)
            S11_phase = np.pi/180*data[2][i]
            S11_real = S11_modulo*np.cos(S11_phase)
            S11_imag = S11_modulo*np.sin(S11_phase)
            Z = SIMULATION_Z0 * (1 + S11_real + 1j*S11_imag) / (1 - S11_real - 1j*S11_imag)
            Z_modulo.append(np.sqrt(Z.real**2 + Z.imag**2))
            Z_phase.append(np.atan(Z.imag/Z.real))

        fig3, ax3 = plt.subplots()
        ax3.errorbar(data[0], Z_modulo, xerr=0, yerr=0, color='blue', linestyle='solid', marker='o', markersize=marker_size, label=label)
        ax3.set_title('|Z| vs Frequency', fontsize=18)
        ax3.grid(True)
        ax3.legend(fontsize='x-large')
        ax3.set_xlabel('Frequency [GHz]')
        ax3.set_ylabel('$|Z| [\\Omega]$')

        fig4, ax4 = plt.subplots()
        ax4.errorbar(data[0], Z_phase, xerr=0, yerr=0, color='blue', linestyle='solid', marker='o', markersize=marker_size, label=label)
        ax4.set_title('arg(Z) vs Frequency', fontsize=18)
        ax4.grid(True)
        ax4.legend(fontsize='x-large')
        ax4.set_xlabel('Frequency [GHz]')
        ax4.set_ylabel('arg(Z) [rad]')

    plt.xticks(rotation=25)
    #plt.tight_layout()
    plt.show()

def plot_sim_vs_VNA(sim_path = '', VNA_path = '', label = '', one_port_only=False):
    
    if sim_path == '' and VNA_path == '':
        return

    sim_data = []
    if sim_path != '':
        if one_port_only:
            sim_data.append(read_csv_data(f'{sim_path}/port-S.csv', [0, 1, 2])) # Freq, S11 [dB], arg(S11) [deg]
        else:
            sim_data.append(read_csv_data(f'{sim_path}/port-S.csv', [0, 1, 2, 3, 4])) # Freq, S11 [dB], arg(S11) [deg], S21 [dB], arg(S21) [deg]

    VNA_data = []
    if one_port_only:
        VNA_data.append(read_csv_data(f'{VNA_path}/TRACE01.CSV', [0, 1]))
    """
    I comment this part because we don't have S21 data from VNA yet
    else:
        VNA_data.append(read_csv_data(f'{VNA_path}/TRACE01.CSV', [0, 1, 3]))
    """
    
    # Calculate skin effect threshold
    #skin_effect_thr = calculate_skin_effect_threshold(SIMULATION_CONDUCTIVITY, SIMULATION_LINE_LENGTH, SIMULATION_LINE_WIDTH, SIMULATION_LINE_THICKNESS)

    marker_size = MARKERSIZE
    if len(VNA_data[0][0]) > 100:
        marker_size = 1

    fig1, ax1 = plt.subplots()
    ax1.set_title('S11 vs Frequency', fontsize=18)
    ax1.grid(True)
    ax1.set_xlabel('Frequency [GHz]')
    ax1.set_ylabel('|S11|')
    if sim_path != '':
        for i, d in enumerate(sim_data):
            ax1.errorbar(d[0], d[1], xerr=0, yerr=0, color=COLORS[i], linestyle='solid', marker='o', markersize=marker_size, label=f'Sim {label}')
    for i, d in enumerate(VNA_data):
        multiplied = []
        for number in d[0]:
            multiplied.append(number * 1e-9)
        ax1.errorbar(multiplied, d[1], xerr=0, yerr=0, color=COLORS[i+1], linestyle='solid', marker='o', markersize=marker_size, label=f'VNA {label}')
    ax1.legend(fontsize='x-large')

    if not one_port_only:
        fig2, ax2 = plt.subplots()
        ax2.set_title('S21 vs Frequency', fontsize=18)
        ax2.grid(True)
        ax2.set_xlabel('Frequency [GHz]')
        ax2.set_ylabel('|S21|')
        for i, d in enumerate(sim_data):
            ax2.errorbar(d[0], d[3], xerr=0, yerr=0, color=COLORS[i], linestyle='solid', marker='o', markersize=marker_size, label=f'Sim {label}')
        for i, d in enumerate(VNA_data):
            multiplied = []
            for number in d[0]:
                multiplied.append(number * 1e-9)
            ax2.errorbar(multiplied, d[2], xerr=0, yerr=0, color=COLORS[i+1], linestyle='solid', marker='o', markersize=marker_size, label=f'VNA {label}')
        ax2.legend(fontsize='x-large')

    if one_port_only:
        Z_modulo = []
        Z_phase = []

        fig3, ax3 = plt.subplots()
        ax3.set_title('|Z| vs Frequency', fontsize=18)
        ax3.grid(True)
        ax3.set_xlabel('Frequency [GHz]')
        ax3.set_ylabel('$|Z| [\\Omega]$')
        
        fig4, ax4 = plt.subplots()
        ax4.set_title('arg(Z) vs Frequency', fontsize=18)
        ax4.grid(True)
        ax4.set_xlabel('Frequency [GHz]')
        ax4.set_ylabel('arg(Z) [deg]')
        
        for i, data in enumerate(sim_data):
            for i_freq in enumerate(data[0]):
                i = i_freq[0]
                S11_modulo = 10**(data[1][i]/20)
                S11_phase = np.pi/180*data[2][i]
                S11_real = S11_modulo*np.cos(S11_phase)
                S11_imag = S11_modulo*np.sin(S11_phase)
                Z = SIMULATION_Z0 * (1 + S11_real + 1j*S11_imag) / (1 - S11_real - 1j*S11_imag)
                Z_modulo.append(np.sqrt(Z.real**2 + Z.imag**2))
                Z_phase.append(np.atan(Z.imag/Z.real)*180/np.pi)
            ax3.errorbar(data[0], Z_modulo, xerr=0, yerr=0, color='green', linestyle='solid', marker='o', markersize=marker_size, label=f'Sim {label}')
            ax4.errorbar(data[0], Z_phase, xerr=0, yerr=0, color='green', linestyle='solid', marker='o', markersize=marker_size, label=f'Sim {label}')

        try:
            VNA_Z_modulo = []
            VNA_Z_modulo.append(read_csv_data(f'{VNA_path}/TRACE02.CSV', [0, 1]))
            multiplied = []
            for number in VNA_Z_modulo[0][0]:
                multiplied.append(number * 1e-9)
            ax3.errorbar(multiplied, VNA_Z_modulo[0][1], xerr=0, yerr=0, color='red', linestyle='solid', marker='o', markersize=marker_size, label=f'VNA {label}')
            VNA_Z_fase = []
            VNA_Z_fase.append(read_csv_data(f'{VNA_path}/TRACE03.CSV', [0, 1]))
            ax4.errorbar(multiplied, VNA_Z_fase[0][1], xerr=0, yerr=0, color='red', linestyle='solid', marker='o', markersize=marker_size, label=f'VNA {label}')
        except FileNotFoundError:
            print("File TRACE02.CSV or TRACE03.CSV not found in VNA path. Skipping Z plot from VNA data.")
        
        ax3.legend(fontsize='x-large')
        ax4.legend(fontsize='x-large')

    plt.xticks(rotation=25)
    #plt.tight_layout()
    plt.show()




# --- MAIN EXECUTION ---
if __name__ == '__main__':

    sim_path = 'postpro/Gold/Lumped/CalKit/LOAD103R/TWOPORT_cpw_lumped_3.33um_50_50ohm_cond_45000000_port0.075_50um_NOPEC_Absorbing_4'
    #'postpro/Gold/Lumped/CalKit/0509/TWOPORT_cpw_lumped_10um_50_99999999ohm_cond_45000000_port0.075_1800um_NOPEC_Absorbing_4'
    #'postpro/Gold/Lumped/CalKit/0509/TWOPORT_cpw_lumped_8um_50_50ohm_cond_45000000_port0.075_1600um_NOPEC_Absorbing_4'
    #'postpro/Gold/Wave/CalKit/0509/ONEPORT_cpw_wave_10um_removeMetal_condSubstrate_port500x785_1800um_PEC_5_6_8_10_Absorbing_7'
    #'postpro/Gold/Wave/CalKit/0509/ONEPORT_cpw_wave_10um_cond_45000000_port500x785_1800um_PEC_5_6_8_10_Absorbing_7'
    #'postpro/Gold/Wave/CalKit/0509/ONEPORT_cpw_wave_10um_removeMetal_port500x785_1800um_PEC_7_9_NOAbsorbing'
    #'postpro/Gold/Wave/CalKit/0509/ONEPORT_cpw_wave_10um_removeMetal_port500x785_1800um_PEC_4_5_7_9_Absorbing_6'
    #'postpro/Gold/Lumped/CalKit/0509/TWOPORT_cpw_lumped_10um_50_99999999ohm_cond_45000000_port0.075_1800um_NOPEC_Absorbing_4'
    #'postpro/Gold/Lumped/CalKit/0509/TWOPORT_cpw_lumped_10um_50_378ohm_cond_45000000_port0.075_1800um_NOPEC_Absorbing_4'
    VNA_path = '/home/temp/Desktop/VNA/10-12-2025/2'
    #'/home/temp/Desktop/VNA/10-12-2025/1'
    #'/home/temp/Desktop/VNA/11-12-2025/0509'
    label = '[CalKitLOAD103R]'
    #'[CalKit0509-Open]'
    plot_sim_vs_VNA(sim_path, VNA_path, label, one_port_only=True)


    #========================================================================================================================================#


    #Open CPW
    plot_data_vs_frequency('postpro/Alluminio/Lumped/OpenCPW/cpw_lumped_22um_50ohm_cond_32894736_port0.075_8000um_NOPEC_Absorbing_4',
                           '[Kapton-Al, OpenCPW-8mm]',
                           one_port_only=True)
    
    #calculate_Z0_from_conductor_losses(trace_width=90e-6, trace_thickness=20e-6, gap_width=100e-6, substrate_thickness=25e-6, dielectric_constant=3.3)
    
    """
    #CPW with Zc = 100.9619 Ohm (Z0 = 126 Ohm, e_eff = 1.57) post matlab results
    plot_data_vs_frequency('postpro/Rame/Lumped/CPW/TestPostMatlab/cpw_lumped_22um_100.9619ohm_cond_59600000_port0.075_1000um_NOPEC_Absorbing_4',
                           '[Kapton-Copper, CPW-1mm]',
                           one_port_only=False)
    """
    
    #print(calculate_conductor_losses(1, CONDUCTIVITY_COPPER, 90e-6, 20e-6, 1, 100e-6, 105))

    #CalKit0509 with dimensions before SEM measurements
    plot_data_vs_frequency('postpro/Gold/Lumped/CalKit/0509/ONEPORT_cpw_lumped_8um_50ohm_cond_45000000_port0.075_1600um_NOPEC_Absorbing_4',
                           '[CalKit0509-Open]',
                           one_port_only=True)

    #CalKit0509 with real dimensions measured with SEM
    plot_data_vs_frequency('postpro/Gold/Lumped/CalKit/0509/ONEPORT_cpw_lumped_10um_50ohm_cond_45000000_port0.075_1800um_NOPEC_Absorbing_4',
                           '[CalKit0509-Open]',
                           one_port_only=True)
    
    
    
    

