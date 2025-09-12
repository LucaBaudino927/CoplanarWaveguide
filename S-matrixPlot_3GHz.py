import matplotlib.pyplot as plt
import csv
from pathlib import Path
import numpy as np

# --- CONFIGURATION ---
COMPARE_SIMULATIONS = True
BASE_PATH_WAVE = Path('postpro/Rame/Wave')
BASE_PATH_LUMPED = Path('postpro/Rame/Lumped')
BASE_PATH_102OHM = Path('postpro/Rame/Lumped/Impedenza_102ohm')
BASE_PATH_102OHM_DIVERSE_CONDUCIBILITA = Path('postpro/Rame/Lumped/Impedenza_102ohm_mesh_22um_diverse_conducibilita')
BASE_PATH_TEST = Path('postpro/Rame/Lumped/Test_diversi_parametri_e_geometrie')

IMPEDANCES = [100, 102, 102.1, 102.5]
#IMPEDANCES = [102, 122]
CONDUCTIVITIES = [33112582, 57471264, 59600000, 62500000, 72500000] #S/m
COLORS = ['green', 'red', 'blue', 'orange', 'cyan', 'black', 'magenta', 'purple', 'brown', 'pink', 'gray']
SIERRA_PATH = BASE_PATH_LUMPED / 'SierraData_122ohm'
TARGET_PORTS = ['S11', 'S21']
TRACE_FILE = BASE_PATH_LUMPED / 'Test_diversi_parametri_e_geometrie/coplanar_waveguide_lumped_22um_102ohm_cond_62500000_porte_50um/port-S.csv'

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


def plot_curve(x, ys, labels, colors, title, ylabel, figure_id, skin_effect_thr=-1, show_loss_thr=True):
    """
    Plots multiple curves on the same figure.
    """
    plt.figure(figure_id)
    for y, label, color in zip(ys, labels, colors):
        if label == 'Sierra':
            plt.plot(x, y, label=label, color=color, linestyle='solid', marker='x')
        else:
            plt.plot(x, y, label=label, color=color, linestyle='solid', marker='o')
    plt.title(title, fontsize=18)
    plt.xlabel('Frequency [GHz]')
    plt.ylabel(ylabel)
    plt.xticks(rotation=25)
    plt.grid(True)
    plt.tight_layout()
    
    if skin_effect_thr > 0:
        plt.axvline(skin_effect_thr, color='green', linestyle='--', label='Skin Effect Threshold')

    plt.legend()

    if "S21" in title:
        if show_loss_thr:
            dielectric_loss_thr = []
            conductor_loss_thr = []
            print("\n===========Calculating dielectric losses===========")
            dielectric_loss_thr.extend(calculate_dielectric_losses(i * 1e9) for i in x)
            print("\n\n===========Calculating conductor losses===========")
            conductor_loss_thr.extend(calculate_conductor_losses(i * 1e9) for i in x)
            for i in range(len(x)):
                if i < len(x) - 1:
                    plt.hlines(-dielectric_loss_thr[i], x[i], x[i+1], color='blue', linestyle='--', label='Dielectric Loss Threshold')
                    plt.fill_betweenx([-dielectric_loss_thr[i], 0], x[i], x[i+1], color='blue', alpha=0.1)
                    plt.hlines(-dielectric_loss_thr[i]-conductor_loss_thr[i], x[i], x[i+1], color='red', linestyle='--', label='Conductor Loss Threshold')
                    plt.fill_betweenx([-dielectric_loss_thr[i]-conductor_loss_thr[i], -dielectric_loss_thr[i]], x[i], x[i+1], color='red', alpha=0.1)
                else:
                    plt.hlines(-dielectric_loss_thr[i], x[i], 4, color='blue', linestyle='--', label='Dielectric Loss Threshold')
                    plt.fill_betweenx([-dielectric_loss_thr[i], 0], x[i], 4, color='blue', alpha=0.1)
                    plt.hlines(-dielectric_loss_thr[i]-conductor_loss_thr[i], x[i], 4, color='red', linestyle='--', label='Conductor Loss Threshold')
                    plt.fill_betweenx([-dielectric_loss_thr[i]-conductor_loss_thr[i], -dielectric_loss_thr[i]], x[i], 4, color='red', alpha=0.1)             
    
    #plt.xscale('log')


def plot_single_trace(path):
    ports_idx = [0, 1, 3]

    data = read_csv_data(path, ports_idx)
    x = data[0]

    port_labels = TARGET_PORTS[:len(data) - 1]
    port_data = data[1:]

    sierra_path = SIERRA_PATH / 'port-S_10MHz_3GHz.csv'
    # Sierra data
    sierra_s11, sierra_s21 = read_csv_data(sierra_path, [1, 3])

    count = 0
    for label, y in zip(port_labels, port_data):
        if "S11" in label:
            plot_curve(x, [y], [label], ['blue'], f"{label} vs Frequency", f"{label} [dB]", count)
            plot_curve(x, [sierra_s11], ['Sierra'], ['red'], f"{label} vs Frequency", f"{label} [dB]", count)
        if "S21" in label:
            plot_curve(x, [y], [label], ['blue'], f"{label} vs Frequency", f"{label} [dB]", count)
            plot_curve(x, [sierra_s21], ['Sierra'], ['red'], f"{label} vs Frequency", f"{label} [dB]", count)
        count += 1

    plt.show()


def plot_multiple_simulations(base_path):
    x = []
    s11_all = []
    s21_all = []
    labels = []
    sierra_path = ''

    if base_path == BASE_PATH_102OHM or base_path == BASE_PATH_TEST:
        sierra_path = SIERRA_PATH / 'port-S_10MHz_3GHz.csv'
        labels = [str(imp) for imp in IMPEDANCES]
        for i, imp in enumerate(IMPEDANCES):
            path = base_path / f'coplanar_waveguide_lumped_22um_{imp}ohm_cond_62500000/port-S.csv'
            data = read_csv_data(path, [1, 3])
            s11_all.append(data[0])
            s21_all.append(data[1])

    if base_path == BASE_PATH_102OHM_DIVERSE_CONDUCIBILITA:
        sierra_path = SIERRA_PATH / 'port-S_500MHz_3GHz.csv'
        labels = [str(cond) for cond in CONDUCTIVITIES]
        for i, cond in enumerate(CONDUCTIVITIES):
            path = base_path / f'coplanar_waveguide_lumped_22um_102ohm_cond_{cond}/port-S.csv'
            data = read_csv_data(path, [1, 3])
            s11_all.append(data[0])
            s21_all.append(data[1])

    # Sierra data
    x, sierra_s11, sierra_s21 = read_csv_data(sierra_path, [0, 1, 3])

    # Calculate skin effect threshold
    skin_effect_thr = calculate_skin_effect_threshold(0.3)
    print(f"Skin effect threshold frequency: {skin_effect_thr} GHz")
    # Plot S11
    plot_curve(x, s11_all + [sierra_s11], labels + ['Sierra'], COLORS + ['r'], 'S11 vs Frequency', 'S11 [dB]', 0, skin_effect_thr)

    # Plot S21
    plot_curve(x, s21_all + [sierra_s21], labels + ['Sierra'], COLORS + ['r'], 'S21 vs Frequency', 'S21 [dB]', 1, skin_effect_thr)

    plt.show()

def calculate_skin_effect_threshold(x = 0.1, conductivity=62500000, trace_length=1000e-6, trace_width=90e-6, trace_thickness=20e-6):
    
    """
    Metodo 1:
    Calculate the skin effect threshold frequency for a given conductivity and trace properties.
    We consider the skin effect relevant when the resistence become x higher than ohmic resistance. Default x=0.1 means 10% increase.
    """
    k = np.sqrt(4 * 1e-7 * np.pi * np.pi * 0.999994)  # Assuming mu_r = 0.999994 for copper
    perimeter = 2 * (trace_thickness + trace_width)
    R = trace_length / (conductivity * trace_thickness * trace_width)
    frequency = conductivity * np.pow(R * (x + 1), 2) * np.pow(perimeter, 2) / (np.pow(trace_length, 2) * np.pow(k, 2))

    #return frequency/1e9  # Convert to GHz

    """
    Metodo 2:
    Calculate the skin effect threshold frequency for a given conductivity and trace properties.
    We consider the skin effect relevant when the Ac become x higher than Ac at low frequency. Default x=0.1 means 10% increase.
    """
    k = np.sqrt(4 * 1e-7 * np.pi * np.pi * 0.999994)  # Assuming mu_r = 0.999994 for copper
    perimeter = 2 * (trace_thickness + trace_width)
    R = trace_length / (conductivity * trace_thickness * trace_width)
    frequency = np.pow(x + 1, 2) / (conductivity * np.pi * 4e-7 * np.pi *np.pow(trace_thickness, 2))

    return frequency/1e9  # Convert to GHz

def calculate_dielectric_losses(frequency=3e9, trace_length=1000e-6, trace_width=90e-6, gap_width=100e-6, trace_thickness=20e-6, dielectric_constant=3.3, loss_tangent=0.0013):
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

    # Dielectric loss in dB/cm
    Ad = (8.686 * np.pi) * np.sqrt(e_eff) * loss_tangent * frequency / 2.99792458e8
    print(f"Frequency: {frequency/1e9:.2f} GHz, Effective Dielectric Constant: {e_eff:.4f}, Dielectric Loss: {1e-2 * Ad} dB/cm,  Dielectric Loss: {Ad * trace_length} dB")

    return Ad * (trace_length)  # Convert to dB for the given trace length in meters

def calculate_conductor_losses(frequency=3e9, conductivity=62500000, trace_width=90e-6, trace_length=1000e-6, Z0=102.1):
    """
    Fraction of power loss in the conductor for a CPW line.
    It calculates Ac in dB/m and return the value in dB.
    """
    mu0 = 4e-7 * np.pi
    delta = np.sqrt(1 / (np.pi * frequency * mu0 * conductivity * 0.999994))  # Skin depth

    Rs = 1 / (conductivity * delta) # Surface resistance in Ohm
    Rstrip = Rs / trace_width # Ohm/m
    # Conductor loss in dB/m
    Ac = (8.686 * Rstrip) / (2 * Z0)
    print(f"Frequency: {frequency/1e9:.2f} GHz, Skin Depth: {delta*1e6:.4f} um, Rs: {Rs:.4f} Ohm, Conductor Loss: {1e-2 * Ac} dB/cm, Conductor Loss: {Ac * trace_length} dB")

    return Ac * (trace_length)  # Convert to dB for the given trace length in meters
        


# --- MAIN EXECUTION ---
if __name__ == '__main__':
    if COMPARE_SIMULATIONS:
        plot_multiple_simulations(BASE_PATH_102OHM)
        #plot_multiple_simulations(BASE_PATH_102OHM_DIVERSE_CONDUCIBILITA)
        #plot_multiple_simulations(BASE_PATH_TEST)
    else:
        plot_single_trace(TRACE_FILE)
        