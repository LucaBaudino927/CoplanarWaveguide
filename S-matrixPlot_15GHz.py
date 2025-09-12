import matplotlib.pyplot as plt
import csv
from pathlib import Path

# --- CONFIGURATION ---
DOUBLE_TRACES = False
COMPARE_SIMULATIONS = True
BASE_PATH_WAVE = Path('postpro/Rame/Wave')
BASE_PATH_LUMPED = Path('postpro/Rame/Lumped')
BASE_PATH_DIVERSE_Z = Path('postpro/Rame/Lumped/Mesh_22um_diverse_impedenze')
BASE_PATH_DIVERSE_MESH_122OHM = Path('postpro/Rame/Lumped/Impedenza_122ohm_diverse_mesh')
BASE_PATH_DIVERSE_MESH_102OHM = Path('postpro/Rame/Lumped/Impedenza_102ohm_diverse_mesh')

IMPEDANCES = [62, 92, 102, 112, 122, 132, 142]
MESHES = [15, 22, 45, 90]
COLORS = ['green', 'red', 'blue', 'orange', 'cyan', 'black', 'magenta']
SIERRA_PATH = BASE_PATH_LUMPED / 'SierraData_122ohm/port-S_1GHz_15GHz.csv'
TARGET_PORTS = ['S11', 'S21', 'S31', 'S41']
TRACE_FILE = BASE_PATH_WAVE / 'coplanar_waveguide_wave/port-S.csv'

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


def plot_curve(x, ys, labels, colors, title, ylabel):
    """
    Plots multiple curves on the same figure.
    """
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
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_single_trace(path, double_traces=False):
    ports_idx = [0, 1, 3]
    if double_traces:
        ports_idx += [5, 7]

    data = read_csv_data(path, ports_idx)
    x = data[0]

    port_labels = TARGET_PORTS[:len(data) - 1]
    port_data = data[1:]

    count = 0
    for label, y in zip(port_labels, port_data):
        plot_curve(x, [y], [label], ['blue'], f"{label} vs Frequency", f"{label} [dB]")
        count += 1


def plot_multiple_simulations(base_path):
    x = []
    s11_all = []
    s21_all = []
    labels = []

    if base_path == BASE_PATH_DIVERSE_Z:
        labels = [str(imp) for imp in IMPEDANCES]
        for i, imp in enumerate(IMPEDANCES):
            path = base_path / f'coplanar_waveguide_lumped_22um_{imp}ohm/port-S.csv'
            data = read_csv_data(path, [0, 1, 3])
            if imp == 62:
                x = data[0]
            s11_all.append(data[1])
            s21_all.append(data[2])
        

    if base_path == BASE_PATH_DIVERSE_MESH_122OHM:
        labels = [str(mesh) + 'um' for mesh in MESHES]
        for i, mesh in enumerate(MESHES):
            path = base_path / f'coplanar_waveguide_lumped_{mesh}um_122ohm/port-S.csv'
            data = read_csv_data(path, [0, 1, 3])
            if mesh == 90:
                x = data[0]
            s11_all.append(data[1])
            s21_all.append(data[2])

    if base_path == BASE_PATH_DIVERSE_MESH_102OHM:
        labels = [str(mesh) + 'um' for mesh in MESHES]
        for i, mesh in enumerate(MESHES):
            path = base_path / f'coplanar_waveguide_lumped_15um_102ohm/port-S.csv'
            data = read_csv_data(path, [0, 1, 3])
            x = data[0]
            s11_all.append(data[1])
            s21_all.append(data[2])

    # Sierra data
    sierra_s11, sierra_s21 = read_csv_data(SIERRA_PATH, [1, 3])

    # Plot S11
    plot_curve(x, s11_all + [sierra_s11], labels + ['Sierra'], COLORS + ['r'], 'S11 vs Frequency', 'S11 [dB]')

    # Plot S21
    plot_curve(x, s21_all + [sierra_s21], labels + ['Sierra'], COLORS + ['r'], 'S21 vs Frequency', 'S21 [dB]')

# --- MAIN EXECUTION ---
if __name__ == '__main__':
    if COMPARE_SIMULATIONS:
        plot_multiple_simulations(BASE_PATH_DIVERSE_Z)
        plot_multiple_simulations(BASE_PATH_DIVERSE_MESH_122OHM)
        plot_multiple_simulations(BASE_PATH_DIVERSE_MESH_102OHM)
    else:
        plot_single_trace(TRACE_FILE, double_traces=DOUBLE_TRACES)