'''MIT License

Copyright (c) 2026 Banu Darius-Matei

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.'''

# ---------------------------------------------------------- #

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})

# ---------------------------------------------------------- #

MAIN_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = MAIN_DIR.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
IMAGE_DIR = PROJECT_ROOT / "output-image"

# ---------------------------------------------------------- #

def plot_energy_spectrum():
    input_file = OUTPUT_DIR / "energy_spectrum.txt"
    output_image = IMAGE_DIR / "energy_spectrum.png"
    
    data = np.loadtxt(input_file)
    energy = data[:, 0]
    strength_function = data[:, 1]
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.plot(energy, strength_function, color='blue', linestyle='-', linewidth=1.5)
    
    ax.set_title("Energy Spectrum")
    ax.set_xlabel("E (MeV)")
    ax.set_ylabel(r"S(E) ($fm^2/MeV$)")
    ax.set_xlim(np.min(energy), 25.0)
    ax.set_ylim(np.min(strength_function), 1.1 * np.max(strength_function))
    
    plt.savefig(output_image, dpi=200, bbox_inches='tight')
    
    plt.close(fig)
    print(f"Created energy spectrum graph.")

# ---------------------------------------------------------- #

def plot_radius():
    input_file = OUTPUT_DIR / "stats.txt"
    output_image = IMAGE_DIR / "radius.png"
    
    data = np.loadtxt(input_file)
    time = data[:, 0]
    y_proton = data[:, 1]
    y_neutron = data[:, 2]
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.plot(time, y_proton, color='blue', linestyle='-', linewidth=1.5)
    ax.plot(time, y_neutron, color='red', linestyle='-', linewidth=1.5)
    
    ax.set_title("Nucleus proton neutron radii")
    ax.set_xlabel("t (MeV/c)")
    ax.set_ylabel("r (fm)")
    ax.set_xlim(0.0, np.max(time))
    ax.set_ylim(0.0, 1.1 * np.max(y_neutron))
    
    plt.savefig(output_image, dpi=200, bbox_inches='tight')
    
    plt.close(fig)
    print(f"Created radius graph.")
    
# ---------------------------------------------------------- #

def plot_center_of_mass():
    input_file = OUTPUT_DIR / "stats.txt"
    output_image = IMAGE_DIR / "center_of_mass.png"
    
    data = np.loadtxt(input_file)
    time = data[:, 0]
    y_proton = data[:, 3]
    y_neutron = data[:, 4]
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.plot(time, y_proton, color='blue', linestyle='-', linewidth=1.5)
    ax.plot(time, y_neutron, color='red', linestyle='-', linewidth=1.5)
    
    ax.set_title("Nucleus proton neutron center of mass")
    ax.set_xlabel("t (MeV/c)")
    ax.set_ylabel(r"$r_{CM}$ (fm)")
    ax.set_xlim(0.0, np.max(time))
    ax.set_ylim(1.3 * np.min(y_proton), 1.3 * np.max(y_neutron))
    
    plt.savefig(output_image, dpi=200, bbox_inches='tight')
    
    plt.close(fig)
    print(f"Created center of mass graph.")

# ---------------------------------------------------------- #
