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
from scipy.fft import fft, fftfreq
from scipy import signal
import subprocess
from pathlib import Path

h_bar_c = 197.33

# ---------------------------------------------------------- #

MAIN_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = MAIN_DIR.parent
BIN_DIR = PROJECT_ROOT / "bin"
OUTPUT_DIR = PROJECT_ROOT / "output"

# ---------------------------------------------------------- #

def output_init_file(sim_parameters):
    input_file = sim_parameters.input_file
    
    with open(input_file, "w") as file:
        file.write(f"a {sim_parameters.a}\n")
        file.write(f"A {sim_parameters.A}\n")
        file.write(f"B {sim_parameters.B}\n")
        file.write(f"C {sim_parameters.C}\n")
        file.write(f"z {sim_parameters.z}\n")
        file.write(f"n {sim_parameters.n}\n")
        file.write(f"nx {sim_parameters.nx}\n")
        file.write(f"V0 {sim_parameters.V0}\n")
        file.write(f"t_f {sim_parameters.t_f}\n")
        file.write(f"t_exc {sim_parameters.t_exc}\n")
        file.write(f"eta_exc {sim_parameters.eta_exc}\n")
        file.write(f"steps {sim_parameters.steps}\n")
        file.write(f"gamma {sim_parameters.gamma}\n")
        file.write(f"r_fwhm {sim_parameters.r_fwhm}\n")
        file.write(f"k_fwhm {sim_parameters.k_fwhm}\n")
        file.write(f"substeps {sim_parameters.substeps}\n")
        file.write(f"use_gpu {int(sim_parameters.use_gpu)}\n")
        file.write(f"epsilon_p {sim_parameters.epsilon_p}\n")
        file.write(f"epsilon_n {sim_parameters.epsilon_n}\n")
        file.write(f"d_max_scale {sim_parameters.d_max_scale}\n")
        file.write(f"num_test_part {sim_parameters.num_test_part}\n")
        
# ---------------------------------------------------------- #

def run_simulation(sim_parameters):
    output_init_file(sim_parameters)
    
    if(sim_parameters.use_floats == False):
        arguments = [BIN_DIR / "hyper_vlasov", sim_parameters.input_file, sim_parameters.output_directory]
    else:
        arguments = [BIN_DIR / "hyper_vlasov", "--float", sim_parameters.input_file, sim_parameters.output_directory]
    
    subprocess.run(arguments, text=True)
    
# ---------------------------------------------------------- #

def compute_energy_spectrum(sim_parameters):
    z = sim_parameters.z
    n = sim_parameters.n
    t_f = sim_parameters.t_f
    steps = sim_parameters.steps
    t_exc = sim_parameters.t_exc
    substeps = sim_parameters.substeps
    
    stats_filename = OUTPUT_DIR / "stats.txt"
    
    data = np.loadtxt(stats_filename, dtype=np.float64)
    
    start_idx = int(t_exc / t_f * steps / substeps)
    time = data[start_idx:, 0]
    cm_protons = data[start_idx:, 3]
    cm_neutrons = data[start_idx:, 4]
    
    dipole = (n * z) / (n + z) * (cm_protons - cm_neutrons)
    signal.detrend(dipole)
    
    t0 = time[0]
    dt = t_f / steps * substeps
    dipole *= np.cos(np.pi * (time - t0) / (2.0 * (t_f - t0))) ** 2.0
    
    freq = fftfreq(np.size(dipole), d=dt)
    omega = 2.0 * np.pi * freq
    
    strength_raw = np.conjugate(fft(dipole)) * dt
    if(sim_parameters.eta_exc > 1e-5):
        strength = np.imag(strength_raw) / (np.pi * sim_parameters.eta_exc * h_bar_c)
    else:
        strength = np.imag(strength_raw)/ (np.pi * h_bar_c)
    
    pos_mask = omega > 0
    energy = omega[pos_mask] * h_bar_c
    strength = strength[pos_mask]
    
    out_file = OUTPUT_DIR / "energy_spectrum.txt"
    np.savetxt(out_file, np.column_stack((energy, strength)), fmt="%e")
    
    print("Computed energy spectrum.")
    
# ---------------------------------------------------------- #
