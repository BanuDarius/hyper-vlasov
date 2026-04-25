/* MIT License

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
SOFTWARE. */

#ifndef PHYSICS_H
#define PHYSICS_H

#include "sim_structs.h"

void initialize_particles(TestParticles *part, Parameters param, WoodsSaxon *ws, Skyrme skm, Fermi *fermi_levels);
void compute_volumetric_coulomb_potentials_sor(ScalarField *coulomb, ScalarField volume, World world);
void compute_volumetric_forces_fdm(VectorField *forces, ScalarField potentials, World world);
void compute_volumetric_skyrme_potentials(ScalarField *potentials, ScalarField volume, Skyrme skm, World world);
void update_momenta_half(TestParticles *part, double dt);
void update_positions_full(TestParticles *part, double dt);
double nuclear_radius(unsigned short a);
int max_particles(double r_max, double k_max, int total_test_part);

static inline double woods_saxon_potential(WoodsSaxon ws, double r) {
	double v = ws.V0 / (1.0 + exp((r - ws.R12) / ws.a));
	return v;
}

static inline double skyrme_potential(Skyrme skm, double rho_p, double rho_n, int type) {
	double tau = (type == PROTONS) ? -1.0 : +1.0;
	double rho = rho_p + rho_n;
	double t = rho / RHO_0;
	double v = skm.A * t + skm.B * pow(t, skm.gamma) + 2.0 * tau * skm.C * ((rho_n - rho_p) / RHO_0);
	return v;
}

static inline double coulomb_potential(WoodsSaxon ws, double z, double r) {
	double R12 = ws.R12, v;
	if(r <= ws.R12)
		v = 1.44 * (z - 1.0) / R12 * (1.5 - 0.5 * (r / R12) * (r / R12));
	else
		v = 1.44 * (z - 1.0) / r;
	return v;
}

#endif