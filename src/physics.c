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

#include <math.h>
#include <stdio.h>

#include "init.h"
#include "tools.h"
#include "physics.h"
#include "math_tools.h"
#include "sim_structs.h"
#include "fit_algorithm.h"

void initialize_particles(struct test_particles *part, struct parameters param, struct woods_saxon *ws, struct skyrme skm, struct fermi *fermi_levels) {
	double r_max = param.r_max, total_delta_epsilon, relax_coef = 0.6;
	int max_part = param.max_test_part, z = param.z, n = param.n, part_per_nucleon = param.test_part_per_nucleon, it = 0;
	int total_p = z * part_per_nucleon, total_n = n * part_per_nucleon;
	
	struct test_particles temp_part;
	create_particles(part, total_p, total_n);
	create_particles(&temp_part, max_part, max_part);
	do {
		generate_random_particles(&temp_part, r_max);
		compute_particle_energies(&temp_part, ws, param);
		int check_less_p = 0, check_equal_p = 0, check_more_p = 0;
		int check_less_n = 0, check_equal_n = 0, check_more_n = 0;
		for(int i = 0; i < max_part; i++) {
			if(temp_part.energy[i] < fermi_levels->epsilon_p)
				check_equal_p += 2;
			if(temp_part.energy[i] < fermi_levels->epsilon_p + 0.5)
				check_more_p += 2;
			if(temp_part.energy[i] < fermi_levels->epsilon_p - 0.5)
				check_less_p += 2;
			
			if(temp_part.energy[i + max_part] < fermi_levels->epsilon_n)
				check_equal_n += 2;
			if(temp_part.energy[i + max_part] < fermi_levels->epsilon_n + 0.5)
				check_more_n += 2;
			if(temp_part.energy[i + max_part] < fermi_levels->epsilon_n - 0.5)
				check_less_n += 2;
		}
		double delta_part_n = total_n - check_equal_n;
		double delta_part_p = total_p - check_equal_p;
		double delta_epsilon_n = 0.5 * delta_part_n / (check_more_n - check_less_n);
		double delta_epsilon_p = 0.5 * delta_part_p / (check_more_p - check_less_p);
		
		if(fabs(delta_epsilon_p) > 0.5) delta_epsilon_p *= relax_coef;
		if(fabs(delta_epsilon_n) > 0.5) delta_epsilon_n *= relax_coef;
		
		fermi_levels->epsilon_p += delta_epsilon_p;
		fermi_levels->epsilon_n += delta_epsilon_n;
		
		generate_checking_particles(part, ws, param, fermi_levels);
		compute_particle_densities(part, param);
		
		struct woods_saxon ws_old[2];
		ws_old[0] = ws[0]; ws_old[1] = ws[1];
		
		minim_woods_saxon(part, ws, skm);
		relax_woods_saxon(ws, ws_old, relax_coef);
		
		total_delta_epsilon = fabs(delta_epsilon_n) + fabs(delta_epsilon_p);
		printf("%i %i %i\n", check_less_p, check_equal_p, check_more_p);
		printf("%i %i %i\n", check_less_n, check_equal_n, check_more_n);
		printf("WS PARAM %0.2lf %0.2lf %0.2lf\n", ws[0].V0, ws[0].R12, ws[0].a);
		printf("WS PARAM %0.2lf %0.2lf %0.2lf\n", ws[1].V0, ws[1].R12, ws[1].a);
		printf("DELTA EPSILON = %0.2lf\nITERATION = %i\n", total_delta_epsilon, it);
		
		it++;
	} while(total_delta_epsilon > DELTA_EPSILON_TOLERANCE && it < MAX_ITERATIONS);	
	compute_particle_energies(part, ws, param);
}

double nuclear_radius(unsigned short a) {
	double radius = 1.5 * pow((double)a, 1.0 / 3.0);
	return radius;
}

int max_particles(double r_max, double k_max, int test_part_per_nucleon) {
	double t = r_max * k_max, ct = 2.0 * M_PI;
	double phase_space_volume = (16.0 / 9.0) * M_PI * M_PI * (t * t * t);
	int max = test_part_per_nucleon * (int)floor(phase_space_volume / (ct * ct * ct) + 0.5);
	return max;
}

double woods_saxon_potential(struct woods_saxon ws, double r) {
	double v = ws.V0 / (1.0 + exp((r - ws.R12) / ws.a));
	return v;
}

double skyrme_potential(struct skyrme skm, double rho_p, double rho_n, int type) {
	double tau = (type == PROTONS) ? -1.0 : +1;
	double rho = rho_p + rho_n;
	double t = rho / RHO_0;
	double v = skm.A * t + skm.B * pow(t, skm.gamma) + tau * skm.C * ((rho_n - rho_p) / RHO_0);
	return v;
}

double coulomb_potential(struct woods_saxon ws, double z, double r) {
	double R12 = ws.R12, v;
	if(r <= ws.R12)
		v = 1.44 * (z - 1.0) / R12 * (1.5 - 0.5 * (r / R12) * (r / R12));
	else
		v = 1.44 * (z - 1.0) / r;
	return v;
}