#include <math.h>
#include <stdio.h>

#include "init.h"
#include "tools.h"
#include "physics.h"
#include "math_tools.h"
#include "sim_structs.h"

void initialize_particles(struct test_particles *part, struct parameters param, struct woods_saxon *ws, struct skyrme skm, struct fermi *fermi_levels) {
	double delta_part_p, delta_part_n, delta_epsilon_p, delta_epsilon_n;
	double r_max = param.r_max, sigma_k = param.sigma_k, sigma_r = param.sigma_r, total_delta_epsilon;
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
		
		delta_part_n = total_n - check_equal_n;
		delta_part_p = total_p - check_equal_p;
		delta_epsilon_n = 0.5 * delta_part_n / (check_more_n - check_less_n);
		delta_epsilon_p = 0.5 * delta_part_p / (check_more_p - check_less_p);
		
		if(fabs(delta_epsilon_p) > 0.5) delta_epsilon_p *= 0.6;
		if(fabs(delta_epsilon_n) > 0.5) delta_epsilon_n *= 0.6;
		
		fermi_levels->epsilon_p += delta_epsilon_p;
		fermi_levels->epsilon_n += delta_epsilon_n;
		total_delta_epsilon = fabs(delta_epsilon_n) + fabs(delta_epsilon_p);
		
		generate_checking_particles(part, ws, param, fermi_levels);
		compute_particle_densities(part, sigma_r, (double)part_per_nucleon);
		
		//minim_woods_saxon(part, ws, skm);
		
		printf("ITERATION = %i DELTA EPSILON = %lf\n", it,  total_delta_epsilon);
		printf("%i %i %i\n", check_less_p, check_equal_p, check_more_p);
		printf("%i %i %i\n", check_less_n, check_equal_n, check_more_n);
		
		it++;
	} while(total_delta_epsilon > DELTA_EPSILON_TOLERANCE && it < MAX_ITERATIONS);
	
	free_particles(part);
	create_particles(part, total_p, total_n);
	generate_checking_particles(part, ws, param, fermi_levels);
	compute_particle_energies(part, ws, param);
	compute_particle_densities(part, sigma_r, (double)part_per_nucleon);
	
	printf("WS PARAM %lf %lf %lf\n", ws[0].V0, ws[1].R12, ws[1].a);
	printf("WS PARAM %lf %lf %lf\n", ws[1].V0, ws[1].R12, ws[1].a);
	
	//fit_woods_saxon_param(&ws_p, part_p, skm, total_p);
	//fit_woods_saxon_param(&ws_n, part_n, skm, total_n);
	
	//printf("WS PARAM %lf %lf %lf\n", ws[0]->V0, ws[1]->R12, ws[1]->a);
	//printf("WS PARAM %lf %lf %lf\n", ws[1]->V0, ws[1]->R12, ws[1]->a);
	
	chi_squared(part, *ws, skm);
}