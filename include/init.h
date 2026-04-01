#ifndef INIT_H
#define INIT_H

#include "sim_structs.h"

void create_particles(struct test_particles *part, struct parameters param);
void generate_random_particles(struct test_particles *part, struct parameters param);
void free_particles(struct test_particles *part);

#endif