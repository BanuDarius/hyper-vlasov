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

#ifndef CUDA_INIT_H
#define CUDA_INIT_H

#include <cuda_runtime.h>

#include "sim_structs.hpp"

template <typename T>
void gpu_create_scalar_field_single(ScalarField<T> *d_field, const World<T> &world) {
	size_t world_size = world.n[0] * world.n[1] * world.n[2];
	cudaMalloc((void**)&d_field->v, world_size * sizeof(T));
	cudaMemset(d_field->v, 0, world_size * sizeof(T));
}

template <typename T>
void gpu_create_scalar_field_double(ScalarField<T> *d_field, const World<T> &world) {
	size_t world_size = world.n[0] * world.n[1] * world.n[2];
	cudaMalloc((void**)&d_field->v, 2 * world_size * sizeof(T));
	cudaMemset(d_field->v, 0, 2 * world_size * sizeof(T));
}

template <typename T>
void gpu_create_vector_field_double(VectorField<T> *d_field, const World<T> &world) {
	size_t world_size = world.n[0] * world.n[1] * world.n[2];
	cudaMalloc((void**)&d_field->x, 2 * world_size * sizeof(T));
	cudaMalloc((void**)&d_field->y, 2 * world_size * sizeof(T));
	cudaMalloc((void**)&d_field->z, 2 * world_size * sizeof(T));
	cudaMemset(d_field->x, 0, 2 * world_size * sizeof(T));
	cudaMemset(d_field->y, 0, 2 * world_size * sizeof(T));
	cudaMemset(d_field->z, 0, 2 * world_size * sizeof(T));
}

template <typename T>
void gpu_create_particles(TestParticles<T> *d_part, int protons, int neutrons) {
	int total = protons + neutrons;
	d_part->protons = protons;
	d_part->neutrons = neutrons;
	
	cudaMalloc((void**)&d_part->x, total * sizeof(T));
	cudaMalloc((void**)&d_part->y, total * sizeof(T));
	cudaMalloc((void**)&d_part->z, total * sizeof(T));
	cudaMalloc((void**)&d_part->kx, total * sizeof(T));
	cudaMalloc((void**)&d_part->ky, total * sizeof(T));
	cudaMalloc((void**)&d_part->kz, total * sizeof(T));
	cudaMalloc((void**)&d_part->fx, total * sizeof(T));
	cudaMalloc((void**)&d_part->fy, total * sizeof(T));
	cudaMalloc((void**)&d_part->fz, total * sizeof(T));
	
	cudaMemset(d_part->x, 0, total * sizeof(T));
	cudaMemset(d_part->y, 0, total * sizeof(T));
	cudaMemset(d_part->z, 0, total * sizeof(T));
	cudaMemset(d_part->kx, 0, total * sizeof(T));
	cudaMemset(d_part->ky, 0, total * sizeof(T));
	cudaMemset(d_part->kz, 0, total * sizeof(T));
	cudaMemset(d_part->fx, 0, total * sizeof(T));
	cudaMemset(d_part->fy, 0, total * sizeof(T));
	cudaMemset(d_part->fz, 0, total * sizeof(T));
}

template <typename T>
void gpu_free_particles(TestParticles<T> *part) {
	cudaFree(part->x); cudaFree(part->y); cudaFree(part->z);
	cudaFree(part->kx); cudaFree(part->ky); cudaFree(part->kz);
	cudaFree(part->fx); cudaFree(part->fy); cudaFree(part->fz);
}

template <typename T>
void gpu_free_vector_field(VectorField<T> *field) {
	cudaFree(field->x); free(field->y); free(field->z);
}

template <typename T>
void gpu_free_scalar_field(ScalarField<T> *field) {
	cudaFree(field->v);
}

template <typename T>
void copy_scalar_field_gpu_to_cpu(ScalarField<T> *h_field, ScalarField<T> *d_field, const World &world) {
	size_t world_size = world.n[0] * world.n[1] * world.n[2];
	cudaMemcpy(h_field->v, d_field->v, world_size * sizeof(T), cudaMemcpyDeviceToHost);
}

template <typename T>
void copy_vector_field_gpu_to_cpu(VectorField<T> *h_field, VectorField<T> *d_field, const World &world) {
	size_t world_size = world.n[0] * world.n[1] * world.n[2];
	cudaMemcpy(h_field->x, d_field->x, world_size * sizeof(T), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_field->y, d_field->y, world_size * sizeof(T), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_field->z, d_field->z, world_size * sizeof(T), cudaMemcpyDeviceToHost);
}

template <typename T>
void copy_particles_cpu_to_gpu(TestParticles<T> *d_part, TestParticles<T> *h_part) {
	size_t total = h_part->protons + h_part->neutrons;
	cudaMemcpy(d_part->x, h_part->x, total * sizeof(T), cudaMemcpyHostToDevice);
	cudaMemcpy(d_part->y, h_part->y, total * sizeof(T), cudaMemcpyHostToDevice);
	cudaMemcpy(d_part->z, h_part->z, total * sizeof(T), cudaMemcpyHostToDevice);
	cudaMemcpy(d_part->kx, h_part->kx, total * sizeof(T), cudaMemcpyHostToDevice);
	cudaMemcpy(d_part->ky, h_part->ky, total * sizeof(T), cudaMemcpyHostToDevice);
	cudaMemcpy(d_part->kz, h_part->kz, total * sizeof(T), cudaMemcpyHostToDevice);
	cudaMemcpy(d_part->fx, h_part->fx, total * sizeof(T), cudaMemcpyHostToDevice);
	cudaMemcpy(d_part->fy, h_part->fy, total * sizeof(T), cudaMemcpyHostToDevice);
	cudaMemcpy(d_part->fz, h_part->fz, total * sizeof(T), cudaMemcpyHostToDevice);
}

template <typename T>
void copy_particles_gpu_to_cpu(TestParticles<T> *h_part, TestParticles<T> *d_part) {
	size_t total = d_part->protons + d_part->neutrons;
	cudaMemcpy(h_part->x, d_part->x, total * sizeof(T), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_part->y, d_part->y, total * sizeof(T), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_part->z, d_part->z, total * sizeof(T), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_part->kx, d_part->kx, total * sizeof(T), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_part->ky, d_part->ky, total * sizeof(T), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_part->kz, d_part->kz, total * sizeof(T), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_part->fx, d_part->fx, total * sizeof(T), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_part->fy, d_part->fy, total * sizeof(T), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_part->fz, d_part->fz, total * sizeof(T), cudaMemcpyDeviceToHost);
}

#endif