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

#ifndef CUDA_PHYSICS_H
#define CUDA_PHYSICS_H

#include <cuda_runtime.h>

#include "sim_structs.hpp"

template <typename T>
__global__ void gpu_update_momenta_half(TestParticles<T> *part, T dt) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int total = part->protons + part->neutrons;
	if(idx >= total)
		return;
	T fact = dt / (T(2.0) * h_bar_c<T>);
	part->kx[idx] += fact * part->fx[idx];
	part->ky[idx] += fact * part->fy[idx];
	part->kz[idx] += fact * part->fz[idx];
}

template <typename T>
__global__ void gpu_update_positions_full(TestParticles<T> *part, T dt) {
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int total = part->protons + part->neutrons;
	if(idx >= total)
		return;
	T fact = (dt * h_bar_c<T>) / mc2<T>;
	part->x[idx] += fact * part->kx[idx];
	part->y[idx] += fact * part->ky[idx];
	part->z[idx] += fact * part->kz[idx];
}

#endif