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

#ifndef CUDA_TOOLS_H
#define CUDA_TOOLS_H

#include <cuda_runtime.h>

#include "sim_structs.hpp"

template <typename T>
void center_momentum(TestParticles<T> *part) {
	int total = part->protons + part->neutrons;
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if(idx >= total)
		return;
	for(int i = 0; i < total; i++) {
		k_sum[0] += part->kx[i];
		k_sum[1] += part->ky[i];
		k_sum[2] += part->kz[i];
	}
	for(int i = 0; i < total; i++) {
		part->kx[i] -= k_sum[0] / total;
		part->ky[i] -= k_sum[1] / total;
		part->kz[i] -= k_sum[2] / total;
	}
}

#endif