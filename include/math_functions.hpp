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

#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include <cmath>
#include <cstdlib>

template <typename T>
static inline T rand_val(T min, T max) {
	T s = std::rand() / T(RAND_MAX);
	return min + s * (max - min);
}

template <typename T>
static inline void random_vec(T *v, T max) {
	for(int i = 0; i < 3; i++)
		v[i] = rand_val(-max, max);
}

template <typename T>
static inline T dot(const T *a, const T *b) {
	T x = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return x;
}

template <typename T>
static inline T magnitude(const T *a) {
	T x = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	return x;
}

template <typename T>
static inline void mult_vec(T *a, const T *b, const T c) {
	for(int i = 0; i < 3; i++)
		a[i] = b[i] * c;
}

template <typename T>
static inline void add_vec(T *a, const T *b, const T *c) {
	for(int i = 0; i < 3; i++)
		a[i] = b[i] + c[i];
}

template <typename T>
static inline void sub_vec(T *a, const T *b, const T *c) {
	for(int i = 0; i < 3; i++)
		a[i] = b[i] - c[i];
}

#endif