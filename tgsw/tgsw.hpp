#ifndef TGSW_HPP
#define TGSW_HPP

#include <iostream>
#include <random>
#include <cstdint>
#include <vector>
#include <string.h>

using Torus32 = int32_t;


void sk_gen(int32_t* sk, int32_t n);
void print_matrix(Torus32** matrix, int32_t rows, int32_t cols);
void external_product(Torus32* result, Torus32** C, Torus32* lwe, int32_t n, int32_t bit_length, int32_t base_bit);
void LWE_encrypt(Torus32* result, int message, int32_t n, int32_t* sk);
void GSW_encrypt(Torus32** C, int32_t message, int32_t n, int32_t bit_length, int32_t Bksbit, int32_t* sk);
int32_t LWE_Decrypt(const Torus32 *sample, const int32_t* sk, int32_t n);

#endif // TGSW_HPP

