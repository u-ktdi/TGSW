// #ifndef TGSW_HPP
// #define TGSW_HPP

// #include <vector>
// #include <iostream>
// #include <cstdint>

// using namespace std;
// using Torus32 = int32_t;

// class TGSW {
//     public:
//         TGSW(int size);
//         void encrypt(const vector<Torus32>& data);
//         vector<Torus32> decrypt() const;
//         void print() const;

//     private:
//         vector<Torus32> encrypted_data;
// };

// #endif // TGSW_HPP

#ifndef TGSW_HPP
#define TGSW_HPP

#include <iostream>
#include <random>
#include <cstdint>
#include <vector>

using Torus32 = int32_t;

// TGSW関連の関数宣言
Torus32 modSwitchToTorus32(int32_t mu, int32_t Msize);
Torus32 dtot32(double d);
Torus32 gaussian32(int32_t message);
void LWE_encrypt(Torus32* result, int message, int32_t n, int32_t* sk);
void sk_gen(int32_t* sk, int32_t n);
void print_matrix(Torus32** matrix, int32_t rows, int32_t cols);
void G_inverse_of_scalar(Torus32* result, int32_t num, int32_t bit_length, int32_t base_bit);
void G_inverse_of_vector(Torus32* result, Torus32* ciphertext, int32_t n, int32_t bit_length, int32_t base_bit);
void G_Ginv_multiplication(Torus32* result, Torus32** G, Torus32* G_inv, int32_t n, int32_t bit_length);
int32_t vector_multiplication(Torus32* v1, Torus32* v2, int32_t size);
void matrix_vector_multiplication(Torus32* result, Torus32** mat, Torus32* vec, int32_t rows, int32_t cols);
void vector_matrix_multiplication(Torus32* result, Torus32* vec, Torus32** mat, int32_t mat_rows, int32_t mat_cols);
void create_G(Torus32** G, int32_t n, int32_t bit_length, int32_t base_bit);
void GSW_encrypt(Torus32** C, int message, int32_t n, int32_t bit_length, int32_t Bksbit, int32_t* sk);
void add(Torus32** C1, Torus32** C2, int32_t rows, int32_t cols);
void external_product(Torus32* result, Torus32** C, int32_t C_rows, int32_t C_cols, Torus32* lwe, int32_t bit_length, int32_t base_bit);
Torus32 lwePhase(const Torus32* sample, const int32_t* sk, int32_t n);
int32_t LWE_Decrypt(const Torus32 *sample, const int32_t* sk, int32_t n);
Torus32* AND_GSW(Torus32* c1, Torus32* c2, int32_t n, int32_t bit_length, int32_t base_bit);
void GSW();

#endif // TGSW_HPP

