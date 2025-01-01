/*
This is an implementation of the GSW scheme using the NTL and GMP C++ libraries.

Sources:
NTL C++ library: https://libntl.org/
GMP C++ library: https://gmplib.org/

Articles and references:
GSW scheme: https://eprint.iacr.org/2013/340.pdf

other related articles:
https://web.eecs.umich.edu/~cpeikert/pubs/polyboot.pdf
https://eprint.iacr.org/2021/691.pdf
https://eprint.iacr.org/2020/086.pdf
*/

#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <iomanip>
#include <cstdlib>

using namespace std;

// print matrix
void print_matrix(int32_t** matrix, size_t rows, size_t cols) {
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// create a matrix G^-1(scalar)
int32_t* G_inverse_of_scalar(int32_t num, int32_t bit_length, int32_t base_bit) {
    const int32_t base = 1 << base_bit;
    // const int32_t prec_offset = 1 << (32 - (1 + base_bit * bit_length)); // precision
    const int32_t mask = base - 1;

    // num += prec_offset;

    int32_t* num_bitdecomposed = (int32_t*)malloc(bit_length * sizeof(int32_t));
    for (size_t i = 0; i < bit_length; i++) {
        num_bitdecomposed[i] = (num >> (32 - (i + 1) * base_bit)) & mask;
    }

    return num_bitdecomposed;
}


// create a matrix G^-1(vector)
int32_t* G_inverse_of_vector(int32_t* ciphertext, size_t n, int32_t bit_length, int32_t base_bit) {
    size_t result_size = (n + 1) * bit_length;
    int32_t* bitdecomposed_ciphertext = (int32_t*)malloc(result_size * sizeof(int32_t));
    size_t index = 0;

    for (size_t j = 0; j < n + 1; j++) {
        int32_t* num_bitdecomposed = G_inverse_of_scalar(ciphertext[j], bit_length, base_bit);
        for (size_t k = 0; k < bit_length; k++) {
            bitdecomposed_ciphertext[index++] = num_bitdecomposed[k];
        }
        free(num_bitdecomposed);
    }

    return bitdecomposed_ciphertext;
}


// G * G^-1 
int32_t* G_Ginv_multiplication(int32_t** G, int32_t* G_inv, size_t n, int32_t bit_length) {
    size_t cols = n + 1;
    int32_t* result = (int32_t*)calloc(cols, sizeof(int32_t));

    size_t j = 0, rows = (n + 1) * bit_length;
    for (size_t i = 0; i < rows; i++) {
        result[j] += G[i][j] * G_inv[i];
        if ((i != 0) && (i % bit_length == 0)) j++;
    }

    return result;
}


// ベクトルの内積 (ポインタ版)
int32_t vector_multiplication(int32_t* v1, int32_t* v2, size_t n) {
    int32_t result = 0, size = n;

    for (size_t i = 0; i < size; i++) {
        result += v1[i] * v2[i];
    }

    return result;
}

// matrix vector multiplication
int32_t* matrix_vector_multiplication(int32_t** mat, int32_t* vec, size_t rows, size_t cols) {
    int32_t* result = (int32_t*)calloc(rows, sizeof(int32_t));
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            result[i] += mat[i][j] * vec[j];
        }
    }
    return result;
}


// vector matrix multiplication
int32_t* vector_matrix_multiplication(int32_t* vec, int32_t** mat, size_t vec_size, size_t mat_cols) {
    int32_t* result = (int32_t*)calloc(mat_cols, sizeof(int32_t));

    for (size_t i = 0; i < mat_cols; i++) {
        for (size_t j = 0; j < vec_size; j++) {
            result[i] += vec[j] * mat[j][i];
        }
    }

    return result;
}


// create matrix G
int32_t** create_G(int32_t n, int32_t bit_length, int32_t base_bit) {
    size_t rows = (n + 1) * bit_length;
    size_t cols = n + 1;
    int32_t** G = (int32_t**)malloc(rows * sizeof(int32_t*));
    for (size_t i = 0; i < rows; i++) {
        G[i] = (int32_t*)calloc(cols, sizeof(int32_t));
    }

    for (size_t i = 0; i < cols; i++) {
        for (size_t j = 0; j < bit_length; j++) {
            G[i * bit_length + j][i] = 1 << (32 - (j + 1) * base_bit);
        }
    }

    return G;
}


// GSW encryption
void GSW_encrypt(int32_t** C, int message, int32_t n, int32_t bit_length, int32_t Bksbit, int32_t* sk, int32_t* error_vec) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int32_t> bit_dist(0, 1);
    // uniform_int_distribution<int32_t> int32_t_dist(numeric_limits<int32_t>::min(), numeric_limits<int32_t>::max());
    uniform_int_distribution<int32_t> int32_t_dist(0, 5);

    size_t rows = (n + 1) * bit_length;
    size_t cols = n + 1;

    // create a matrix A
    for (size_t i = 0; i < rows; i++) {
        error_vec[i] = bit_dist(gen);
        for (size_t j = 1; j < n + 1; j++) {
            C[i][j] = int32_t_dist(gen);
        }
    }

    // creating the public key, i.e., matrix A = (s*A + e, A):
    for (size_t i = 0; i < rows; i++) {
        for (size_t k = 1; k < n + 1; k++) {
            C[i][0] += sk[k] * C[i][k];
        }
        C[i][0] += error_vec[i];
    }



    int32_t m = int32_t(message);
    if (m == 0) return;

    int32_t** G = create_G(n, bit_length, Bksbit);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            C[i][j] += G[i][j];
        }
        free(G[i]);
    }

    free(G);

    return;
}


// implement the external product
int32_t* external_product(int32_t** C, size_t C_rows, size_t C_cols, int32_t* lwe, int32_t bit_length, int32_t base_bit) {

    int32_t* G_inv_lwe = G_inverse_of_vector(lwe, C_cols, bit_length, base_bit);
    int32_t* external_product_result = vector_matrix_multiplication(G_inv_lwe, C, C_rows, C_cols);

    free(G_inv_lwe);
    return external_product_result;
}


// GSW implementation
void myGSW() {
    int32_t bit_length = 32;
    int32_t base_bit = 1;
    int32_t n = 560;
    cout << "Please wait......., some computations might take some time..." << endl;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int32_t> bit_dist(0, 1);
    uniform_int_distribution<int32_t> int32_t_dist(numeric_limits<int32_t>::min(), numeric_limits<int32_t>::max());
    // uniform_int_distribution<int32_t> int32_t_dist(0, 5);

    int32_t* sk = (int32_t*)malloc((n + 1) * sizeof(int32_t));
    sk[0] = 1;
    for (size_t i = 1; i < n + 1; i++) sk[i] = bit_dist(gen);

    // create an error vector
    int32_t* error_vec = (int32_t*)malloc((n + 1) * bit_length * sizeof(int32_t));

    // // create a matrix A
    size_t rows = (n + 1) * bit_length;
    int32_t** C = (int32_t**)malloc(rows * sizeof(int32_t*));
    for (size_t i = 0; i < rows; i++) {
        C[i] = (int32_t*)calloc(n + 1, sizeof(int32_t));
    }


    int32_t* lwe = (int32_t*)calloc(n + 1, sizeof(int32_t));
    for (size_t i = 1; i < n + 1; i++) {
        lwe[i] = int32_t_dist(gen);
        lwe[0] += sk[i] * lwe[i];
    }
    lwe[0] += 8;



    int32_t* G_inv_lwe = G_inverse_of_vector(lwe, n + 1, bit_length, base_bit);
    GSW_encrypt(C, 1, n, bit_length, base_bit, sk, error_vec);
    int32_t* result = vector_matrix_multiplication(G_inv_lwe, C, rows, n + 1);


    for (size_t i = 1; i < n + 1; i++) sk[i] = -sk[i];
    cout << vector_multiplication(lwe, sk, n + 1) << endl;    
    cout << vector_multiplication(error_vec, G_inv_lwe, (n + 1) * bit_length) << endl;
    cout << vector_multiplication(result, sk , n + 1) << endl;

    int32_t* c = external_product(C, rows, n + 1, lwe, bit_length, base_bit);
    cout << vector_multiplication(c, sk, n + 1) << endl;


    // int32_t q[3] = {123123123, 11115555, 18181818};
    // int32_t* ginv = G_inverse_of_vector(q, 2, bit_length, base_bit);

    // int32_t** g = create_G(2, bit_length, base_bit);
    // int32_t* d = vector_matrix_multiplication(ginv, g, 3*bit_length, 3);

    // cout << endl;
    // for (size_t i = 0; i < 3; i++) {
    //     cout << d[i] << " ";
    // }
    // cout << endl;

    free(sk);
    free(lwe);
    free(error_vec);
    for (size_t i = 0; i < rows; i++) {
        free(C[i]);
    }
    free(C);
    free(G_inv_lwe);
    free(result);
    // free(d);
    // free(c);
}



int main() {
        
    cout << "################################################################" << endl;
    cout << "#-------------  GSW Implementation ! -------------#" << endl;
    cout << "################################################################" << endl;

    auto start = clock();
    myGSW();
    auto end = clock();
    double time = ((double) end - start)/CLOCKS_PER_SEC;
    cout << "The program took " << time << " s to run." << endl;

    cout << endl;
    cout << "################################################################" << endl;
    cout << "#----------------------- End of Program -----------------------#" << endl;
    cout << "################################################################" << endl;

    return 0;
}