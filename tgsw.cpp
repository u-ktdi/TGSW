/*
This is an implementation of the GSW scheme (security bit : At least 110 bit)

Articles and references:
GSW scheme: https://eprint.iacr.org/2013/340.pdf
*/

#include <iostream>
#include <random>
#include <cstdint>
#include <chrono>
#include <cstdlib>

using namespace std;
using Torus32 = int32_t;

const double stdev = 3.05e-5;

// random number generator
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<int32_t> bit_dist(0, 1);
uniform_int_distribution<int32_t> int32_t_dist(numeric_limits<int32_t>::min(), numeric_limits<int32_t>::max());
normal_distribution<double> error_distribution(0., stdev);



Torus32 modSwitchToTorus32(int32_t mu, int32_t Msize){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2; 
    uint64_t phase64 = mu*interv;
    return phase64>>32;
}


// from double to Torus32
Torus32 dtot32(double d) {
    return int32_t(int64_t((d - int64_t(d))*(INT64_C(1) << 32)));
}


Torus32 gaussian32(int32_t message){
    Torus32 _1s8 = modSwitchToTorus32(1, 8);
    Torus32 mu = message ? _1s8 : -_1s8;

    double err = error_distribution(gen);
    return mu + dtot32(err);
}


// b = <a_i,s_i> + m + e
// m comes with a scaling factor 
void LWE_encrypt(Torus32* result, int message, int32_t n, int32_t* sk) {
    
    result[0] = gaussian32(message); 

    for (int j = 1; j < n + 1; ++j) {
        result[j] = int32_t_dist(gen);
        result[0] += result[j]*sk[j];
    } 
    
    return;
}


// secret key generation
void sk_gen(int32_t* sk, int32_t n) {
    for (int32_t i = 1; i < n + 1; i++) {
        sk[i] = bit_dist(gen);
    }
    sk[0] = 1;
    return;
}


// print matrix
void print_matrix(Torus32** matrix, int32_t rows, int32_t cols) {
    for (int32_t i = 0; i < rows; i++) {
        for (int32_t j = 0; j < cols; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// create a matrix G^-1(scalar)
void G_inverse_of_scalar(Torus32* result, int32_t num, int32_t bit_length, int32_t base_bit) {
    const int32_t base = 1 << base_bit;
    const int32_t prec_offset = 1 << (32 - (1 + base_bit * bit_length)); // precision
    const int32_t mask = base - 1;

    num += prec_offset;

    for (int32_t i = 0; i < bit_length; i++) {
        result[i] = (num >> (32 - (i + 1) * base_bit)) & mask;
    }

    return;
}


// create a matrix G^-1(vector)
void G_inverse_of_vector(Torus32* result, Torus32* ciphertext, int32_t n, int32_t bit_length, int32_t base_bit) {
    int32_t index = 0;

    for (int32_t j = 0; j < n + 1; j++) {
        Torus32* num_bitdecomposed = (Torus32*)malloc(bit_length * sizeof(Torus32));
        G_inverse_of_scalar(num_bitdecomposed, ciphertext[j], bit_length, base_bit);
        for (int32_t k = 0; k < bit_length; k++) {
            result[index++] = num_bitdecomposed[k];
        }
        free(num_bitdecomposed);
    }

    return;
}


// G * G^-1 (calloc)
void G_Ginv_multiplication(Torus32* result, Torus32** G, Torus32* G_inv, int32_t n, int32_t bit_length) {
    int32_t cols = n + 1;

    int32_t j = 0, rows = (n + 1) * bit_length;
    for (int32_t i = 0; i < rows; i++) {
        result[j] += G[i][j] * G_inv[i];
        if ((i != 0) && (i % bit_length == 0)) j++;
    }

    return;
}


// internal product
int32_t vector_multiplication(Torus32* v1, Torus32* v2, int32_t size) {
    int32_t result = 0;

    for (int32_t i = 0; i < size; i++) {
        result += v1[i] * v2[i];
    }

    return result;
}

// matrix vector multiplication (need to prepare the vector with calloc)
void matrix_vector_multiplication(Torus32* result, Torus32** mat, Torus32* vec, int32_t rows, int32_t cols) {
    for (int32_t i = 0; i < rows; i++) {
        for (int32_t j = 0; j < cols; j++) {
            result[i] += mat[i][j] * vec[j];
        }
    }
    return;
}


// vector matrix multiplication
void vector_matrix_multiplication(Torus32* result, Torus32* vec, Torus32** mat, int32_t mat_rows, int32_t mat_cols) {

    for (int32_t i = 0; i < mat_cols; i++) {
        for (int32_t j = 0; j < mat_rows; j++) {
            result[i] += vec[j] * mat[j][i];
        }
    }

    return;
}


// create matrix G
void create_G(Torus32** G, int32_t n, int32_t bit_length, int32_t base_bit) {
    int32_t rows = (n + 1) * bit_length;
    int32_t cols = n + 1;

    for (int32_t i = 0; i < cols; i++) {
        for (int32_t j = 0; j < bit_length; j++) {
            G[i * bit_length + j][i] = 1 << (32 - (j + 1) * base_bit);
        }
    }

    return;
}


// GSW encryption
void GSW_encrypt(Torus32** C, int message, int32_t n, int32_t bit_length, int32_t Bksbit, int32_t* sk) {

    int32_t rows = (n + 1) * bit_length;
    int32_t cols = n + 1;

    // create a error vector
    Torus32* error_vec = (Torus32*)malloc(rows * sizeof(Torus32));

    // create a matrix A
    for (int32_t i = 0; i < rows; i++) {
        error_vec[i] = dtot32(error_distribution(gen));
        for (int32_t j = 1; j < n + 1; j++) {
            C[i][j] = int32_t_dist(gen);
        }
    }

    // creating the key, i.e., matrix A = (s*A + e, A):
    for (int32_t i = 0; i < rows; i++) {
        for (int32_t k = 1; k < n + 1; k++) {
            C[i][0] += sk[k] * C[i][k];
        }
        C[i][0] += error_vec[i];
    }

    int32_t m = int32_t(message);
    if (m == 0) return;


    // create the matrix G
    int32_t** G = (int32_t**)malloc(rows * sizeof(int32_t*));
    for (int32_t i = 0; i < rows; i++) {
        G[i] = (int32_t*)calloc(cols, sizeof(int32_t));
    }

    // C += G
    create_G(G, n, bit_length, Bksbit);
    for (int32_t i = 0; i < rows; i++) {
        for (int32_t j = 0; j < cols; j++) {
            C[i][j] += G[i][j];
        }
        free(G[i]);
    }

    free(G);
    free(error_vec);

    return;
}


// Add computation by tgsw (add two matrices)
void add(Torus32** C1, Torus32** C2, int32_t rows, int32_t cols) {
    for (int32_t i = 0; i < rows; i++) {
        for (int32_t j = 0; j < cols; j++) {
            C1[i][j] += C2[i][j];
        }
    }
    return;
}


// implement the external product
void external_product(Torus32* result, Torus32** C, int32_t C_rows, int32_t C_cols, Torus32* lwe, int32_t bit_length, int32_t base_bit) {
    Torus32* G_inv_lwe = (Torus32*)malloc(C_rows * sizeof(Torus32));
    G_inverse_of_vector(G_inv_lwe, lwe, C_cols - 1, bit_length, base_bit);
    vector_matrix_multiplication(result, G_inv_lwe, C, C_rows, C_cols);

    free(G_inv_lwe);
    return;
}


Torus32 lwePhase(const Torus32* sample, const int32_t* sk, int32_t n) {
    Torus32 axs = 0;

    for (int j = 1; j < n + 1; ++j) {
        axs += sample[j]*sk[j];
    } 
    return sample[0] - axs;
}


int32_t LWE_Decrypt(const Torus32 *sample, const int32_t* sk, int32_t n) {

    Torus32 mu = lwePhase(sample, sk, n);
    return (mu > 0 ? 1 : 0); //we have to do that because of the C binding
}

Torus32* AND_GSW(Torus32* c1, Torus32* c2, int32_t n, int32_t bit_length, int32_t base_bit) {

    Torus32* c = (Torus32*)calloc(n + 1, sizeof(Torus32));
    for (int32_t i = 0; i < n + 1; i++) {
        c[i] = c1[i] * c2[i];
    }
    return c;
} 


// GSW implementation
void GSW() {
    int32_t bit_length = 8;
    int32_t base_bit = 2;
    int32_t n = 560;
    int32_t message = 1;

    cout << "Please wait......., some computations might take some time..." << endl;

    // create a secret key
    int32_t* sk = (int32_t*)malloc((n + 1) * sizeof(int32_t));

    // // create a matrix A
    int32_t rows = (n + 1) * bit_length;
    Torus32** C = (Torus32**)malloc(rows * sizeof(Torus32*));
    for (int32_t i = 0; i < rows; i++) {
        C[i] = (Torus32*)calloc(n + 1, sizeof(Torus32));
    }

    // create a lwe
    Torus32* lwe = (Torus32*)calloc(n + 1, sizeof(Torus32));

    // create the G^-1(lwe)
    Torus32* G_inv_lwe = (Torus32*)malloc((n + 1) * bit_length * sizeof(Torus32));

    sk_gen(sk, n);
    LWE_encrypt(lwe, message, n, sk);
    G_inverse_of_vector(G_inv_lwe, lwe, n, bit_length, base_bit);
    GSW_encrypt(C, 1, n, bit_length, base_bit, sk);

    Torus32* result = (Torus32*)calloc(n + 1, sizeof(Torus32));
    vector_matrix_multiplication(result, G_inv_lwe, C, rows, n + 1);


    for (int32_t i = 1; i < n + 1; i++) sk[i] = -sk[i];
    // cout << vector_multiplication(lwe, sk, n + 1) << endl;    
    // cout << vector_multiplication(error_vec, G_inv_lwe, (n + 1) * bit_length) << endl;
    // cout << vector_multiplication(result, sk , n + 1) << endl;
    
    Torus32* c = (Torus32*)calloc(n + 1, sizeof(Torus32));
    external_product(c, C, rows, n + 1, lwe, bit_length, base_bit);
    // cout << vector_multiplication(c, sk, n + 1) << endl;

    cout << "The message is: " << LWE_Decrypt(c, sk, n) << endl;

    free(sk);
    free(lwe);
    for (int32_t i = 0; i < rows; i++) {
        free(C[i]);
    }
    free(C);
    free(G_inv_lwe);
    free(result);
    // free(d);
    free(c);
}



int main() {
        
    cout << "#--------------------  GSW Implementation ! -------------------#" << endl;
    cout << endl;

    auto start = clock();
    GSW();
    auto end = clock();
    double time = ((double) end - start)/CLOCKS_PER_SEC;
    cout << "The program took " << time << " s to run." << endl;

    cout << endl;
    cout << "#----------------------- End of Program -----------------------#" << endl;

    return 0;
}