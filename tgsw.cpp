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
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <iomanip>

using namespace std;

void print_matrix(vector<vector<int32_t>>& matrix) {
    for (size_t i = 0; i < matrix.size(); i++) {
        int32_t* row = matrix[i].data();
        for (size_t j = 0; j < matrix[0].size(); j++) {
            cout << row[j] << " ";
        }
        cout << endl;
    }
}

vector<int32_t> G_inverse_of_scalar(int32_t num, int32_t bit_length, int32_t Bksbit) {
    const int32_t base = 1 << Bksbit;
    const int32_t prec_offset = 1 << (32 - (1 + Bksbit * bit_length));
    const int32_t mask = base - 1;

    num += prec_offset;
    vector<int32_t> num_bitdecomposed(bit_length);
    int32_t* ptr = num_bitdecomposed.data();

    for (size_t i = 0; i < bit_length; i++) {
        ptr[i] = (num >> (32 - (i + 1) * Bksbit)) & mask;
    }

    return num_bitdecomposed;
}

vector<int32_t> G_inverse_of_vector(vector<int32_t>& ciphertext, int32_t bit_length, int32_t Bksbit) {
    vector<int32_t> bitdecomposed_ciphertext, num_bitdecomposed(bit_length);

    int32_t* cipher_ptr = ciphertext.data();
    size_t cipher_size = ciphertext.size();

    for (size_t j = 0; j < cipher_size; j++) {
        num_bitdecomposed = G_inverse_of_scalar(cipher_ptr[j], bit_length, Bksbit);
        bitdecomposed_ciphertext.insert(bitdecomposed_ciphertext.end(), num_bitdecomposed.begin(), num_bitdecomposed.end());
    }

    return bitdecomposed_ciphertext;
}

vector<int32_t> G_Ginv_multiplication(vector<vector<int32_t>>& G, vector<int32_t>& G_inv, int32_t bit_length) {
    vector<int32_t> result(G[0].size(), 0);
    int32_t* result_ptr = result.data();

    size_t j = 0;
    for (size_t i = 0; i < G.size(); i++) {
        result_ptr[j] += G[i][j] * G_inv[i];
        if ((i != 0) && (i % bit_length == 0)) j++;
    }

    return result;
}

int32_t vector_multiplication(vector<int32_t>& v1, vector<int32_t>& v2) {
    if (v1.size() != v2.size()) {
        cout << "The lengths of vectors are different!" << endl;
    }

    int32_t result = 0;
    int32_t* v1_ptr = v1.data();
    int32_t* v2_ptr = v2.data();

    for (size_t i = 0; i < v1.size(); i++) {
        result += v1_ptr[i] * v2_ptr[i];
    }

    return result;
}

vector<int32_t> matrix_vector_multiplication(vector<vector<int32_t>>& mat, vector<int32_t>& vec) {
    vector<int32_t> result(mat.size(), 0);
    int32_t* res_ptr = result.data();
    int32_t* vec_ptr = vec.data();

    for (size_t i = 0; i < mat.size(); i++) {
        int32_t* mat_ptr = mat[i].data();
        for (size_t j = 0; j < mat[0].size(); j++) {
            res_ptr[i] += mat_ptr[j] * vec_ptr[j];
        }
    }

    return result;
}

vector<vector<int32_t>> create_G(int32_t n, int32_t bit_length, int32_t Bksbit) {
    vector<vector<int32_t>> G((n + 1) * bit_length, vector<int32_t>(n + 1, 0));

    for (size_t i = 0; i < n + 1; i++) {
        for (size_t j = 0; j < bit_length; j++) {
            G[i * bit_length + j][i] = 1 << (32 - (j + 1) * Bksbit);
        }
    }

    return G;
}

void myGSW() {
    int32_t bit_length = 31;
    int32_t Bksbit = 1;
    int32_t n = 560;
    cout << "Please wait......., some computations might take some time..." << endl;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int32_t> bit_dist(0, 1);
    uniform_int_distribution<int32_t> int32_t_dist(numeric_limits<int32_t>::min(), numeric_limits<int32_t>::max());

    vector<int32_t> sk(n + 1);
    sk[0] = 1;
    int32_t* sk_ptr = sk.data();
    for (size_t i = 1; i < n + 1; i++) sk_ptr[i] = bit_dist(gen);

    vector<vector<int32_t>> A((n + 1) * bit_length, vector<int32_t>(n + 1, 0));
    for (int i = 0; i < A.size(); i++) {
        int32_t* row = A[i].data();
        for (int j = 1; j < A[0].size(); j++) {
            row[j] = int32_t_dist(gen);
        }
    }

    vector<int32_t> lwe(n + 1, 0);
    int32_t* lwe_ptr = lwe.data();
    for (int i = 1; i < n + 1; i++) {
        lwe_ptr[i] = int32_t_dist(gen);
        lwe_ptr[0] += sk_ptr[i] * lwe_ptr[i];
    }

    lwe_ptr[0] += 8;

    vector<int32_t> G_inv = G_inverse_of_vector(lwe, bit_length, Bksbit);

    cout << "GSW Implementation completed." << endl;
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