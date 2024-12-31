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


//////////////////////////////////////////////////////////////////////////////
// forward declaration:
vector<int32_t> G_inverse(int32_t num, long bit_length);

//////////////////////////////////////////////////////////////////////////////


void print_matrix(vector<vector<int32_t>>& matrix) {
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[0].size(); j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// vector<vector<int32_t>> G_inverse_of_matrix(vector<vector<int32_t>> ciphertext, long n, long bit_length) {
//     // this function  gets a ciphertext as input, which is a matrix of dimension (2 ⨯ 2*bit_length).
//     // The function returns the bitdecomposition of the ciphertext.
//     // the output is a matrix of dimension (2*bit_length ⨯ 2*bit_length) with binary entries.

//     vector<vector<int32_t>> bitdecomposed_ciphertext((n + 1) * bit_length, vector<int32_t>((n + 1) * bit_length, 0));
//     vector<int32_t> num_bitdecomposed(bit_length);

//     int32_t temp_num;
//     for (size_t i = 0; i < bitdecomposed_ciphertext.size(); i++) {
//         for (size_t j = 0; j < n + 1; j++) {
//             temp_num = ciphertext[j][i];
//             num_bitdecomposed = G_inverse(temp_num, bit_length);
//             for (size_t k = 0; k < bit_length; k++) {
//                 bitdecomposed_ciphertext[j * bit_length + k][i] = num_bitdecomposed[k];
//             }
//         }
//     }

//     return bitdecomposed_ciphertext;
// }

vector<int32_t> G_inverse_of_scalar(int32_t num, long bit_length) {
    // this function gets an integer number and bit-decomposes the number the way that G^-1 works:

    vector<int32_t> num_bitdecomposed(bit_length);
    for (size_t i = 0; i < bit_length; i++) {
        num_bitdecomposed[i] = (num >> i) & 1;
    }

    return num_bitdecomposed;
}


vector<int32_t> G_inverse_of_vector(vector<int32_t>& ciphertext, long bit_length) {
    // this function  gets a ciphertext as input, which is a matrix of dimension (2 ⨯ 2*bit_length).
    // The function returns the bitdecomposition of the ciphertext.
    // the output is a matrix of dimension (2*bit_length ⨯ 2*bit_length) with binary entries.

    vector<int32_t> bitdecomposed_ciphertext, num_bitdecomposed(bit_length);

    int32_t temp_num;
    for (size_t j = 0; j < ciphertext.size(); j++) {
        temp_num = ciphertext[j];
        num_bitdecomposed = G_inverse_of_scalar(temp_num, bit_length);
        bitdecomposed_ciphertext.insert(bitdecomposed_ciphertext.end(), num_bitdecomposed.begin(), num_bitdecomposed.end());
    }

    return bitdecomposed_ciphertext;
}

vector<int32_t> G_Ginv_multiplication(vector<vector<int32_t>>& G, vector<int32_t>& G_inv, long bit_length) {
    vector<int32_t> result(G[0].size(), 0);
    size_t j = 0;
    for (size_t i = 0; i < G.size(); i++) {
        result[j] += G[i][j] * G_inv[i];
        if ((i != 0) && (i % bit_length == 0)) j++;
    }
    return result;
}


int32_t vector_multiplication(vector<int32_t>& v1, vector<int32_t>& v2) {
    // this function gets two vectors and returns the dot product of the two vectors.
    if (v1.size() != v2.size()) {
        cout << "The lengths of vectors are different!" << endl;
        cout << "v1.length(): " << v1.size() << "\tv2.length(): " << v2.size() << endl;
    }

    int32_t result = int32_t(0);
    for (size_t i = 0; i < v1.size(); i++) result += v1[i] * v2[i];

    return result;
}


vector<int32_t> matrix_vector_multiplication(vector<vector<int32_t>>& mat, vector<int32_t>& vec) {
    // this function gets a matrix and a vector and returns the result of the matrix-vector multiplication.
    vector<int32_t> result(mat.size(), 0);
    if (mat[0].size() != vec.size()) {
        cout << "Error in the matrix_vector_multiplication function: the dimensions of the matrix and the vector do not match!" << endl;
    }
    for (size_t i = 0; i < mat.size(); i++) {
        for (size_t j = 0; j < mat[0].size(); j++) {
            result[i] += mat[i][j] * vec[j];
        }
    }
    return result;
}

vector<int32_t> vector_matrix_multiplication(vector<int32_t>& vec, vector<vector<int32_t>>& mat) {
    vector<int32_t> result(mat[0].size(), 0);
    for (size_t i = 0; i < mat[0].size(); i++) {
        for (size_t j = 0; j < vec.size(); j++) {
            result[i] += vec[j] * mat[j][i];
        }
    }
    return result;
}



// Creating a gaget matrix G:
vector<vector<int32_t>> create_G(long n, long bit_length) {
    vector<vector<int32_t>> G((n + 1) * bit_length, vector<int32_t>(n + 1, 0));

    for (size_t i = 0; i < n + 1; i++) {
        int32_t temp = int32_t(1);
        for (size_t j = 0; j < bit_length; j++) {
            G[i * bit_length + j][i] = temp;
            temp *= int32_t(2);
        }
    }

    return G;
}

// Implementation of the encryption algorithm for the GSW scheme
vector<vector<int32_t>> encrypt(int message, vector<vector<int32_t>>& A, long n, long bit_length) {
    
    int32_t m = int32_t(message);

    if (m == 0) return A;

    // // calculating A * R:
    // vector<vector<int32_t>> A_milttipliedBy_R(n, vector<int32_t>(n * bit_length, 0));
    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < n * bit_length; j++) {
    //         for (int k = 0; k < n * bit_length; k++) {
    //             A_milttipliedBy_R[i][j] += A[i][k] * R[k][j];
    //         }
    //     }
    // }

    // creating matrix G:
    vector<vector<int32_t>> G = create_G(n, bit_length);


    // Encrypting a message (where m can be 0 or 1) using the GSW scheme: ciphertext = A + m*G:
    vector<vector<int32_t>> ciphertext((n + 1) * bit_length, vector<int32_t>(n + 1, 0));
    for (size_t i = 0; i < (n + 1) * bit_length; i++) {
        for (size_t j = 0; j < n + 1; j++) {
            ciphertext[i][j] = A[i][j] + G[i][j];
        }
    }

    return ciphertext;
}

// Implementation of the decryption algorithm for the RGSW scheme
// int decrypt(vector<int32_t> sk, vector<vector<int32_t>> ciphertext, long n, long bit_length, int32_t modulous) {
//     int decrypted_message;

//     decrypted_message = 0; // initializing the variabele.
//     vector<int32_t> penultimate_column_of_c, last_column_of_c;
//     vector<int32_t> first_column_of_c, second_column_of_c;

//     for (size_t i = 0; i < n; i++) {
//         first_column_of_c.push_back(ciphertext[i][0]);
//         second_column_of_c.push_back(ciphertext[i][1]);

//         penultimate_column_of_c.push_back(ciphertext[i][2 * bit_length - 2]);
//         last_column_of_c.push_back(ciphertext[i][2 * bit_length - 1]);
//     }

//     int32_t result1, result2, result3, result4;
//     int32_t result3int32_t;
//     long comparison_result;

//     result1 = dot_product(sk, first_column_of_c);
//     result2 = dot_product(sk, second_column_of_c);
//     result3 = dot_product(sk, penultimate_column_of_c);
//     result4 = dot_product(sk, last_column_of_c);

//     int32_t p_dividedBy_2;;
//     div(p_dividedBy_2, modulous, int32_t(2));
//     conv(result3int32_t, result3);
//     comparison_result = compare(result3int32_t, p_dividedBy_2);

//     if (comparison_result > 0)
//         decrypted_message = 1;
//     else
//         decrypted_message = 0;

//     return decrypted_message;
// }


// Implementation of the addition gate for the GSW scheme
vector<vector<int32_t>> add(vector<vector<int32_t>>& c1, vector<vector<int32_t>>& c2) {
    vector<vector<int32_t>> addtion_result(c1.size(), vector<int32_t>(c1[0].size()));

    if ((c1.size() != c2.size()) || (c1[0].size() != c2[0].size())) {
        cout << "Error in the add function: the dimensions of the ciphertext do not matach!" << endl;
    }

    for (size_t i = 0; i < c1.size(); i++) {
        for (size_t j = 0; j < c1[0].size(); j++) {
            addtion_result[i][j] = c1[i][j] + c2[i][j];
        }
    }

    return addtion_result;
}


// Implementation of the external product gate for the GSW scheme
vector<int32_t> external_product(vector<vector<int32_t>>& C, vector<int32_t>& lwe) {
    vector<int32_t> external_product_result(C[0].size());

    if (C.size() != lwe.size()) {
        cout << "Error in the external_product function: the dimensions of the ciphertext do not matach!" << endl;
    }

    for (size_t i = 0; i < C.size(); i++) {
        external_product_result[i] = int32_t(0);
        for (size_t j = 0; j < C[0].size(); j++) {
            external_product_result[i] += C[i][j] * lwe[j];
        }
    }

    return external_product_result;
}


//////////////////////////////////////////////////////////////////////////////
// Implementation of the multiplication gate for the GSW scheme
// vector<vector<int32_t>> multiply(vector<vector<int32_t>> c1, vector<vector<int32_t>> c2, long n, long bit_length) {

//     if ((c1.size() != c2.size()) || (c1[0].size() != c2[0].size())) {
//         cout << "Error in the multiply function: the dimensions of the ciphertext do not matach!" << endl;
//     }

//     vector<vector<int32_t>> multiplication_result(c1.size(), vector<int32_t>(c2[0].size()));

//     // calculating G_inverse (c2):
//     vector<vector<int32_t>> G_inverseOf_c2 = G_inverse_of_ciphertext(c2, n, bit_length);

//     // initializing all elements of multiplication_result to zero.
//     for (size_t i = 0; i < c1.size(); i++) {
//         for (size_t j = 0; j < c2[0].size(); j++) {
//             multiplication_result[i][j] = int32_t(0);
//         }
//     }

//     // calculating c1 . G_inverse(c2):
//     int32_t temp;
//     for (size_t i = 0; i < c1.size(); i++) {
//         for (size_t j = 0; j < c1[0].size(); j++) {
//             for (size_t k = 0; k < G_inverseOf_c2.size(); k++) {
//                 //cout << "i: " << i << ", j: " << j << ", k: " << k << endl;
//                 temp = c1[i][k] * G_inverseOf_c2[k][j];
//                 multiplication_result[i][j] += temp;
//             }
//         }
//     }

//     return multiplication_result;
// }


void myGSW() {
    long bit_length = 16;
    long n = 560; 
    cout << "Please wait......., some computations might take some time...";

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int32_t> bit_dist(0, 1);
    // uniform_int_distribution<int32_t> int32_t_dist(numeric_limits<int32_t>::min(), numeric_limits<int32_t>::max());
    uniform_int_distribution<int32_t> int32_t_dist(0, 30);

    // generating the secret key:
    vector<int32_t> sk(n+1);
    sk[0] = int32_t(1);
    for (size_t i = 1; i < n+1; i++) sk[i] = bit_dist(gen);
    

    // generating an error vector:
    vector<int32_t> error_vec((n + 1) * bit_length);
    for (size_t i = 0; i < error_vec.size(); i++) error_vec[i] = bit_dist(gen);

    // creating the public key, i.e., matrix A = (s*A + e, A):
    vector<vector<int32_t>> A((n + 1) * bit_length, vector<int32_t>(n + 1, 0));

    for (int i = 0; i < A.size(); i++) 
        for (int j = 1; j < A[0].size(); j++) 
            A[i][j] = int32_t_dist(gen);

    for (size_t i = 0; i < A.size(); i++) {
        for (size_t k = 1; k < A[0].size(); k++) {
            A[i][0] += sk[k] * A[i][k];
        }
        A[i][0] += error_vec[i];
    }

    vector<int32_t> lwe(n+1, 0);
    for (int i = 1; i < n + 1; i++) {
        lwe[i] = int32_t_dist(gen);
        lwe[0] += sk[i] * lwe[i];
    }
    // lwe[0] += bit_dist(gen);
    lwe[0] += 8;

    vector<int32_t> G_inv = G_inverse_of_vector(lwe, bit_length);
    vector<vector<int32_t>> C = encrypt(1, A, n, bit_length);
    vector<int32_t> result = vector_matrix_multiplication(G_inv, C);

    for (size_t i = 1; i < sk.size(); i++) sk[i] = -sk[i];

    // cout << endl;
    // for (int i = 0; i < lwe.size(); i++) cout << lwe[i] << " ";
    // cout << endl;
    // vector<int32_t> lwe2 = vector_matrix_multiplication(G_inv, create_G(n, bit_length));
    // for (int i = 0; i < lwe2.size(); i++) cout << lwe2[i] << " ";
    // cout << endl;

    // bool ok = true;
    // for (int i = 0; i < lwe.size(); i++) {
    //     if (lwe[i] != lwe2[i]) {
    //         ok = false;
    //         break;
    //     }
    // } 
    // if (ok) cout << "The two vectors are equal!" << endl;
    // else cout << "The two vectors are not equal!" << endl;

    cout << endl;
    cout << vector_multiplication(sk, result) << endl;

    cout << vector_multiplication(G_inv, error_vec) << endl;
    cout << vector_multiplication(sk, lwe) << endl;

    // print_matrix(A);

//     // constructing the random matrix R (as defined in the GSW scheme) with binary entries:
//     // setting the modulous to 2, for generating a small random numbers:
//     vector<vector<int32_t>> R(n * bit_length, vector<int32_t>(n * bit_length));
//     int32_t temp_num;
//     for (size_t i = 0; i < n * bit_length; i++) {
//         for (size_t j = 0; j < n * bit_length; j++) {
//             R[i][j] = dist(gen) % 2;
//         }
//     }

    // creating two ciphertexts and testing the add and multiply operations on ciphertexts:
    // int message1 = 1, message2 = 0;
    // vector<vector<int32_t>> ciphertext1, ciphertext2;
    // vector<int32_t> ciphertext3, external_product_result;
    // vector<vector<int32_t>> addition_result_ct, multiplication_result_ct;

    // ciphertext1 = encrypt(message1, A, n, bit_length);
    // ciphertext2 = encrypt(message2, A, n, bit_length);

    //cout << "ciphertext1: " << endl;
    //cout << ciphertext1 << endl;

    //cout << "ciphertext2: " << endl;
    //cout << ciphertext2 << endl;

    // addition_result_ct = add(ciphertext1, ciphertext2);


    // external_product_result = external_product(ciphertext1, ciphertext3);



//     multiplication_result_ct = multiply(ciphertext1, ciphertext2, n, bit_length);

//     //cout << "\naddition_result: " << endl;
//     //cout << addition_result << endl;

//     //cout << "multiplication_result: " << endl;
//     //cout << multiplication_result << endl;

//     //////////////////////////////////////////////////////////////////////////////
//     // decrypting the result:
//     int decryption_result1, decryption_result2;

//     decryption_result1 = decrypt(sk, ciphertext1, n, bit_length, p);
//     decryption_result2 = decrypt(sk, ciphertext2, n, bit_length, p);

//     cout << "\n\nmessage 1 is: " << message1 << endl;
//     cout << "message 2 is: " << message2 << endl;
//     cout << "decryption_result1 is: " << decryption_result1 << endl;
//     cout << "decryption_result2 is: " << decryption_result2 << endl;

// }
}

int main() {
    
    cout << "################################################################" << endl;
    cout << "#------------- Welcome to GSW Implementation! -------------#" << endl;
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