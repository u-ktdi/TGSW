#include "tgsw.hpp"

// GSW implementation
int main() {
    int32_t bit_length = 8;
    int32_t base_bit = 2;
    int32_t n = 560;
    int32_t trial = 100;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int32_t> bit_dist(0, 1);

    std::cout << "Please wait......., some computations might take some time..." << std::endl;

    int32_t cnt = 0;
    for (int32_t i = 0; i < trial; i++) {
        // create a message
        int32_t lwe_message = bit_dist(gen);
        int32_t gsw_message = 1;

        // create a secret key
        int32_t* sk = (int32_t*)malloc((n + 1) * sizeof(int32_t));

        // create a matrix A
        int32_t rows = (n + 1) * bit_length;
        Torus32** C = (Torus32**)malloc(rows * sizeof(Torus32*));
        for (int32_t i = 0; i < rows; i++) {
            C[i] = (Torus32*)calloc(n + 1, sizeof(Torus32));
        }

        // create a lwe ciphertext
        Torus32* lwe = (Torus32*)calloc(n + 1, sizeof(Torus32));


        // key generation
        sk_gen(sk, n);

        // LWE encryption
        LWE_encrypt(lwe, lwe_message, n, sk);

        // GSW encryption
        GSW_encrypt(C, gsw_message, n, bit_length, base_bit, sk);

        // result = C *** lwe (*** is the external product)
        Torus32* result = (Torus32*)calloc(n + 1, sizeof(Torus32));
        external_product(result, C, lwe, n, bit_length, base_bit);

        // LWE decryption
        int32_t m = LWE_Decrypt(result, sk, n);  

        if (m != (lwe_message * gsw_message)) {
            std::cout << "####### trial : " << i << " #######" << std::endl;
            if (lwe_message != LWE_Decrypt(lwe, sk, n)) std::cout << "plaintext : " << lwe_message << " , decryption result : " << LWE_Decrypt(lwe, sk, n) << std::endl;
            std::cout << "plaintext : " << gsw_message*lwe_message << " , decryption result : " << m << std::endl << std::endl;
            
            cnt++;
        }
        
        // free memory
        free(sk);
        free(lwe);
        for (int32_t i = 0; i < rows; i++) {
            free(C[i]);
        }
        free(C);
        free(result);
    }
    std::cout << "Failure : " << cnt << " / " << trial << " ! " << std::endl;
}