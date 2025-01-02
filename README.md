# Homomorphic Encryption Implementation of TGSW (GSW) Scheme

Homomorphic encryption has gained significant attention in recent years as a groundbreaking cryptographic method.  
Various schemes such as TFHE, CKKS, GSW, BGV, and FV have been proposed.  
Among these, libraries for TFHE, CKKS, FV, and BGV are already available, while GSW remains unimplemented in any public library.  
To address this gap, I have developed an implementation of GSW, specifically focusing on TGSW (Torus GSW).  
I hope this code proves to be a valuable resource for you.

---

## How to Execute the Code

Since most widely-used homomorphic encryption libraries are implemented in C++, this project is also written in C++.  
Please ensure that your environment is set up with tools such as `g++` or `gcc` to compile and execute C++ code.  
Once your environment is ready, follow the steps below to run the code:

```bash
g++ -o tgsw tgsw.cpp
./tgsw
```

> [!NOTE]
> Multiplication between two GSW ciphertexts is not yet implemented.
This feature is planned for a future update.

