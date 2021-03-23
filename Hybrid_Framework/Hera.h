#ifndef __HERA_H__
#define __HERA_H__

#include <cstring>
#include <vector>
#include "ShakeAVX2.h"
#include "parms.h"

using namespace std;

typedef uint64_t* block_t;

block_t block_init(size_t sz);


/*
Provides functionality of cipher HERA. It takes nonce and counter as
inputs, and the nonce is assumed not to be reused. This class only
offers encryption part on the client side of CKKS-FV transciphering
framework.
 */
class Hera
{
    public:
        /*
         Create a HERA instance initialized with secret key

         @param[in] key The secret key
         */
        Hera(uint64_t *key)
        {
            key_ = block_init(BLOCKSIZE);
#if MONTGOMERY
            for (size_t i = 0; i < BLOCKSIZE; i++)
            {
                key_[i] = (key[i] << MONT_MOD_BIT) % MODULUS;
            }
#else
            for (int i = 0; i < BLOCKSIZE; i++)
            {
                key_[i] = key[i] % MODULUS;
            }
#endif
            rand_vectors_ = block_init(XOF_ELEMENT_COUNT + 8);
            round_keys_ = block_init(XOF_ELEMENT_COUNT);
            is_shake_init_ = false;
        }

        // Destruct a HERA instance
        ~Hera()
        {
            free(key_);
            free(rand_vectors_);
            if (is_shake_init_)
            {
                delete shake_;
            }
        }

        // Re-keying function
        void set_key(uint64_t *key);

        /*
        Both init and update function compute round keys from
        extendable output function. The difference is, the init
        function creates a new ShakeAVX2 object while the update
        function does not.

        @param[in] nonce Distinct nonce
        @param[in] counter Counter, but may be used as an integrated nonce
         */
        void init(uint64_t nonce, uint64_t counter);
        void update(uint64_t nonce, uint64_t counter);

        /*
        Both crypt and crypt_naive compute a block of HERA. crypt
        function is an AVX2 implementation while crypt_naive function
        is a reference code of HERA.

        @param[out] out Keystream of HERA
         */
        void crypt(block_t out);
        void crypt_naive(block_t out);

        // Copying outputs of XOF
        void get_rand_vectors(uint64_t *output);

    private:
        // Secret key
        block_t key_;

        // Round constants
        block_t rand_vectors_;

        // Key multiplied by round constant
        block_t round_keys_;

        // Shake object
        ShakeAVX2 *shake_;

        bool is_shake_init_;

        // The inner key schedule function in init and update
        void keyschedule();
};

#endif
