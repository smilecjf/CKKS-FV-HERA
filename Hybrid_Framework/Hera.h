#ifndef __HERA_H__
#define __HERA_H__

#include <cstring>
#include <vector>
#include "ShakeAVX2.h"
#include "parms.h"

using namespace std;

typedef uint64_t* block_t;

block_t block_init(size_t sz);

class Hera
{
    public:
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

        ~Hera()
        {
            free(key_);
            free(rand_vectors_);
            if (is_shake_init_)
            {
                delete shake_;
            }
        }

        // Hera public functions
        void set_key(uint64_t *key);
        void init(uint64_t nonce, uint64_t counter);
        void update(uint64_t nonce, uint64_t counter);
        void crypt(block_t out);
        void crypt_naive(block_t out);
        void get_rand_vectors(uint64_t *output);
    private:
        // Hera private data
        block_t key_;
        block_t rand_vectors_;
        block_t round_keys_;
        ShakeAVX2 *shake_;
        bool is_shake_init_;

        // Hera private functions
        void keyschedule();
};

#endif
