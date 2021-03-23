#include <iostream>
#include <iomanip>
#include <random>
#include "Hera.h"

int main()
{
    int elem_length = MOD_BIT_COUNT / 4;
    if (MOD_BIT_COUNT % 4)
        elem_length += 1;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> dis(0, MODULUS - 1);

    uint64_t key[16];
    for (int i = 0; i < 16; i++)
    {
        key[i] = dis(gen);
    }
    cout << "Key: " << flush;
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        cout << hex << setfill('0') << setw(elem_length) << key[i] << " ";
    }
    cout << dec << endl;

    uint64_t nonce = 0x123456789abcdef;
    uint64_t counter = 0;

    Hera cipher(key);
    cout << "Hera init" << endl;
    cipher.init(nonce, counter);

    // Check validity of rand vector
    cout << "Check validity" << endl;
    uint64_t buf[XOF_ELEMENT_COUNT];
    for (size_t i = 0; i < (1ULL << 16); i++)
    {
        cipher.update(nonce, i);
        cipher.get_rand_vectors(buf);
#if MONTGOMERY
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            uint64_t m = buf[i] & ALL_ONE_MOD_BIT;
            m = (m * MONT_PSEUDO_INV) & ALL_ONE_MOD_BIT;
            buf[i] = (buf[i] + m * MODULUS) >> MONT_MOD_BIT;
            if (buf[i] >= MODULUS)
                buf[i] -= MODULUS;
        }
#endif

        for (int j = 0; j < XOF_ELEMENT_COUNT; j++)
        {
            if (buf[j] >= MODULUS)
            {
                cout << "Counter_" << i << "[" << j << "]: "
                    << buf[j] << endl;
                return -1;
            }
        }

        block_t keystream_naive = block_init(BLOCKSIZE);
        cipher.crypt_naive(keystream_naive);

        block_t keystream_avx2 = block_init(BLOCKSIZE);
        cipher.crypt(keystream_avx2);

        for (int j = 0; j < BLOCKSIZE; j++)
        {
            if (keystream_naive[j] != keystream_avx2[j])
            {
                cout << "Counter " << hex << i << dec << endl;
                cout << "  - AVX2 : " << flush;
                for (int k = 0; k < BLOCKSIZE; k++)
                {
                    cout << setfill('0') << setw(elem_length) << hex << keystream_avx2[k] << " ";
                }
                cout << dec << endl;
                cout << "  - Naive: " << flush;
                for (int k = 0; k < BLOCKSIZE; k++)
                {
                    cout << setfill('0') << setw(elem_length) << hex << keystream_naive[k] << " ";
                }
                cout << dec << endl;

                break;
            }
        }
    }

    cout << "Get random vectors" << endl;
    uint64_t rand_vectors[XOF_ELEMENT_COUNT];
    cipher.get_rand_vectors(rand_vectors);
#if MONTGOMERY
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        uint64_t m = rand_vectors[i] & ALL_ONE_MOD_BIT;
        m = (m * MONT_PSEUDO_INV) & ALL_ONE_MOD_BIT;
        rand_vectors[i] = (rand_vectors[i] + m * MODULUS) >> MONT_MOD_BIT;
        if (rand_vectors[i] >= MODULUS)
            rand_vectors[i] -= MODULUS;
    }
#endif

    // Print random vectors
    for (int r = 0; r < ROUNDS; r++)
    {
        cout << "Round " << r << "  : ";
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            cout << hex << setfill('0') << setw(elem_length) << rand_vectors[r * BLOCKSIZE + i] << dec << " ";
        }
        cout << endl;
    }
    cout << "Final ARK: ";
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        cout << hex << setfill('0') << setw(elem_length) << rand_vectors[ROUNDS * BLOCKSIZE + i] << dec << " ";
    }
    cout << endl;

    // Print encryption result
    cout << "Result" << endl;
    block_t keystream_avx2 = block_init(BLOCKSIZE);
    cipher.crypt(keystream_avx2);
    cout << "  - AVX2 : ";
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        cout << hex << setfill('0') << setw(elem_length) << keystream_avx2[i] << dec << " ";
    }
    cout << endl;

    block_t keystream_naive = block_init(BLOCKSIZE);
    cipher.crypt_naive(keystream_naive);
    cout << "  - naive: ";
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        cout << hex << setfill('0') << setw(elem_length) << keystream_naive[i] << dec << " ";
    }
    cout << endl;

    return 0;
}

