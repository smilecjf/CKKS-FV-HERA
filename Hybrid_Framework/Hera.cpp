#include <stdlib.h>
#include <x86intrin.h>
#include "Hera.h"

#define MONT_REDUCTION (MOD_BIT_COUNT >= 27)

__attribute__((aligned(32))) uint64_t INPUT_CONSTANT[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
constexpr int unpack_order[8] = {2, 3, 6, 7, 0, 1, 4, 5};

block_t block_init(size_t sz)
{
    block_t block = static_cast<uint64_t*>
        (aligned_alloc(32, sizeof(uint64_t) * sz));
    memset(block, 0, sizeof(uint64_t) * sz);
    return block;
}

// Hera public functions
void Hera::set_key(uint64_t *key)
{
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
}

void Hera::init(uint64_t nonce, uint64_t counter)
{
    shake_ = new ShakeAVX2(nonce, counter);
    is_shake_init_ = true;
    keyschedule();
}

void Hera::update(uint64_t nonce, uint64_t counter)
{
    shake_->update(nonce, counter);
    keyschedule();
}

void Hera::crypt_naive(block_t out)
{
    uint64_t buf[BLOCKSIZE];

    for (int i = 0; i < 16; i++)
        out[i] = i + 1;

    for (int r = 0; r < ROUNDS; r++)
    {
        // AddRoundKey
        uint64_t round_key[BLOCKSIZE];
        memcpy(round_key, round_keys_ + r * BLOCKSIZE, sizeof(uint64_t) * BLOCKSIZE);

        for (int i = 0; i < BLOCKSIZE; i++)
        {
#if MONTGOMERY
            uint64_t m = round_key[i] & ALL_ONE_MOD_BIT;
            m = (m * MONT_PSEUDO_INV) & ALL_ONE_MOD_BIT;
            round_key[i] = (round_key[i] + m * MODULUS) >> MONT_MOD_BIT;
            if (round_key[i] >= MODULUS)
                round_key[i] -= MODULUS;
#endif
            out[i] += round_key[i];
            out[i] %= MODULUS;
        }

        // MixColumns
        for (int row = 0; row < 4; row++)
        {
            for (int col = 0; col < 4; col++)
            {
                buf[row * 4 + col] = 2 * out[row * 4 + col];
                buf[row * 4 + col] += 3 * out[((row + 1) % 4) * 4 + col];
                buf[row * 4 + col] += out[((row + 2) % 4) * 4 + col];
                buf[row * 4 + col] += out[((row + 3) % 4) * 4 + col];
                buf[row * 4 + col] %= MODULUS;
            }
        }

        // MixRows
        for (int row = 0; row < 4; row++)
        {
            for (int col = 0; col < 4; col++)
            {
                out[row * 4 + col] = 2 * buf[row * 4 + col];
                out[row * 4 + col] += 3 * buf[row * 4 + (col + 1) % 4];
                out[row * 4 + col] += buf[row * 4 + (col + 2) % 4];
                out[row * 4 + col] += buf[row * 4 + (col + 3) % 4];
                out[row * 4 + col] %= MODULUS;
            }
        }

        // Cube
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            buf[i] = (out[i] * out[i]) % MODULUS;
            out[i] = (out[i] * buf[i]) % MODULUS;
        }
    }

    // Finalize
    {
        // MixColumns
        for (int row = 0; row < 4; row++)
        {
            for (int col = 0; col < 4; col++)
            {
                buf[row * 4 + col] = 2 * out[row * 4 + col];
                buf[row * 4 + col] += 3 * out[((row + 1) % 4) * 4 + col];
                buf[row * 4 + col] += out[((row + 2) % 4) * 4 + col];
                buf[row * 4 + col] += out[((row + 3) % 4) * 4 + col];
                buf[row * 4 + col] %= MODULUS;
            }
        }

        // MixRows
        for (int row = 0; row < 4; row++)
        {
            for (int col = 0; col < 4; col++)
            {
                out[row * 4 + col] = 2 * buf[row * 4 + col];
                out[row * 4 + col] += 3 * buf[row * 4 + (col + 1) % 4];
                out[row * 4 + col] += buf[row * 4 + (col + 2) % 4];
                out[row * 4 + col] += buf[row * 4 + (col + 3) % 4];
                out[row * 4 + col] %= MODULUS;
            }
        }

        // AddRoundKey
        uint64_t round_key[BLOCKSIZE];
        memcpy(round_key, round_keys_ + ROUNDS * BLOCKSIZE, sizeof(uint64_t) * BLOCKSIZE);
        for (int i = 0; i < BLOCKSIZE; i++)
        {
#if MONTGOMERY
            uint64_t m = round_key[i] & ALL_ONE_MOD_BIT;
            m = (m * MONT_PSEUDO_INV) & ALL_ONE_MOD_BIT;
            round_key[i] = (round_key[i] + m * MODULUS) >> MONT_MOD_BIT;
            if (round_key[i] >= MODULUS)
                round_key[i] -= MODULUS;
#endif
            out[i] += round_key[i];
            out[i] %= MODULUS;
        }
    }
}

void Hera::crypt(block_t out)
{
#if MONTGOMERY
    // Montgomery form change
    for (size_t i = 0; i < BLOCKSIZE; i++)
    {
        out[i] = ((i + 1) << MONT_MOD_BIT) % MODULUS;
    }

    __m256i u0, u1, u2, u3, v0, v1, v2, v3, p0, p1, p2, p3, q0, q1, q2, q3;

    // Set constant for montgomery reduction
    q0 = _mm256_set1_epi64x(ALL_ONE_MOD_BIT);
    q1 = _mm256_set1_epi64x(MONT_PSEUDO_INV);
    q2 = _mm256_set1_epi64x(MODULUS);
    q3 = _mm256_set1_epi64x(MODULUS - 1);

    // Load
    u0 = _mm256_load_si256((__m256i *) out);
    u1 = _mm256_load_si256((__m256i *) (out + 4));
    u2 = _mm256_load_si256((__m256i *) (out + 8));
    u3 = _mm256_load_si256((__m256i *) (out + 12));


    // Round function
    for (size_t r = 0; r < ROUNDS; r++)
    {
        // AddRoundKey
        v0 = _mm256_load_si256((__m256i *) (round_keys_ + r * BLOCKSIZE));
        v1 = _mm256_load_si256((__m256i *) (round_keys_ + r * BLOCKSIZE + 4));
        v2 = _mm256_load_si256((__m256i *) (round_keys_ + r * BLOCKSIZE + 8));
        v3 = _mm256_load_si256((__m256i *) (round_keys_ + r * BLOCKSIZE + 12));

        u0 = _mm256_add_epi32(u0, v0);
        u1 = _mm256_add_epi32(u1, v1);
        u2 = _mm256_add_epi32(u2, v2);
        u3 = _mm256_add_epi32(u3, v3);

        // MixColumns
        p1 = _mm256_slli_epi32(u1, 1);
        p0 = _mm256_slli_epi32(u0, 1);
        p1 = _mm256_add_epi32(p1, u1);
        p3 = _mm256_add_epi32(u2, u3);
        p2 = _mm256_add_epi32(p1, p0);
        v0 = _mm256_add_epi32(p2, p3);

        p2 = _mm256_slli_epi32(u2, 1);
        p1 = _mm256_slli_epi32(u1, 1);
        p2 = _mm256_add_epi32(p2, u2);
        p3 = _mm256_add_epi32(u0, u3);
        p0 = _mm256_add_epi32(p2, p1);
        v1 = _mm256_add_epi32(p0, p3);

        p3 = _mm256_slli_epi32(u3, 1);
        p2 = _mm256_slli_epi32(u2, 1);
        p3 = _mm256_add_epi32(p3, u3);
        p0 = _mm256_add_epi32(u0, u1);
        p1 = _mm256_add_epi32(p2, p3);
        v2 = _mm256_add_epi32(p0, p1);

        p0 = _mm256_slli_epi32(u0, 1);
        p3 = _mm256_slli_epi32(u3, 1);
        p0 = _mm256_add_epi32(p0, u0);
        p1 = _mm256_add_epi32(u1, u2);
        p2 = _mm256_add_epi32(p0, p3);
        v3 = _mm256_add_epi32(p1, p2);

        // Matrix Transpose
        u0 = _mm256_unpacklo_epi64(v0, v1);
        u1 = _mm256_unpackhi_epi64(v0, v1);
        u2 = _mm256_unpacklo_epi64(v2, v3);
        u3 = _mm256_unpackhi_epi64(v2, v3);

        v0 = _mm256_permute2x128_si256(u0, u2, 0b00100000);
        v1 = _mm256_permute2x128_si256(u1, u3, 0b00100000);
        v2 = _mm256_permute2x128_si256(u0, u2, 0b00110001);
        v3 = _mm256_permute2x128_si256(u1, u3, 0b00110001);

#if MONT_REDUCTION
        p0 = _mm256_set1_epi64x(MONT_ONE);
        v0 = _mm256_mul_epu32(v0, p0);
        v1 = _mm256_mul_epu32(v1, p0);
        v2 = _mm256_mul_epu32(v2, p0);
        v3 = _mm256_mul_epu32(v3, p0);

        p0 = _mm256_and_si256(v0, q0);
        p1 = _mm256_and_si256(v1, q0);
        p2 = _mm256_and_si256(v2, q0);
        p3 = _mm256_and_si256(v3, q0);

        p0 = _mm256_mul_epu32(p0, q1);
        p1 = _mm256_mul_epu32(p1, q1);
        p2 = _mm256_mul_epu32(p2, q1);
        p3 = _mm256_mul_epu32(p3, q1);

        p0 = _mm256_and_si256(p0, q0);
        p1 = _mm256_and_si256(p1, q0);
        p2 = _mm256_and_si256(p2, q0);
        p3 = _mm256_and_si256(p3, q0);

        p0 = _mm256_mul_epu32(p0, q2);
        p1 = _mm256_mul_epu32(p1, q2);
        p2 = _mm256_mul_epu32(p2, q2);
        p3 = _mm256_mul_epu32(p3, q2);

        v0 = _mm256_add_epi64(p0, v0);
        v1 = _mm256_add_epi64(p1, v1);
        v2 = _mm256_add_epi64(p2, v2);
        v3 = _mm256_add_epi64(p3, v3);

        v0 = _mm256_srli_epi64(v0, MONT_MOD_BIT);
        v1 = _mm256_srli_epi64(v1, MONT_MOD_BIT);
        v2 = _mm256_srli_epi64(v2, MONT_MOD_BIT);
        v3 = _mm256_srli_epi64(v3, MONT_MOD_BIT);
#endif

        // MixRows
        p1 = _mm256_slli_epi32(v1, 1);
        p0 = _mm256_slli_epi32(v0, 1);
        p1 = _mm256_add_epi32(p1, v1);
        p3 = _mm256_add_epi32(v2, v3);
        p2 = _mm256_add_epi32(p1, p0);
        u0 = _mm256_add_epi32(p2, p3);

        p2 = _mm256_slli_epi32(v2, 1);
        p1 = _mm256_slli_epi32(v1, 1);
        p2 = _mm256_add_epi32(p2, v2);
        p3 = _mm256_add_epi32(v0, v3);
        p0 = _mm256_add_epi32(p2, p1);
        u1 = _mm256_add_epi32(p0, p3);

        p3 = _mm256_slli_epi32(v3, 1);
        p2 = _mm256_slli_epi32(v2, 1);
        p3 = _mm256_add_epi32(p3, v3);
        p0 = _mm256_add_epi32(v0, v1);
        p1 = _mm256_add_epi32(p2, p3);
        u2 = _mm256_add_epi32(p0, p1);

        p0 = _mm256_slli_epi32(v0, 1);
        p3 = _mm256_slli_epi32(v3, 1);
        p0 = _mm256_add_epi32(p0, v0);
        p1 = _mm256_add_epi32(v1, v2);
        p2 = _mm256_add_epi32(p0, p3);
        u3 = _mm256_add_epi32(p1, p2);

        // Matrix Transpose
        v0 = _mm256_unpacklo_epi64(u0, u1);
        v1 = _mm256_unpackhi_epi64(u0, u1);
        v2 = _mm256_unpacklo_epi64(u2, u3);
        v3 = _mm256_unpackhi_epi64(u2, u3);

        u0 = _mm256_permute2x128_si256(v0, v2, 0b00100000);
        u1 = _mm256_permute2x128_si256(v1, v3, 0b00100000);
        u2 = _mm256_permute2x128_si256(v0, v2, 0b00110001);
        u3 = _mm256_permute2x128_si256(v1, v3, 0b00110001);

#if MONT_REDUCTION
        p0 = _mm256_set1_epi64x(MONT_ONE);
        u0 = _mm256_mul_epu32(u0, p0);
        u1 = _mm256_mul_epu32(u1, p0);
        u2 = _mm256_mul_epu32(u2, p0);
        u3 = _mm256_mul_epu32(u3, p0);

        p0 = _mm256_and_si256(u0, q0);
        p1 = _mm256_and_si256(u1, q0);
        p2 = _mm256_and_si256(u2, q0);
        p3 = _mm256_and_si256(u3, q0);

        p0 = _mm256_mul_epu32(p0, q1);
        p1 = _mm256_mul_epu32(p1, q1);
        p2 = _mm256_mul_epu32(p2, q1);
        p3 = _mm256_mul_epu32(p3, q1);

        p0 = _mm256_and_si256(p0, q0);
        p1 = _mm256_and_si256(p1, q0);
        p2 = _mm256_and_si256(p2, q0);
        p3 = _mm256_and_si256(p3, q0);

        p0 = _mm256_mul_epu32(p0, q2);
        p1 = _mm256_mul_epu32(p1, q2);
        p2 = _mm256_mul_epu32(p2, q2);
        p3 = _mm256_mul_epu32(p3, q2);

        u0 = _mm256_add_epi64(p0, u0);
        u1 = _mm256_add_epi64(p1, u1);
        u2 = _mm256_add_epi64(p2, u2);
        u3 = _mm256_add_epi64(p3, u3);

        u0 = _mm256_srli_epi64(u0, MONT_MOD_BIT);
        u1 = _mm256_srli_epi64(u1, MONT_MOD_BIT);
        u2 = _mm256_srli_epi64(u2, MONT_MOD_BIT);
        u3 = _mm256_srli_epi64(u3, MONT_MOD_BIT);
#endif

        // Cube map
        // 1. Square
        v0 = _mm256_mul_epu32(u0, u0);
        v1 = _mm256_mul_epu32(u1, u1);
        v2 = _mm256_mul_epu32(u2, u2);
        v3 = _mm256_mul_epu32(u3, u3);

        // 2. Montgomery reduction
        p0 = _mm256_and_si256(v0, q0);
        p1 = _mm256_and_si256(v1, q0);
        p2 = _mm256_and_si256(v2, q0);
        p3 = _mm256_and_si256(v3, q0);

        p0 = _mm256_mul_epu32(p0, q1);
        p1 = _mm256_mul_epu32(p1, q1);
        p2 = _mm256_mul_epu32(p2, q1);
        p3 = _mm256_mul_epu32(p3, q1);

        p0 = _mm256_and_si256(p0, q0);
        p1 = _mm256_and_si256(p1, q0);
        p2 = _mm256_and_si256(p2, q0);
        p3 = _mm256_and_si256(p3, q0);

        p0 = _mm256_mul_epu32(p0, q2);
        p1 = _mm256_mul_epu32(p1, q2);
        p2 = _mm256_mul_epu32(p2, q2);
        p3 = _mm256_mul_epu32(p3, q2);

        v0 = _mm256_add_epi64(p0, v0);
        v1 = _mm256_add_epi64(p1, v1);
        v2 = _mm256_add_epi64(p2, v2);
        v3 = _mm256_add_epi64(p3, v3);

        v0 = _mm256_srli_epi64(v0, MONT_MOD_BIT);
        v1 = _mm256_srli_epi64(v1, MONT_MOD_BIT);
        v2 = _mm256_srli_epi64(v2, MONT_MOD_BIT);
        v3 = _mm256_srli_epi64(v3, MONT_MOD_BIT);

        // 3. Cube
        v0 = _mm256_mul_epu32(v0, u0);
        v1 = _mm256_mul_epu32(v1, u1);
        v2 = _mm256_mul_epu32(v2, u2);
        v3 = _mm256_mul_epu32(v3, u3);

        // 4. Montgomery reduction
        p0 = _mm256_and_si256(v0, q0);
        p1 = _mm256_and_si256(v1, q0);
        p2 = _mm256_and_si256(v2, q0);
        p3 = _mm256_and_si256(v3, q0);

        p0 = _mm256_mul_epu32(p0, q1);
        p1 = _mm256_mul_epu32(p1, q1);
        p2 = _mm256_mul_epu32(p2, q1);
        p3 = _mm256_mul_epu32(p3, q1);

        p0 = _mm256_and_si256(p0, q0);
        p1 = _mm256_and_si256(p1, q0);
        p2 = _mm256_and_si256(p2, q0);
        p3 = _mm256_and_si256(p3, q0);

        p0 = _mm256_mul_epu32(p0, q2);
        p1 = _mm256_mul_epu32(p1, q2);
        p2 = _mm256_mul_epu32(p2, q2);
        p3 = _mm256_mul_epu32(p3, q2);

        v0 = _mm256_add_epi64(p0, v0);
        v1 = _mm256_add_epi64(p1, v1);
        v2 = _mm256_add_epi64(p2, v2);
        v3 = _mm256_add_epi64(p3, v3);

        v0 = _mm256_srli_epi64(v0, MONT_MOD_BIT);
        v1 = _mm256_srli_epi64(v1, MONT_MOD_BIT);
        v2 = _mm256_srli_epi64(v2, MONT_MOD_BIT);
        v3 = _mm256_srli_epi64(v3, MONT_MOD_BIT);

        p0 = _mm256_cmpgt_epi32(v0, q3);
        p1 = _mm256_cmpgt_epi32(v1, q3);
        p2 = _mm256_cmpgt_epi32(v2, q3);
        p3 = _mm256_cmpgt_epi32(v3, q3);

        p0 = _mm256_and_si256(p0, q2);
        p1 = _mm256_and_si256(p1, q2);
        p2 = _mm256_and_si256(p2, q2);
        p3 = _mm256_and_si256(p3, q2);

        u0 = _mm256_sub_epi32(v0, p0);
        u1 = _mm256_sub_epi32(v1, p1);
        u2 = _mm256_sub_epi32(v2, p2);
        u3 = _mm256_sub_epi32(v3, p3);
    }

    // Fin: MixColumns -> MixRows -> AddRoundKey

    // MixColumns
    p1 = _mm256_slli_epi32(u1, 1);
    p0 = _mm256_slli_epi32(u0, 1);
    p1 = _mm256_add_epi32(p1, u1);
    p3 = _mm256_add_epi32(u2, u3);
    p2 = _mm256_add_epi32(p1, p0);
    v0 = _mm256_add_epi32(p2, p3);

    p2 = _mm256_slli_epi32(u2, 1);
    p1 = _mm256_slli_epi32(u1, 1);
    p2 = _mm256_add_epi32(p2, u2);
    p3 = _mm256_add_epi32(u0, u3);
    p0 = _mm256_add_epi32(p2, p1);
    v1 = _mm256_add_epi32(p0, p3);

    p3 = _mm256_slli_epi32(u3, 1);
    p2 = _mm256_slli_epi32(u2, 1);
    p3 = _mm256_add_epi32(p3, u3);
    p0 = _mm256_add_epi32(u0, u1);
    p1 = _mm256_add_epi32(p2, p3);
    v2 = _mm256_add_epi32(p0, p1);

    p0 = _mm256_slli_epi32(u0, 1);
    p3 = _mm256_slli_epi32(u3, 1);
    p0 = _mm256_add_epi32(p0, u0);
    p1 = _mm256_add_epi32(u1, u2);
    p2 = _mm256_add_epi32(p0, p3);
    v3 = _mm256_add_epi32(p1, p2);

    // Matrix Transpose
    u0 = _mm256_unpacklo_epi64(v0, v1);
    u1 = _mm256_unpackhi_epi64(v0, v1);
    u2 = _mm256_unpacklo_epi64(v2, v3);
    u3 = _mm256_unpackhi_epi64(v2, v3);

    v0 = _mm256_permute2x128_si256(u0, u2, 0b00100000);
    v1 = _mm256_permute2x128_si256(u1, u3, 0b00100000);
    v2 = _mm256_permute2x128_si256(u0, u2, 0b00110001);
    v3 = _mm256_permute2x128_si256(u1, u3, 0b00110001);

    // MixRows
#if MONT_REDUCTION
    p1 = _mm256_slli_epi64(v1, 1);
    p0 = _mm256_slli_epi64(v0, 1);
    p1 = _mm256_add_epi64(p1, v1);
    p3 = _mm256_add_epi64(v2, v3);

    p2 = _mm256_add_epi64(p0, p3);
    u0 = _mm256_add_epi64(p1, p2);
#else
    p1 = _mm256_slli_epi32(v1, 1);
    p0 = _mm256_slli_epi32(v0, 1);
    p1 = _mm256_add_epi32(p1, v1);
    p3 = _mm256_add_epi32(v2, v3);

    p2 = _mm256_add_epi32(p0, p3);
    u0 = _mm256_add_epi32(p1, p2);
#endif

#if MONT_REDUCTION
    p2 = _mm256_slli_epi64(v2, 1);
    p1 = _mm256_slli_epi64(v1, 1);
    p2 = _mm256_add_epi64(p2, v2);
    p0 = _mm256_add_epi64(v0, v3);

    p3 = _mm256_add_epi64(p1, p0);
    u1 = _mm256_add_epi64(p2, p3);
#else
    p2 = _mm256_slli_epi32(v2, 1);
    p1 = _mm256_slli_epi32(v1, 1);
    p2 = _mm256_add_epi32(p2, v2);
    p0 = _mm256_add_epi32(v0, v3);

    p3 = _mm256_add_epi32(p1, p0);
    u1 = _mm256_add_epi32(p2, p3);
#endif

#if MONT_REDUCTION
    p3 = _mm256_slli_epi64(v3, 1);
    p2 = _mm256_slli_epi64(v2, 1);
    p3 = _mm256_add_epi64(p3, v3);
    p1 = _mm256_add_epi64(v0, v1);

    p0 = _mm256_add_epi64(p2, p1);
    u2 = _mm256_add_epi64(p3, p0);
#else
    p3 = _mm256_slli_epi32(v3, 1);
    p2 = _mm256_slli_epi32(v2, 1);
    p3 = _mm256_add_epi32(p3, v3);
    p1 = _mm256_add_epi32(v0, v1);

    p0 = _mm256_add_epi32(p2, p1);
    u2 = _mm256_add_epi32(p3, p0);
#endif

#if MONT_REDUCTION
    p0 = _mm256_slli_epi64(v0, 1);
    p3 = _mm256_slli_epi64(v3, 1);
    p0 = _mm256_add_epi64(p0, v0);
    p2 = _mm256_add_epi64(v1, v2);

    p1 = _mm256_add_epi64(p3, p2);
    u3 = _mm256_add_epi64(p0, p1);
#else
    p0 = _mm256_slli_epi32(v0, 1);
    p3 = _mm256_slli_epi32(v3, 1);
    p0 = _mm256_add_epi32(p0, v0);
    p2 = _mm256_add_epi32(v1, v2);

    p1 = _mm256_add_epi32(p3, p2);
    u3 = _mm256_add_epi32(p0, p1);
#endif

    // Matrix Transpose
    v0 = _mm256_unpacklo_epi64(u0, u1);
    v1 = _mm256_unpackhi_epi64(u0, u1);
    v2 = _mm256_unpacklo_epi64(u2, u3);
    v3 = _mm256_unpackhi_epi64(u2, u3);

    u0 = _mm256_permute2x128_si256(v0, v2, 0b00100000);
    u1 = _mm256_permute2x128_si256(v1, v3, 0b00100000);
    u2 = _mm256_permute2x128_si256(v0, v2, 0b00110001);
    u3 = _mm256_permute2x128_si256(v1, v3, 0b00110001);

    // AddRoundKey
    v0 = _mm256_load_si256((__m256i *) (round_keys_ + ROUNDS * BLOCKSIZE));
    v1 = _mm256_load_si256((__m256i *) (round_keys_ + ROUNDS * BLOCKSIZE + 4));
    v2 = _mm256_load_si256((__m256i *) (round_keys_ + ROUNDS * BLOCKSIZE + 8));
    v3 = _mm256_load_si256((__m256i *) (round_keys_ + ROUNDS * BLOCKSIZE + 12));

    u0 = _mm256_add_epi64(u0, v0);
    u1 = _mm256_add_epi64(u1, v1);
    u2 = _mm256_add_epi64(u2, v2);
    u3 = _mm256_add_epi64(u3, v3);

    // Montgomery reduction
    v0 = _mm256_and_si256(u0, q0);
    v1 = _mm256_and_si256(u1, q0);
    v2 = _mm256_and_si256(u2, q0);
    v3 = _mm256_and_si256(u3, q0);

    v0 = _mm256_mul_epu32(v0, q1);
    v1 = _mm256_mul_epu32(v1, q1);
    v2 = _mm256_mul_epu32(v2, q1);
    v3 = _mm256_mul_epu32(v3, q1);

    v0 = _mm256_and_si256(v0, q0);
    v1 = _mm256_and_si256(v1, q0);
    v2 = _mm256_and_si256(v2, q0);
    v3 = _mm256_and_si256(v3, q0);

    v0 = _mm256_mul_epu32(v0, q2);
    v1 = _mm256_mul_epu32(v1, q2);
    v2 = _mm256_mul_epu32(v2, q2);
    v3 = _mm256_mul_epu32(v3, q2);

    u0 = _mm256_add_epi64(u0, v0);
    u1 = _mm256_add_epi64(u1, v1);
    u2 = _mm256_add_epi64(u2, v2);
    u3 = _mm256_add_epi64(u3, v3);

    u0 = _mm256_srli_epi64(u0, MONT_MOD_BIT);
    u1 = _mm256_srli_epi64(u1, MONT_MOD_BIT);
    u2 = _mm256_srli_epi64(u2, MONT_MOD_BIT);
    u3 = _mm256_srli_epi64(u3, MONT_MOD_BIT);

    v0 = _mm256_cmpgt_epi32(u0, q3);
    v1 = _mm256_cmpgt_epi32(u1, q3);
    v2 = _mm256_cmpgt_epi32(u2, q3);
    v3 = _mm256_cmpgt_epi32(u3, q3);

    v0 = _mm256_and_si256(v0, q2);
    v1 = _mm256_and_si256(v1, q2);
    v2 = _mm256_and_si256(v2, q2);
    v3 = _mm256_and_si256(v3, q2);

    u0 = _mm256_sub_epi32(u0, v0);
    u1 = _mm256_sub_epi32(u1, v1);
    u2 = _mm256_sub_epi32(u2, v2);
    u3 = _mm256_sub_epi32(u3, v3);

    // Store
    _mm256_store_si256((__m256i *) out, u0);
    _mm256_store_si256((__m256i *) (out + 4), u1);
    _mm256_store_si256((__m256i *) (out + 8), u2);
    _mm256_store_si256((__m256i *) (out + 12), u3);


#else // Not MONTGOMERY
    __m256i u0, u1, u2, u3, v0, v1, v2, v3, p0, p1, p2, p3, q0, q1, q2, q3;

    // Load
    u0 = _mm256_load_si256((__m256i *) INPUT_CONSTANT);
    u1 = _mm256_load_si256((__m256i *) (INPUT_CONSTANT + 4));
    u2 = _mm256_load_si256((__m256i *) (INPUT_CONSTANT + 8));
    u3 = _mm256_load_si256((__m256i *) (INPUT_CONSTANT + 12));


    // Round function
    for (size_t r = 0; r < ROUNDS; r++)
    {
        // AddRoundKey
        v0 = _mm256_load_si256((__m256i *) (round_keys_ + r * BLOCKSIZE));
        v1 = _mm256_load_si256((__m256i *) (round_keys_ + r * BLOCKSIZE + 4));
        v2 = _mm256_load_si256((__m256i *) (round_keys_ + r * BLOCKSIZE + 8));
        v3 = _mm256_load_si256((__m256i *) (round_keys_ + r * BLOCKSIZE + 12));

        u0 = _mm256_add_epi64(u0, v0);
        u1 = _mm256_add_epi64(u1, v1);
        u2 = _mm256_add_epi64(u2, v2);
        u3 = _mm256_add_epi64(u3, v3);


        // MixColumns
        p1 = _mm256_slli_epi64(u1, 1);
        p0 = _mm256_slli_epi64(u0, 1);
        p1 = _mm256_add_epi64(p1, u1);
        p3 = _mm256_add_epi64(u2, u3);
        p2 = _mm256_add_epi64(p1, p0);
        v0 = _mm256_add_epi64(p2, p3);

        p2 = _mm256_slli_epi64(u2, 1);
        p1 = _mm256_slli_epi64(u1, 1);
        p2 = _mm256_add_epi64(p2, u2);
        p3 = _mm256_add_epi64(u0, u3);
        p0 = _mm256_add_epi64(p2, p1);
        v1 = _mm256_add_epi64(p0, p3);

        p3 = _mm256_slli_epi64(u3, 1);
        p2 = _mm256_slli_epi64(u2, 1);
        p3 = _mm256_add_epi64(p3, u3);
        p0 = _mm256_add_epi64(u0, u1);
        p1 = _mm256_add_epi64(p2, p3);
        v2 = _mm256_add_epi64(p0, p1);

        p0 = _mm256_slli_epi64(u0, 1);
        p3 = _mm256_slli_epi64(u3, 1);
        p0 = _mm256_add_epi64(p0, u0);
        p1 = _mm256_add_epi64(u1, u2);
        p2 = _mm256_add_epi64(p0, p3);
        v3 = _mm256_add_epi64(p1, p2);

        // Matrix Transpose
        u0 = _mm256_unpacklo_epi64(v0, v1);
        u1 = _mm256_unpackhi_epi64(v0, v1);
        u2 = _mm256_unpacklo_epi64(v2, v3);
        u3 = _mm256_unpackhi_epi64(v2, v3);

        v0 = _mm256_permute2x128_si256(u0, u2, 0b00100000);
        v1 = _mm256_permute2x128_si256(u1, u3, 0b00100000);
        v2 = _mm256_permute2x128_si256(u0, u2, 0b00110001);
        v3 = _mm256_permute2x128_si256(u1, u3, 0b00110001);

        // MixRows
        p1 = _mm256_slli_epi64(v1, 1);
        p0 = _mm256_slli_epi64(v0, 1);
        p1 = _mm256_add_epi64(p1, v1);
        p3 = _mm256_add_epi64(v2, v3);
        p2 = _mm256_add_epi64(p1, p0);
        u0 = _mm256_add_epi64(p2, p3);

        p2 = _mm256_slli_epi64(v2, 1);
        p1 = _mm256_slli_epi64(v1, 1);
        p2 = _mm256_add_epi64(p2, v2);
        p3 = _mm256_add_epi64(v0, v3);
        p0 = _mm256_add_epi64(p2, p1);
        u1 = _mm256_add_epi64(p0, p3);

        p3 = _mm256_slli_epi64(v3, 1);
        p2 = _mm256_slli_epi64(v2, 1);
        p3 = _mm256_add_epi64(p3, v3);
        p0 = _mm256_add_epi64(v0, v1);
        p1 = _mm256_add_epi64(p2, p3);
        u2 = _mm256_add_epi64(p0, p1);

        p0 = _mm256_slli_epi64(v0, 1);
        p3 = _mm256_slli_epi64(v3, 1);
        p0 = _mm256_add_epi64(p0, v0);
        p1 = _mm256_add_epi64(v1, v2);
        p2 = _mm256_add_epi64(p0, p3);
        u3 = _mm256_add_epi64(p1, p2);

        // Matrix Transpose
        v0 = _mm256_unpacklo_epi64(u0, u1);
        v1 = _mm256_unpackhi_epi64(u0, u1);
        v2 = _mm256_unpacklo_epi64(u2, u3);
        v3 = _mm256_unpackhi_epi64(u2, u3);

        u0 = _mm256_permute2x128_si256(v0, v2, 0b00100000);
        u1 = _mm256_permute2x128_si256(v1, v3, 0b00100000);
        u2 = _mm256_permute2x128_si256(v0, v2, 0b00110001);
        u3 = _mm256_permute2x128_si256(v1, v3, 0b00110001);

/*
Whne using PRIME17, we did not use full reduction which makes element
always less than t. Instead, we use semi-reduction as follows.

When x = a * 2^32 + b * 2^16 + c, then semi-red(x) = a - b + c + t.
 */
#if MODULUS == PRIME17
        p0 = _mm256_set1_epi64x(0x10001);
        p1 = _mm256_set1_epi64x(0xffff);

        // semi-reduction
        // u0 and u1
        v0 = _mm256_srli_epi64(u0, 16);
        v1 = _mm256_srli_epi64(u1, 16);
        v2 = _mm256_srli_epi64(u0, 32);
        v3 = _mm256_srli_epi64(u1, 32);

        q0 = _mm256_and_si256(u0, p1);
        q1 = _mm256_and_si256(u1, p1);
        v0 = _mm256_and_si256(v0, p1);
        v1 = _mm256_and_si256(v1, p1);
        v2 = _mm256_and_si256(v2, p1);
        v3 = _mm256_and_si256(v3, p1);

        u0 = _mm256_add_epi64(p0, q0);
        u1 = _mm256_add_epi64(p0, q1);

        u0 = _mm256_add_epi64(u0, v2);
        u1 = _mm256_add_epi64(u1, v3);

        u0 = _mm256_sub_epi64(u0, v0);
        u1 = _mm256_sub_epi64(u1, v1);

        // u2 and u3
        v0 = _mm256_srli_epi64(u2, 16);
        v1 = _mm256_srli_epi64(u3, 16);
        v2 = _mm256_srli_epi64(u2, 32);
        v3 = _mm256_srli_epi64(u3, 32);

        q0 = _mm256_and_si256(u2, p1);
        q1 = _mm256_and_si256(u3, p1);
        v0 = _mm256_and_si256(v0, p1);
        v1 = _mm256_and_si256(v1, p1);
        v2 = _mm256_and_si256(v2, p1);
        v3 = _mm256_and_si256(v3, p1);

        u2 = _mm256_add_epi64(p0, q0);
        u3 = _mm256_add_epi64(p0, q1);

        u2 = _mm256_add_epi64(u2, v2);
        u3 = _mm256_add_epi64(u3, v3);

        u2 = _mm256_sub_epi64(u2, v0);
        u3 = _mm256_sub_epi64(u3, v1);

#else
        _mm256_store_si256((__m256i *) out, u0);
        _mm256_store_si256((__m256i *) (out + 4), u1);
        _mm256_store_si256((__m256i *) (out + 8), u2);
        _mm256_store_si256((__m256i *) (out + 12), u3);

        for (size_t i = 0; i < BLOCKSIZE; i++)
        {
            out[i] %= MODULUS;
        }

        u0 = _mm256_load_si256((__m256i *) out);
        u1 = _mm256_load_si256((__m256i *) (out + 4));
        u2 = _mm256_load_si256((__m256i *) (out + 8));
        u3 = _mm256_load_si256((__m256i *) (out + 12));

#endif
        // Cube map
        v0 = _mm256_mul_epu32(u0, u0);
        v1 = _mm256_mul_epu32(u1, u1);
        v2 = _mm256_mul_epu32(u2, u2);
        v3 = _mm256_mul_epu32(u3, u3);

#if MODULUS == PRIME17
        // semi-reduction
        // v0 and v1
        q0 = _mm256_srli_epi64(v0, 16);
        q1 = _mm256_srli_epi64(v1, 16);
        q2 = _mm256_srli_epi64(v0, 32);
        q3 = _mm256_srli_epi64(v1, 32);

        p2 = _mm256_and_si256(v0, p1);
        p3 = _mm256_and_si256(v1, p1);
        q0 = _mm256_and_si256(q0, p1);
        q1 = _mm256_and_si256(q1, p1);
        q2 = _mm256_and_si256(q2, p1);
        q3 = _mm256_and_si256(q3, p1);

        v0 = _mm256_add_epi64(p0, p2);
        v1 = _mm256_add_epi64(p0, p3);

        v0 = _mm256_add_epi64(v0, q2);
        v1 = _mm256_add_epi64(v1, q3);

        v0 = _mm256_sub_epi64(v0, q0);
        v1 = _mm256_sub_epi64(v1, q1);

        // v2 and v3
        q0 = _mm256_srli_epi64(v2, 16);
        q1 = _mm256_srli_epi64(v3, 16);
        q2 = _mm256_srli_epi64(v2, 32);
        q3 = _mm256_srli_epi64(v3, 32);

        p2 = _mm256_and_si256(v2, p1);
        p3 = _mm256_and_si256(v3, p1);
        q0 = _mm256_and_si256(q0, p1);
        q1 = _mm256_and_si256(q1, p1);
        q2 = _mm256_and_si256(q2, p1);
        q3 = _mm256_and_si256(q3, p1);

        v2 = _mm256_add_epi64(p0, p2);
        v3 = _mm256_add_epi64(p0, p3);

        v2 = _mm256_add_epi64(v2, q2);
        v3 = _mm256_add_epi64(v3, q3);

        v2 = _mm256_sub_epi64(v2, q0);
        v3 = _mm256_sub_epi64(v3, q1);

#else
        _mm256_store_si256((__m256i *) out, v0);
        _mm256_store_si256((__m256i *) (out + 4), v1);
        _mm256_store_si256((__m256i *) (out + 8), v2);
        _mm256_store_si256((__m256i *) (out + 12), v3);

        for (size_t i = 0; i < BLOCKSIZE; i++)
        {
            out[i] %= MODULUS;
        }

        v0 = _mm256_load_si256((__m256i *) out);
        v1 = _mm256_load_si256((__m256i *) (out + 4));
        v2 = _mm256_load_si256((__m256i *) (out + 8));
        v3 = _mm256_load_si256((__m256i *) (out + 12));

#endif
        u0 = _mm256_mul_epu32(v0, u0);
        u1 = _mm256_mul_epu32(v1, u1);
        u2 = _mm256_mul_epu32(v2, u2);
        u3 = _mm256_mul_epu32(v3, u3);

#if MODULUS == PRIME17
        // semi-reduction
        // u0 and u1
        q0 = _mm256_srli_epi64(u0, 16);
        q1 = _mm256_srli_epi64(u1, 16);
        q2 = _mm256_srli_epi64(u0, 32);
        q3 = _mm256_srli_epi64(u1, 32);

        p2 = _mm256_and_si256(u0, p1);
        p3 = _mm256_and_si256(u1, p1);
        q0 = _mm256_and_si256(q0, p1);
        q1 = _mm256_and_si256(q1, p1);
        q2 = _mm256_and_si256(q2, p1);
        q3 = _mm256_and_si256(q3, p1);

        u0 = _mm256_add_epi64(p0, p2);
        u1 = _mm256_add_epi64(p0, p3);

        u0 = _mm256_add_epi64(u0, q2);
        u1 = _mm256_add_epi64(u1, q3);

        u0 = _mm256_sub_epi64(u0, q0);
        u1 = _mm256_sub_epi64(u1, q1);

        // u2 and u3
        q0 = _mm256_srli_epi64(u2, 16);
        q1 = _mm256_srli_epi64(u3, 16);
        q2 = _mm256_srli_epi64(u2, 32);
        q3 = _mm256_srli_epi64(u3, 32);

        p2 = _mm256_and_si256(u2, p1);
        p3 = _mm256_and_si256(u3, p1);
        q0 = _mm256_and_si256(q0, p1);
        q1 = _mm256_and_si256(q1, p1);
        q2 = _mm256_and_si256(q2, p1);
        q3 = _mm256_and_si256(q3, p1);

        u2 = _mm256_add_epi64(p0, p2);
        u3 = _mm256_add_epi64(p0, p3);

        u2 = _mm256_add_epi64(u2, q2);
        u3 = _mm256_add_epi64(u3, q3);

        u2 = _mm256_sub_epi64(u2, q0);
        u3 = _mm256_sub_epi64(u3, q1);

#else
        _mm256_store_si256((__m256i *) out, u0);
        _mm256_store_si256((__m256i *) (out + 4), u1);
        _mm256_store_si256((__m256i *) (out + 8), u2);
        _mm256_store_si256((__m256i *) (out + 12), u3);

        for (size_t i = 0; i < BLOCKSIZE; i++)
        {
            out[i] %= MODULUS;
        }

        u0 = _mm256_load_si256((__m256i *) out);
        u1 = _mm256_load_si256((__m256i *) (out + 4));
        u2 = _mm256_load_si256((__m256i *) (out + 8));
        u3 = _mm256_load_si256((__m256i *) (out + 12));
#endif
    }

    // Fin: MixColumns -> MixRows -> AddRoundKey

    // MixColumns
    p1 = _mm256_slli_epi64(u1, 1);
    p0 = _mm256_slli_epi64(u0, 1);
    p1 = _mm256_add_epi64(p1, u1);
    p3 = _mm256_add_epi64(u2, u3);
    p2 = _mm256_add_epi64(p1, p0);
    v0 = _mm256_add_epi64(p2, p3);

    p2 = _mm256_slli_epi64(u2, 1);
    p1 = _mm256_slli_epi64(u1, 1);
    p2 = _mm256_add_epi64(p2, u2);
    p3 = _mm256_add_epi64(u0, u3);
    p0 = _mm256_add_epi64(p2, p1);
    v1 = _mm256_add_epi64(p0, p3);

    p3 = _mm256_slli_epi64(u3, 1);
    p2 = _mm256_slli_epi64(u2, 1);
    p3 = _mm256_add_epi64(p3, u3);
    p0 = _mm256_add_epi64(u0, u1);
    p1 = _mm256_add_epi64(p2, p3);
    v2 = _mm256_add_epi64(p0, p1);

    p0 = _mm256_slli_epi64(u0, 1);
    p3 = _mm256_slli_epi64(u3, 1);
    p0 = _mm256_add_epi64(p0, u0);
    p1 = _mm256_add_epi64(u1, u2);
    p2 = _mm256_add_epi64(p0, p3);
    v3 = _mm256_add_epi64(p1, p2);

    // Matrix Transpose
    u0 = _mm256_unpacklo_epi64(v0, v1);
    u1 = _mm256_unpackhi_epi64(v0, v1);
    u2 = _mm256_unpacklo_epi64(v2, v3);
    u3 = _mm256_unpackhi_epi64(v2, v3);

    v0 = _mm256_permute2x128_si256(u0, u2, 0b00100000);
    v1 = _mm256_permute2x128_si256(u1, u3, 0b00100000);
    v2 = _mm256_permute2x128_si256(u0, u2, 0b00110001);
    v3 = _mm256_permute2x128_si256(u1, u3, 0b00110001);

    // MixRows
    p1 = _mm256_slli_epi64(v1, 1);
    p0 = _mm256_slli_epi64(v0, 1);
    p1 = _mm256_add_epi64(p1, v1);
    p3 = _mm256_add_epi64(v2, v3);
    p2 = _mm256_add_epi64(p1, p0);
    u0 = _mm256_add_epi64(p2, p3);

    p2 = _mm256_slli_epi64(v2, 1);
    p1 = _mm256_slli_epi64(v1, 1);
    p2 = _mm256_add_epi64(p2, v2);
    p3 = _mm256_add_epi64(v0, v3);
    p0 = _mm256_add_epi64(p2, p1);
    u1 = _mm256_add_epi64(p0, p3);

    p3 = _mm256_slli_epi64(v3, 1);
    p2 = _mm256_slli_epi64(v2, 1);
    p3 = _mm256_add_epi64(p3, v3);
    p0 = _mm256_add_epi64(v0, v1);
    p1 = _mm256_add_epi64(p2, p3);
    u2 = _mm256_add_epi64(p0, p1);

    p0 = _mm256_slli_epi64(v0, 1);
    p3 = _mm256_slli_epi64(v3, 1);
    p0 = _mm256_add_epi64(p0, v0);
    p1 = _mm256_add_epi64(v1, v2);
    p2 = _mm256_add_epi64(p0, p3);
    u3 = _mm256_add_epi64(p1, p2);

    // Matrix Transpose
    v0 = _mm256_unpacklo_epi64(u0, u1);
    v1 = _mm256_unpackhi_epi64(u0, u1);
    v2 = _mm256_unpacklo_epi64(u2, u3);
    v3 = _mm256_unpackhi_epi64(u2, u3);

    u0 = _mm256_permute2x128_si256(v0, v2, 0b00100000);
    u1 = _mm256_permute2x128_si256(v1, v3, 0b00100000);
    u2 = _mm256_permute2x128_si256(v0, v2, 0b00110001);
    u3 = _mm256_permute2x128_si256(v1, v3, 0b00110001);

    // AddRoundKey
    v0 = _mm256_load_si256((__m256i *) (round_keys_ + ROUNDS * BLOCKSIZE));
    v1 = _mm256_load_si256((__m256i *) (round_keys_ + ROUNDS * BLOCKSIZE + 4));
    v2 = _mm256_load_si256((__m256i *) (round_keys_ + ROUNDS * BLOCKSIZE + 8));
    v3 = _mm256_load_si256((__m256i *) (round_keys_ + ROUNDS * BLOCKSIZE + 12));

    u0 = _mm256_add_epi64(u0, v0);
    u1 = _mm256_add_epi64(u1, v1);
    u2 = _mm256_add_epi64(u2, v2);
    u3 = _mm256_add_epi64(u3, v3);

#if MODULUS == PRIME17
    p0 = _mm256_set1_epi64x(0x10001);
    p1 = _mm256_set1_epi64x(0xffff);

    // full-reduction
    // u0 and u1
    v0 = _mm256_srli_epi64(u0, 16);
    v1 = _mm256_srli_epi64(u1, 16);
    v2 = _mm256_srli_epi64(u0, 32);
    v3 = _mm256_srli_epi64(u1, 32);

    q0 = _mm256_and_si256(u0, p1);
    q1 = _mm256_and_si256(u1, p1);
    v0 = _mm256_and_si256(v0, p1);
    v1 = _mm256_and_si256(v1, p1);
    v2 = _mm256_and_si256(v2, p1);
    v3 = _mm256_and_si256(v3, p1);

    u0 = _mm256_add_epi64(q0, v2);
    u1 = _mm256_add_epi64(q1, v3);

    p2 = _mm256_cmpgt_epi32(v0, u0);
    p3 = _mm256_cmpgt_epi32(v1, u1);

    q0 = _mm256_and_si256(p2, p0);
    q1 = _mm256_and_si256(p3, p0);

    u0 = _mm256_add_epi64(u0, q0);
    u1 = _mm256_add_epi64(u1, q1);

    u0 = _mm256_sub_epi64(u0, v0);
    u1 = _mm256_sub_epi64(u1, v1);

    // u2 and u3
    v0 = _mm256_srli_epi64(u2, 16);
    v1 = _mm256_srli_epi64(u3, 16);
    v2 = _mm256_srli_epi64(u2, 32);
    v3 = _mm256_srli_epi64(u3, 32);

    q0 = _mm256_and_si256(u2, p1);
    q1 = _mm256_and_si256(u3, p1);
    v0 = _mm256_and_si256(v0, p1);
    v1 = _mm256_and_si256(v1, p1);
    v2 = _mm256_and_si256(v2, p1);
    v3 = _mm256_and_si256(v3, p1);

    u2 = _mm256_add_epi64(q0, v2);
    u3 = _mm256_add_epi64(q1, v3);

    p2 = _mm256_cmpgt_epi32(v0, u2);
    p3 = _mm256_cmpgt_epi32(v1, u3);

    q0 = _mm256_and_si256(p2, p0);
    q1 = _mm256_and_si256(p3, p0);

    u2 = _mm256_add_epi64(u2, q0);
    u3 = _mm256_add_epi64(u3, q1);

    u2 = _mm256_sub_epi64(u2, v0);
    u3 = _mm256_sub_epi64(u3, v1);

    // Store
    _mm256_store_si256((__m256i *) out, u0);
    _mm256_store_si256((__m256i *) (out + 4), u1);
    _mm256_store_si256((__m256i *) (out + 8), u2);
    _mm256_store_si256((__m256i *) (out + 12), u3);
#else
    _mm256_store_si256((__m256i *) out, u0);
    _mm256_store_si256((__m256i *) (out + 4), u1);
    _mm256_store_si256((__m256i *) (out + 8), u2);
    _mm256_store_si256((__m256i *) (out + 12), u3);

    for (size_t i = 0; i < BLOCKSIZE; i++)
    {
        out[i] = out[i] % MODULUS;
    }
#endif
#endif

}

void Hera::get_rand_vectors(uint64_t *output)
{
    memcpy(output, rand_vectors_, sizeof(uint64_t) * XOF_ELEMENT_COUNT);
}

// Hera private functions
void Hera::keyschedule()
{
    /*
    Since PRIME17 == 2 ** 16 + 1, random element in Z_t is chosen by
    sampling 16-bit random string and add 1. Using this method, we
    can sample uniformly random nonzero sample in Z_t.
     */
#if MODULUS == PRIME17
    // Set random vector from SHAKE256
    __attribute__((aligned(32))) uint8_t buf[4 * RATE_IN_BYTE * NUM_SQUEEZE];
    memset(buf, 0, 4 * RATE_IN_BYTE * NUM_SQUEEZE);
    shake_->squeeze(buf, 4 * RATE_IN_BYTE * NUM_SQUEEZE);

    int offset = 0;
    uint64_t elem;
    for (int i = 0; i < XOF_ELEMENT_COUNT; i++)
    {
        elem = *(uint16_t*)(buf + offset) + 1;
        offset += 2;
        rand_vectors_[i] = elem;
    }

    // Compute round keys
    // round_keys[i] = (key_[i] * rand_vectors_[i]) % MODULUS

    __m256i u0, u1, u2, u3, v0, v1, v2, v3, p0, p1, w0, w1, w2, w3, w4, w5;

    p0 = _mm256_set1_epi64x(0x10001);
    p1 = _mm256_set1_epi64x(0xffff);

    u0 = _mm256_load_si256((__m256i *) key_);
    u1 = _mm256_load_si256((__m256i *) (key_ + 4));
    u2 = _mm256_load_si256((__m256i *) (key_ + 8));
    u3 = _mm256_load_si256((__m256i *) (key_ + 12));

    for (int i = 0; i < XOF_ELEMENT_COUNT; i += 16)
    {
        // Index 0 to 7
        v0 = _mm256_load_si256((__m256i *) (rand_vectors_ + i + 0x0));
        v1 = _mm256_load_si256((__m256i *) (rand_vectors_ + i + 0x4));

        v0 = _mm256_mul_epu32(u0, v0);
        v1 = _mm256_mul_epu32(u1, v1);

        w2 = _mm256_srli_epi64(v0, 16);
        w3 = _mm256_srli_epi64(v1, 16);
        w4 = _mm256_srli_epi64(v0, 32);
        w5 = _mm256_srli_epi64(v1, 32);

        w0 = _mm256_and_si256(v0, p1);
        w1 = _mm256_and_si256(v1, p1);
        w2 = _mm256_and_si256(w2, p1);
        w3 = _mm256_and_si256(w3, p1);
        w4 = _mm256_and_si256(w4, p1);
        w5 = _mm256_and_si256(w5, p1);

        v0 = _mm256_add_epi64(w0, w4);
        v1 = _mm256_add_epi64(w1, w5);

        v2 = _mm256_cmpgt_epi32(w2, v0);
        v3 = _mm256_cmpgt_epi32(w3, v1);

        w0 = _mm256_and_si256(v2, p0);
        w1 = _mm256_and_si256(v3, p0);

        v0 = _mm256_add_epi64(v0, w0);
        v1 = _mm256_add_epi64(v1, w1);

        v2 = _mm256_sub_epi64(v0, w2);
        v3 = _mm256_sub_epi64(v1, w3);

        _mm256_store_si256((__m256i *) (round_keys_ + i + 0x0), v2);
        _mm256_store_si256((__m256i *) (round_keys_ + i + 0x4), v3);

        // Index 8 to 15
        v0 = _mm256_load_si256((__m256i *) (rand_vectors_ + i + 0x8));
        v1 = _mm256_load_si256((__m256i *) (rand_vectors_ + i + 0xc));

        v0 = _mm256_mul_epu32(u2, v0);
        v1 = _mm256_mul_epu32(u3, v1);

        w2 = _mm256_srli_epi64(v0, 16);
        w3 = _mm256_srli_epi64(v1, 16);
        w4 = _mm256_srli_epi64(v0, 32);
        w5 = _mm256_srli_epi64(v1, 32);

        w0 = _mm256_and_si256(v0, p1);
        w1 = _mm256_and_si256(v1, p1);
        w2 = _mm256_and_si256(w2, p1);
        w3 = _mm256_and_si256(w3, p1);
        w4 = _mm256_and_si256(w4, p1);
        w5 = _mm256_and_si256(w5, p1);

        v0 = _mm256_add_epi64(w0, w4);
        v1 = _mm256_add_epi64(w1, w5);

        v2 = _mm256_cmpgt_epi32(w2, v0);
        v3 = _mm256_cmpgt_epi32(w3, v1);

        w0 = _mm256_and_si256(v2, p0);
        w1 = _mm256_and_si256(v3, p0);

        v0 = _mm256_add_epi64(v0, w0);
        v1 = _mm256_add_epi64(v1, w1);

        v2 = _mm256_sub_epi64(v0, w2);
        v3 = _mm256_sub_epi64(v1, w3);

        _mm256_store_si256((__m256i *) (round_keys_ + i + 0x8), v2);
        _mm256_store_si256((__m256i *) (round_keys_ + i + 0xc), v3);
    }

#else
    // Set random vector from SHAKE256
    const __m256i zero = _mm256_setzero_si256();
    const __m256i modulus32 = _mm256_set1_epi32(MODULUS);
    const __m256i all_one_mask32 = _mm256_set1_epi32((1ULL << MOD_BIT_COUNT) - 1);
    __m256i u0, u1, u2, u3, v0, v1, v2, v3, idx;

    __attribute__((aligned(32))) uint8_t buf[4 * RATE_IN_BYTE];

    unsigned int offset = 0;
    unsigned int squeeze_byte = 0;
    __attribute__((aligned(32))) uint32_t idx_arr[8] = {0};

    shake_->squeeze(buf, 4 * RATE_IN_BYTE);

    int ctr = (4 * RATE_IN_BYTE) / 32;
    while (offset < XOF_ELEMENT_COUNT)
    {
        u0 = _mm256_loadu_si256((__m256i *)(buf + squeeze_byte));
#if MODULUS == PRIME32
        u1 = _mm256_cmpgt_epi32(u0, zero);
        u2 = _mm256_cmpgt_epi32(modulus32, u0);
        u1 = _mm256_or_si256(u1, u2);
#else
        u0 = _mm256_and_si256(u0, all_one_mask32);
        u1 = _mm256_cmpgt_epi32(modulus32, u0);
#endif
        int good_idx = _mm256_movemask_ps((__m256) u1);

        int good_cnt = 0;
        int bad_cnt = 0;
        for (int i = 0; i < 8; i++)
        {
            if (good_idx & (1 << i))
            {
                idx_arr[unpack_order[good_cnt]] = i;
                good_cnt++;
            }
            else
            {
                idx_arr[unpack_order[7 - bad_cnt]] = i;
                bad_cnt++;
            }
        }
        idx = _mm256_load_si256((__m256i *) idx_arr);
        u2 = _mm256_permutevar8x32_epi32(u0, idx);

        u0 = _mm256_unpackhi_epi32(u2, zero);
        u1 = _mm256_unpacklo_epi32(u2, zero);
        _mm256_storeu_si256((__m256i*) (rand_vectors_ + offset), u0);
        _mm256_storeu_si256((__m256i*) (rand_vectors_ + offset + 4), u1);

        offset += good_cnt;
        squeeze_byte += 32;
        ctr--;
        if (ctr == 0)
        {
            shake_->squeeze(buf, 4 * RATE_IN_BYTE);
            squeeze_byte = 0;
            ctr = (4 * RATE_IN_BYTE) / 32;
        }
    }

    // Compute round keys
    // round_keys[i] = (key_[i] * rand_vectors_[i]) % MODULUS;
    u0 = _mm256_load_si256((__m256i *) key_);
    u1 = _mm256_load_si256((__m256i *) (key_ + 4));
    u2 = _mm256_load_si256((__m256i *) (key_ + 8));
    u3 = _mm256_load_si256((__m256i *) (key_ + 12));

    for (int r = 0; r <= ROUNDS; r++)
    {
        v0 = _mm256_load_si256((__m256i *) (rand_vectors_ + BLOCKSIZE * r));
        v1 = _mm256_load_si256((__m256i *) (rand_vectors_ + BLOCKSIZE * r + 4));
        v2 = _mm256_load_si256((__m256i *) (rand_vectors_ + BLOCKSIZE * r + 8));
        v3 = _mm256_load_si256((__m256i *) (rand_vectors_ + BLOCKSIZE * r + 12));

        v0 = _mm256_mul_epu32(u0, v0);
        v1 = _mm256_mul_epu32(u1, v1);
        v2 = _mm256_mul_epu32(u2, v2);
        v3 = _mm256_mul_epu32(u3, v3);

        _mm256_store_si256((__m256i *) (round_keys_ + BLOCKSIZE * r), v0);
        _mm256_store_si256((__m256i *) (round_keys_ + BLOCKSIZE * r + 4), v1);
        _mm256_store_si256((__m256i *) (round_keys_ + BLOCKSIZE * r + 8), v2);
        _mm256_store_si256((__m256i *) (round_keys_ + BLOCKSIZE * r + 12), v3);
    }

    for (int i = 0; i < (ROUNDS + 1) * BLOCKSIZE; i++)
    {
        round_keys_[i] %= MODULUS;
    }
#endif
}

