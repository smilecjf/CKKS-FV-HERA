/*
-----------------------------------------------------------------------
This source code is excerpted and modified from Microsoft SEAL library
version 3.4.5.

Copyright (c) Microsoft Corporation. All rights reserved.
Licensed under the MIT license.
-----------------------------------------------------------------------
 */

#include "FV_Encoder.h"
#include "seal/util/uintarithsmallmod.h"

using namespace std;

/* FV_Encoder Public Functions */
void FV_Encoder::init()
{
    // Compute coeff_count_power
    coeff_count_power_ = seal::util::get_power_of_two(poly_modulus_degree_fv);

    if (!seal::util::try_minimal_primitive_root(2 * poly_modulus_degree_fv, plain_modulus, root_))
    {
        // error
    }
    if (!seal::util::try_mod_inverse(root_, plain_modulus, inverse_root_))
    {
        // error
    }

    // Populate the tables storing (scaled version of) powers of root
    // mod q in bit-scrambled order.
    ntt_powers_of_primitive_root(root_, root_powers_);
    ntt_scale_powers_of_primitive_root(root_powers_,
        scaled_root_powers_);

    // Populate the tables storing (scaled version of) powers of
    // (root)^{-1} mod q in bit-scrambled order.
    ntt_powers_of_primitive_root(inverse_root_, inv_root_powers_);
    ntt_scale_powers_of_primitive_root(inv_root_powers_,
        scaled_inv_root_powers_);

    // Populate the tables storing (scaled version of ) 2 times
    // powers of roots^-1 mod q  in bit-scrambled order.
    for (size_t i = 0; i < poly_modulus_degree_fv; i++)
    {
        inv_root_powers_div_two_[i] =
            seal::util::div2_uint_mod(inv_root_powers_[i], plain_modulus);
    }
    ntt_scale_powers_of_primitive_root(inv_root_powers_div_two_,
        scaled_inv_root_powers_div_two_);

    populate_matrix_reps_index_map();
    populate_roots_of_unity_vector();
}

vector<uint64_t> FV_Encoder::encode(vector<uint64_t> data)
{
    uint64_t data_arr[poly_modulus_degree_fv];
    uint64_t destination[poly_modulus_degree_fv];
    size_t num_data = data.size();

    for (size_t i = 0; i < num_data; i++)
    {
        data_arr[i] = data[i];
    }
    for (size_t i = num_data; i < poly_modulus_degree_fv; i++)
    {
        data_arr[i] = size_t(0);
    }

    for (size_t i = 0; i < poly_modulus_degree_fv; i++)
    {
        *(destination + matrix_reps_index_map_[i]) = data_arr[i];
    }

    inverse_ntt_negacyclic_harvey(destination);

    vector<uint64_t> encoded;
    for (size_t i = 0; i < poly_modulus_degree_fv; i++)
    {
        encoded.push_back(destination[i]);
    }

    return encoded;
}

void FV_Encoder::encode(uint64_t *data, uint64_t *destination)
{
    for (size_t i = 0; i < poly_modulus_degree_fv; i++)
    {
        *(destination + matrix_reps_index_map_[i]) = data[i];
    }

    inverse_ntt_negacyclic_harvey(destination);
}

/* FV_Encoder Private Functions */
// Initialization functions
void FV_Encoder::populate_matrix_reps_index_map()
{
    int logn = seal::util::get_power_of_two(poly_modulus_degree_fv);

    size_t row_size = poly_modulus_degree_fv >> 1;
    size_t m = poly_modulus_degree_fv << 1;
    uint64_t gen = 3;
    uint64_t pos = 1;
    for (size_t i = 0; i < row_size; i++)
    {
        // Position in normal bit order
        uint64_t index1 = (pos - 1) >> 1;
        uint64_t index2 = (m - pos - 1) >> 1;

        // Set the bit-reversed locations
        matrix_reps_index_map_[i] = seal::util::safe_cast<uint64_t>(seal::util::reverse_bits(index1, logn));
        matrix_reps_index_map_[row_size | i] = seal::util::safe_cast<size_t>(seal::util::reverse_bits(index2, logn));

        // Next primitive root
        pos *= gen;
        pos &= (m - 1);
    }
}

void FV_Encoder::populate_roots_of_unity_vector()
{
    uint64_t generator_sq = seal::util::multiply_uint_uint_mod(root_, root_, plain_modulus);
    roots_of_unity_[0] = root_;

    for (size_t i = 1; i < poly_modulus_degree_fv; i++)
    {
        roots_of_unity_[i] = seal::util::multiply_uint_uint_mod(roots_of_unity_[i - 1],
                                                                generator_sq, plain_modulus);
    }
}

// Utility functions
void FV_Encoder::inverse_ntt_negacyclic_harvey_lazy(uint64_t *operand)
{
    uint64_t two_times_modulus = plain_modulus * 2;

    // return the bit-reversed order of NTT.
    size_t n = poly_modulus_degree_fv;
    size_t t = 1;

    for (size_t m = n; m > 1; m >>= 1)
    {
        size_t j1 = 0;
        size_t h = m >> 1;
        if (t >= 4)
        {
            for (size_t i = 0; i < h; i++)
            {
                size_t j2 = j1 + t;
                const uint64_t W = inv_root_powers_div_two_[h + i];
                const uint64_t Wprime = scaled_inv_root_powers_div_two_[h + i];


                uint64_t *U = operand + j1;
                uint64_t *V = U + t;
                uint64_t currU;
                uint64_t T;
                unsigned long long H;
                for (size_t j = j1; j < j2; j += 4)
                {
                    T = two_times_modulus - *V + *U;
                    currU = *U + *V - (two_times_modulus & static_cast<uint64_t>(-static_cast<int64_t>((*U << 1) >= T)));
                    *U++ = (currU + (plain_modulus & static_cast<uint64_t>(-static_cast<int64_t>(T & 1)))) >> 1;
                    seal::util::multiply_uint64_hw64(Wprime, T, &H);
                    *V++ = T * W - H * plain_modulus;

                    T = two_times_modulus - *V + *U;
                    currU = *U + *V - (two_times_modulus & static_cast<uint64_t>(-static_cast<int64_t>((*U << 1) >= T)));
                    *U++ = (currU + (plain_modulus & static_cast<uint64_t>(-static_cast<int64_t>(T & 1)))) >> 1;
                    seal::util::multiply_uint64_hw64(Wprime, T, &H);
                    *V++ = T * W - H * plain_modulus;

                    T = two_times_modulus - *V + *U;
                    currU = *U + *V - (two_times_modulus & static_cast<uint64_t>(-static_cast<int64_t>((*U << 1) >= T)));
                    *U++ = (currU + (plain_modulus & static_cast<uint64_t>(-static_cast<int64_t>(T & 1)))) >> 1;
                    seal::util::multiply_uint64_hw64(Wprime, T, &H);
                    *V++ = T * W - H * plain_modulus;

                    T = two_times_modulus - *V + *U;
                    currU = *U + *V - (two_times_modulus & static_cast<uint64_t>(-static_cast<int64_t>((*U << 1) >= T)));
                    *U++ = (currU + (plain_modulus & static_cast<uint64_t>(-static_cast<int64_t>(T & 1)))) >> 1;
                    seal::util::multiply_uint64_hw64(Wprime, T, &H);
                    *V++ = T * W - H * plain_modulus;
                }
                j1 += (t << 1);
            }
        }
        else
        {
            for (size_t i = 0; i < h; i++)
            {
                size_t j2 = j1 + t;
                const uint64_t W = inv_root_powers_div_two_[h + i];
                const uint64_t Wprime = scaled_inv_root_powers_div_two_[h + i];

                uint64_t *U = operand + j1;
                uint64_t *V = U + t;
                uint64_t currU;
                uint64_t T;
                unsigned long long H;
                for (size_t j = j1; j < j2; j++)
                {
                    T = two_times_modulus - *V + *U;
                    currU = *U + *V - (two_times_modulus & static_cast<uint64_t>(-static_cast<int64_t>((*U << 1) >= T)));
                    *U++ = (currU + (plain_modulus & static_cast<uint64_t>(-static_cast<int64_t>(T & 1)))) >> 1;

                    seal::util::multiply_uint64_hw64(Wprime, T, &H);
                    *V++ = W * T - H * plain_modulus;
                }
                j1 += (t << 1);
            }
        }
        t <<= 1;
    }

}

inline void FV_Encoder::inverse_ntt_negacyclic_harvey(uint64_t *operand)
{
    inverse_ntt_negacyclic_harvey_lazy(operand);

    size_t n = poly_modulus_degree_fv;

    for (; n--; operand++)
    {
        if (*operand >= plain_modulus)
        {
            *operand -= plain_modulus;
        }
    }

}

void FV_Encoder::ntt_powers_of_primitive_root(uint64_t root, uint64_t *destination)
{
    uint64_t *destination_start = destination;
    *destination_start = 1;
    for (size_t i = 1; i < poly_modulus_degree_fv; i++)
    {
        uint64_t *next_destination = destination_start + seal::util::reverse_bits(i, coeff_count_power_);
        *next_destination = seal::util::multiply_uint_uint_mod(*destination, root, plain_modulus);
        destination = next_destination;
    }
}

void FV_Encoder::ntt_scale_powers_of_primitive_root(const uint64_t *input, uint64_t *destination)
{
    for (size_t i = 0; i < poly_modulus_degree_fv; i++, input++, destination++)
    {
        uint64_t wide_quotient[2]{ 0, 0 };
        uint64_t wide_coeff[2]{ 0, *input };
        seal::util::divide_uint128_uint64_inplace(wide_coeff, plain_modulus, wide_quotient);
        *destination = wide_quotient[0];
    }
}
