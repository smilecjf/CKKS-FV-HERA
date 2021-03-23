/*
-----------------------------------------------------------------------
This source code is excerpted and modified from Microsoft SEAL library
version 3.4.5.

Copyright (c) Microsoft Corporation. All rights reserved.
Licensed under the MIT license.
-----------------------------------------------------------------------
 */


#include "CKKS_Encoder.h"

/* CKKS_Encoder Public Functions */
void CKKS_Encoder::init()
{
    slot_count_ = poly_modulus_degree / 2;
    slot_count_power_of_two_ = get_power_of_two(slot_count_);
    coeff_count_power_of_two_ = slot_count_power_of_two_ + 1;
    populate_reps_index_map();
}

void CKKS_Encoder::encode(const complex<double> *values, double scale,
    uint64_t *destination, size_t plain_modulus)
{
    complex<double> conj_values[poly_modulus_degree];
    for (size_t i = 0; i < slot_count_; i++)
    {
        conj_values[matrix_reps_index_map_[i]] = values[i];
        conj_values[matrix_reps_index_map_[i + slot_count_]] = conj(values[i]);
    }

    size_t tt = 1;
    for (int i = 0; i < coeff_count_power_of_two_; i++)
    {
        size_t m = size_t(1) << (coeff_count_power_of_two_ - i);
        size_t k_start = 0;
        size_t h = m / 2;

        for (size_t j = 0; j < h; j++)
        {
            size_t k_end = k_start + tt;
            auto s = inv_roots_[h + j];

            for (size_t k = k_start; k < k_end; k++)
            {
                auto u = conj_values[k];
                auto v = conj_values[k + tt];
                conj_values[k] = u + v;
                conj_values[k + tt] = (u - v) * s;
            }

            k_start += 2 * tt;
        }
        tt *= 2;
    }

    double n_inv = double(1.0) / static_cast<double>(poly_modulus_degree);

    // Put the scale in at this point
    n_inv *= scale;

    for (size_t i = 0; i < poly_modulus_degree; i++)
    {
        // multiply by scale and n_inv (see above)
        conj_values[i] *= n_inv;

        double coeffd = round(conj_values[i].real());
        bool is_negative = signbit(coeffd);

        uint64_t coeffu = static_cast<uint64_t>(fabs(coeffd));
        if (is_negative)
        {
            destination[i] = plain_modulus - (coeffu % plain_modulus);
        }
        else
        {
            destination[i] = coeffu;
        }
    }

}

/* CKKS_Encoder Private Functions */
// Initialization functions
void CKKS_Encoder::populate_reps_index_map()
{
    // Copy from the matrix to the value vectors
    uint64_t gen = 3;
    uint64_t pos = 1;
    uint64_t m = static_cast<uint64_t>(poly_modulus_degree) << 1;
    for (size_t i = 0; i < slot_count_; i++)
    {
        // Position in normal bit order
        uint64_t index1 = (pos - 1) >> 1;
        uint64_t index2 = ( - pos - 1) >> 1;

        // Set the bit-reversed locations
        matrix_reps_index_map_[i] = reverse_bits(index1, coeff_count_power_of_two_);
        matrix_reps_index_map_[slot_count_ | i] = reverse_bits(index2,
            coeff_count_power_of_two_);

        // Next primitive root
        pos *= gen;
        pos &= (m - 1);
    }

    for (size_t i = 0; i < poly_modulus_degree; i++)
    {
        roots_[i] = get_root(reverse_bits(i, coeff_count_power_of_two_), static_cast<size_t>(m));
        inv_roots_[i] = conj(roots_[i]);
    }
}

complex<double> CKKS_Encoder::get_root(size_t index)
{
    if (index <= root_fidelity_ / 8)
    {
        return std::polar<double>(
            1.0,
            2 * PI_ * static_cast<double>(index) / root_fidelity_);
    }
    else if (index <= root_fidelity_ / 4)
    {
        auto mirror_root = get_root(root_fidelity_ / 4 - index);
        return {mirror_root.imag(), mirror_root.real()};
    }
    else if (index <= root_fidelity_ / 2)
    {
        return -conj(get_root(root_fidelity_ / 2 - index));
    }
    else if (index <= 3 * (root_fidelity_ / 2))
    {
        return -get_root(index - root_fidelity_ / 2);
    }
    else
    {
        return conj(get_root(root_fidelity_ - index));
    }
}

inline complex<double> CKKS_Encoder::get_root(size_t index, size_t of_roots)
{
    index &= of_roots - 1;
    return get_root(index * root_fidelity_ / of_roots);
}


// Utility functions
int CKKS_Encoder::get_power_of_two(size_t value)
{
    return floor(log2(value));
}

inline uint32_t CKKS_Encoder::reverse_bits(uint32_t operand)
{
    operand = (((operand & (uint32_t)(0xaaaaaaaa)) >> 1) | ((operand & (uint32_t)(0x55555555)) << 1));
    operand = (((operand & (uint32_t)(0xcccccccc)) >> 2) | ((operand & (uint32_t)(0x33333333)) << 2));
    operand = (((operand & (uint32_t)(0xf0f0f0f0)) >> 4) | ((operand & (uint32_t)(0x0f0f0f0f)) << 4));
    operand = (((operand & (uint32_t)(0xff00ff00)) >> 8) | ((operand & (uint32_t)(0x00ff00ff)) << 8));
    return static_cast<uint32_t>(operand >> 16) | static_cast<uint32_t>(operand << 16);
}

inline size_t CKKS_Encoder::reverse_bits(size_t operand)
{
    return static_cast<size_t>(reverse_bits(static_cast<uint32_t>(operand >> 32))) |
           (static_cast<size_t>(reverse_bits(static_cast<uint32_t>(operand & (size_t)(0xFFFFFFFF)))) << 32);
}

inline size_t CKKS_Encoder::reverse_bits(size_t operand, int bit_count)
{
    return (bit_count == 0) ? size_t(0) : reverse_bits(operand) >> (
        sizeof(size_t) * static_cast<size_t>(bits_per_byte)
            - static_cast<size_t>(bit_count));
}
