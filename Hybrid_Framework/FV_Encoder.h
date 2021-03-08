/*
    Code from https://github.com/microsoft/SEAL.git (version 3.4.5)
*/
#ifndef __FV_ENCODER_H__
#define __FV_ENCODER_H__

#include <cinttypes>
#include <vector>
#include "parms.h"

constexpr size_t poly_modulus_degree_fv = POLY_MOD_DEG;
constexpr size_t plain_modulus = MODULUS;

using namespace std;

class FV_Encoder
{
public:
    void init();
    vector<uint64_t> encode(vector<uint64_t> data);
    void encode(uint64_t* data, uint64_t* destination);

    uint64_t get_poly_modulus_degree()
    {
        return poly_modulus_degree_fv;
    }
    uint64_t get_plain_modulus()
    {
        return plain_modulus;
    }

private:
    // FV_Encoder Private Data
    uint64_t coeff_count_power_;
    uint64_t root_;
    uint64_t inverse_root_;
    uint64_t root_powers_[poly_modulus_degree_fv];
    uint64_t inv_root_powers_[poly_modulus_degree_fv];
    uint64_t roots_of_unity_[poly_modulus_degree_fv];
    uint64_t inv_root_powers_div_two_[poly_modulus_degree_fv];
    uint64_t scaled_inv_root_powers_div_two_[poly_modulus_degree_fv];
    uint64_t scaled_root_powers_[poly_modulus_degree_fv];
    uint64_t scaled_inv_root_powers_[poly_modulus_degree_fv];
    uint64_t matrix_reps_index_map_[poly_modulus_degree_fv];

    // FV_Encoder initialization functions
    void populate_matrix_reps_index_map();
    void populate_roots_of_unity_vector();

    // FV_Encoder utility functions
    void inverse_ntt_negacyclic_harvey_lazy(uint64_t *operand);
    inline void inverse_ntt_negacyclic_harvey(uint64_t *operand);
    void ntt_powers_of_primitive_root(uint64_t root, uint64_t *destination);
    void ntt_scale_powers_of_primitive_root(const uint64_t *input, uint64_t *destination);
};

#endif
