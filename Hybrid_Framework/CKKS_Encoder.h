/*
    Code from https://github.com/microsoft/SEAL.git (version 3.4.5)
*/
#ifndef __CKKS_ENCODER_H__
#define __CKKS_ENCODER_H__

#include <vector>
#include <complex>
#include "parms.h"

using namespace std;

constexpr size_t root_fidelity_ = 131072;
constexpr double PI_ = 3.1415926535897932384626433832795028842;
constexpr int bits_per_byte = 8;
constexpr size_t poly_modulus_degree = POLY_MOD_DEG;

class CKKS_Encoder
{
public:
    CKKS_Encoder() {};
    void init();
    void encode(const complex<double> *values, double scale,
        uint64_t *destination, size_t plain_modulus);

    size_t get_poly_modulus_degree()
    {
        return poly_modulus_degree;
    }

private:
    // Private data
    size_t slot_count_;
    size_t slot_count_power_of_two_;
    size_t coeff_count_power_of_two_;
    complex<double> roots_[poly_modulus_degree];
    complex<double> inv_roots_[poly_modulus_degree];
    size_t matrix_reps_index_map_[poly_modulus_degree];

    // Initialization functions
    complex<double> get_root(size_t index);
    inline complex<double> get_root(std::size_t index, std::size_t of_roots);
    void populate_reps_index_map();

    // Utility functions
    int get_power_of_two(size_t value);
    inline uint32_t reverse_bits(uint32_t operand);
    inline size_t reverse_bits(size_t operand);
    inline size_t reverse_bits(size_t operand, int bit_count);
};

#endif
