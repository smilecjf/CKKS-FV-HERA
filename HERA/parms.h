#ifndef __PARMS_H__
#define __PARMS_H__

#define TRUE 1
#define FALSE 0

// BATCHABLE PRIMES
#define PRIME17 65537 // 2^16 + 1
#define PRIME22 3604481
#define PRIME23 7340033
#define PRIME24 16580609
#define PRIME25 33292289
#define PRIME26 65929217
#define PRIME27 130809857
#define PRIME28 268238849
#define PRIME29 536608769
#define PRIME30 1069219841
#define PRIME31 2146041857
#define PRIME32 4292804609

// Parameter Preset
#define PARM_PRESET 1
#define BLOCKSIZE 16

// Par-I
// 80-bit security, 10-bit precision
#if PARM_PRESET == 1
    #define POLY_MOD_DEG 32768
    #define ROUNDS 4
    #define MODULUS PRIME28
    #define MOD_BIT_COUNT 28
    #define MONTGOMERY FALSE
    #define MONT_MOD_BIT 32
    #define MONT_MOD_INV 16752649
    #define MONT_PSEUDO_INV 268238847
    #define MONT_ONE 3145712
    #define ALL_ONE_MOD_BIT ((1ULL << MONT_MOD_BIT) - 1)
// Par-II
// 80-bit security, 14-bit precision
#elif PARM_PRESET == 2
    #define POLY_MOD_DEG 32768
    #define ROUNDS 4
    #define MODULUS PRIME32
    #define MOD_BIT_COUNT 32
    #define MONTGOMERY FALSE
// Par-III
// 128-bit security, 10-bit precision
#elif PARM_PRESET == 3
    #define POLY_MOD_DEG 32768
    #define ROUNDS 5
    #define MODULUS PRIME28
    #define MOD_BIT_COUNT 28
    #define MONTGOMERY TRUE
    #define MONT_MOD_BIT 32
    #define MONT_MOD_INV 16752649
    #define MONT_PSEUDO_INV 268238847
    #define MONT_ONE 3145712
    #define ALL_ONE_MOD_BIT ((1ULL << MONT_MOD_BIT) - 1)
// Par-IV
// 128-bit security, 14-bit precision
#elif PARM_PRESET == 4
    #define POLY_MOD_DEG 32768
    #define ROUNDS 5
    #define MODULUS PRIME32
    #define MOD_BIT_COUNT 32
    #define MONTGOMERY FALSE
// Par-A
// 128-bit security, 9-bit precision
#elif PARM_PRESET == 5
    #define POLY_MOD_DEG 16
    #define ROUNDS 5
    #define MODULUS PRIME17
    #define MOD_BIT_COUNT 17
    #define NUM_SQUEEZE 1
    #define MONTGOMERY FALSE
    #define MONT_MOD_BIT 32
    #define MONT_MOD_INV 1
    #define MONT_PSEUDO_INV 65535
    #define MONT_ONE 1
    #define ALL_ONE_MOD_BIT ((1ULL << MONT_MOD_BIT) - 1)
// Par-B
// 128-bit security, 13-bit precision
#elif PARM_PRESET == 6
    #define POLY_MOD_DEG 64
    #define ROUNDS 5
    #define MODULUS PRIME22
    #define MOD_BIT_COUNT 22
    #define MONTGOMERY TRUE
    #define MONT_MOD_BIT 32
    #define MONT_MOD_INV 3025
    #define MONT_PSEUDO_INV 3604479
    #define MONT_ONE 2030425
    #define ALL_ONE_MOD_BIT ((1ULL << MONT_MOD_BIT) - 1)
// Par-C
// 128-bit security, 11-bit precision
#elif PARM_PRESET == 7
    #define POLY_MOD_DEG 256
    #define ROUNDS 5
    #define MODULUS PRIME22
    #define MOD_BIT_COUNT 22
    #define MONTGOMERY TRUE
    #define MONT_MOD_BIT 32
    #define MONT_MOD_INV 3025
    #define MONT_PSEUDO_INV 3604479
    #define MONT_ONE 2030425
    #define ALL_ONE_MOD_BIT ((1ULL << MONT_MOD_BIT) - 1)
// Par-D
// 128-bit security, 10-bit precision
#elif PARM_PRESET == 8
    #define POLY_MOD_DEG 1024
    #define ROUNDS 5
    #define MODULUS PRIME23
    #define MOD_BIT_COUNT 23
    #define MONTGOMERY TRUE
    #define MONT_MOD_BIT 32
    #define MONT_MOD_INV 12544
    #define MONT_PSEUDO_INV 7340031
    #define MONT_ONE 1047991
    #define ALL_ONE_MOD_BIT ((1ULL << MONT_MOD_BIT) - 1)
#endif

#define XOF_ELEMENT_COUNT ((ROUNDS + 1) * BLOCKSIZE)
#endif
