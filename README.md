# CKKS-FV-HERA
This is an implementation of the HERA cipher and its application to the CKKS-FV hybrid transciphering framework proposed in [Hybrid Framework for Approximate Computation over Encrypted Data]([https://eprint.iacr.org/2020/1335](https://eprint.iacr.org/archive/2020/1335/20210304:021449)).

## Client-side
[HERA](./HERA) contains a client-side implementation of the hybrid framework:
an AVX2-optimized C++ implementation of the HERA cipher whose XOF is instantiated by SHAKE256 (https://github.com/XKCP/XKCP), and encoders of CKKS and FV schemes to encode complex-number data and HERA keystream.
Parameter presets are given in [parms.h](./HERA/parms.h).

## Server-side
[Hybrid_Framework](./Hybrid_Framework) is a server-side implementation of the hybrid framework.
The hybrid framework is implemented using a modified SEAL library (https://github.com/smilecjf/SEAL) that supports CKKS-FV scheme.
The server-side homomorphically evaluates HERA keystream under FV scheme to recover CKKS ciphertext of original data from the client-side.
Parameter presets are given in [parms.h](./Hybrid_Framework/parms.h).
