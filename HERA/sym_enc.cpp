#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include "FV_Encoder.h"
#include "CKKS_Encoder.h"
#include "Hera.h"
#include "parms.h"

using namespace std;
using namespace chrono;


int main()
{
    system_clock::time_point time_start, time_end;
    uint64_t ntimes = 1000;
    block_t ctxt = block_init(POLY_MOD_DEG);
    block_t encoded_keystream = block_init(POLY_MOD_DEG);

    FV_Encoder fv_encoder;
    fv_encoder.init();
    CKKS_Encoder ckks_encoder;
    ckks_encoder.init();

    // Message generation
    size_t slot_count = POLY_MOD_DEG / 2;
    complex<double> values[slot_count];
    size_t max = 4;
    size_t scale = static_cast<double>(MODULUS) / max;
    for (size_t i = 0; i < slot_count; i++)
    {
        double i_double = static_cast<double> (i % max);
        if (i % 2 == 0)
            i_double *= -1;
        values[i] = complex<double> (i_double, 0);
    }

    block_t encoded_msg = block_init(POLY_MOD_DEG);

    // Key generation
    block_t key = block_init(BLOCKSIZE);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> dis(0, MODULUS - 1);
    for (int i = 0; i < 16; i++)
    {
        key[i] = dis(gen);
    }

    Hera cipher(key);

    uint64_t nonce = 0;
    uint64_t counter = 1;


    cipher.init(nonce, counter);
    time_start = system_clock::now();

    for (size_t j = 0; j < ntimes; j++)
    {
        counter = 1;
        for (size_t i = 0; i < POLY_MOD_DEG / BLOCKSIZE; i++)
        {
            // Keystream generation
            cipher.update(nonce, counter++);
            cipher.crypt(ctxt + BLOCKSIZE * i);
        }

        // FV-encoding of keystream
        fv_encoder.encode(ctxt, encoded_keystream);
    }


    time_end = system_clock::now();
    duration<double> time_offline = (time_end - time_start) / ntimes;


    time_start = system_clock::now();


    for (size_t j = 0; j < ntimes; j++)
    {
        // CKKS-encoding of message
        ckks_encoder.encode(values, scale, encoded_msg, MODULUS);

        for (size_t i = 0; i < POLY_MOD_DEG; i++)
        {
            encoded_msg[i] += encoded_keystream[i];
        }
    }

    time_end = system_clock::now();
    duration<double> time_online = (time_end - time_start) / ntimes;

    cout << "Offline time is " << time_offline.count() << " sec" << endl;
    cout << "Online time is " << time_online.count() << " sec" << endl;

    return 0;
}
