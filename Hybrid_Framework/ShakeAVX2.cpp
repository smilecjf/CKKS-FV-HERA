#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include "ShakeAVX2.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

using namespace std;

void ShakeAVX2::update(uint64_t nonce, uint64_t counter)
{
    nonce_ = nonce;
    counter_ = counter;

    KeccakP1600times4_InitializeAll((void*)states_);
    instance_index_ = 0;
    instance_offset_ = 0;


    uint8_t data[48];
    memset(data, 0, sizeof(data));
    for (int i = 0; i < 8; i++)
    {
        data[i] = (rounds_ >> (8 * i)) & 0xff;
        data[i + 24] = (nonce >> (8 * i)) & 0xff;
        data[i + 40] = (counter >> (8 * i)) & 0xff;
        data[32] = 0x1f;
    }
    KeccakP1600times4_AddLanesAll((void*)states_, data, 6, 0);

    uint8_t instance_seed[8];
    for (int index = 0; index < 4; index++)
    {
        for (int i = 0; i < 8; i++)
        {
            instance_seed[i] = index + 1;
        }
        KeccakP1600times4_AddBytes((void*)states_, index, instance_seed, 8, 8);
        KeccakP1600times4_AddByte((void*)states_, index, 0x80, RATE_IN_BYTE);
    }

    // Permute 
    KeccakP1600times4_PermuteAll_24rounds((void*)states_);
}

void ShakeAVX2::squeeze(uint8_t *output, unsigned int length)
{
    unsigned offset = 0;
    while (length > 0)
    {
        if (instance_index_ == 3 && instance_offset_ == RATE_IN_BYTE)
        {
            KeccakP1600times4_PermuteAll_24rounds((void*)states_);
            instance_index_ = 0;
            instance_offset_ = 0;
        }
        else if (instance_offset_ == RATE_IN_BYTE)
        {
            instance_index_++;
            instance_offset_ = 0;
        }

        unsigned int squeeze_size = MIN(length, RATE_IN_BYTE - instance_offset_);
        KeccakP1600times4_ExtractBytes((void*)states_, instance_index_, output + offset, instance_offset_, squeeze_size);
        offset += squeeze_size;
        length -= squeeze_size;
        instance_offset_ += squeeze_size;
    }
}

void ShakeAVX2::print_state()
{
    for (int index = 0; index < 4; index++)
    {
        cout << "State " << index << ": " << flush;
        for (int i = 0; i < 200; i++)
        {
            uint8_t data = 0;
            KeccakP1600times4_ExtractBytes((void*)states_, index, &data, i, 1);
            cout << hex << setw(2) << setfill('0') << (unsigned int)data << dec << " " << flush;
        }
        cout << endl;
    }
}

