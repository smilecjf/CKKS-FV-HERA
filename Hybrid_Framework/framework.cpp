#include <iostream>
#include <chrono>
#include <random>
#include "seal/seal.h"
#include "seal/util/pointer.h"
#include "seal/util/polycore.h"
#include "seal/util/polyarith.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/scalingvariant.h"
#include "Hera.h"
#include "FV_Encoder.h"
#include "CKKS_Encoder.h"

#define PRINT_RESULT 1
#define PRINT_FULL 0

using namespace std;
using namespace seal;

// Code from examples.h in SEAL library
inline void print_parameters(std::shared_ptr<seal::SEALContext> context)
{
    // Verify parameters
    if (!context)
    {
        throw std::invalid_argument("context is not set");
    }
    auto &context_data = *context->key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::BFV:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::CKKS:
        scheme_name = "CKKS";
        break;
    case seal::scheme_type::CKKS_FV:
        scheme_name = "CKKS-FV";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " <<
        context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_mod_count = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_mod_count - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::BFV)
    {
        std::cout << "|   plain_modulus: " << context_data.
            parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}

Plaintext make_plain_const(uint64_t value, BatchEncoder &batch_encoder, int slot_count)
{
    Plaintext plain;
    vector<uint64_t> data(slot_count, value);
    batch_encoder.encode(data, plain);
    return plain;
}

void linear_layer(Ciphertext *cipher, vector<uint64_t*> buf, vector<SmallModulus> coeff_modulus, size_t coeff_count)
{
    size_t coeff_mod_count = coeff_modulus.size();

    // Check size of ciphertexts
    for (int l = 0; l < BLOCKSIZE; l++)
    {
        if (cipher[l].size() != 2)
        {
            throw invalid_argument("invalid ciphertext size");
        }
    }

    for (int j = 0; j < 2; j++)
    {
        // MixColumns matrix multiplication
        for (int col = 0; col < 4; col++)
        {
            uint64_t *cipher0_ptr = cipher[col].data(j);
            uint64_t *cipher1_ptr = cipher[col + 4].data(j);
            uint64_t *cipher2_ptr = cipher[col + 8].data(j);
            uint64_t *cipher3_ptr = cipher[col + 12].data(j);

            uint64_t *buf0_ptr = buf[col];
            uint64_t *buf1_ptr = buf[col + 4];
            uint64_t *buf2_ptr = buf[col + 8];
            uint64_t *buf3_ptr = buf[col + 12];

            // Compute addition in MixColumns
            size_t total_coeff_count = coeff_count * coeff_mod_count;
            for (size_t i = 0; i < total_coeff_count; i++)
            {
                *(buf0_ptr + i) = (*(cipher0_ptr + i) << 1) + (*(cipher1_ptr + i) << 1) +
                    *(cipher1_ptr + i) + *(cipher2_ptr + i) + *(cipher3_ptr + i);
                *(buf1_ptr + i) = (*(cipher1_ptr + i) << 1) + (*(cipher2_ptr + i) << 1) +
                    *(cipher2_ptr + i) + *(cipher3_ptr + i) + *(cipher0_ptr + i);
                *(buf2_ptr + i) = (*(cipher2_ptr + i) << 1) + (*(cipher3_ptr + i) << 1) +
                    *(cipher3_ptr + i) + *(cipher0_ptr + i) + *(cipher1_ptr + i);
                *(buf3_ptr + i) = (*(cipher3_ptr + i) << 1) + (*(cipher0_ptr + i) << 1) +
                    *(cipher0_ptr + i) + *(cipher1_ptr + i) + *(cipher2_ptr + i);
            }
        }

        // MixRows matrix multiplication
        for (int row = 0; row < 4; row++)
        {
            uint64_t *buf0_ptr = buf[4 * row];
            uint64_t *buf1_ptr = buf[4 * row + 1];
            uint64_t *buf2_ptr = buf[4 * row + 2];
            uint64_t *buf3_ptr = buf[4 * row + 3];

            uint64_t *cipher0_ptr = cipher[4 * row].data(j);
            uint64_t *cipher1_ptr = cipher[4 * row + 1].data(j);
            uint64_t *cipher2_ptr = cipher[4 * row + 2].data(j);
            uint64_t *cipher3_ptr = cipher[4 * row + 3].data(j);

            // Compute addition in MixRows
            size_t total_coeff_count = coeff_count * coeff_mod_count;
            for (size_t i = 0; i < total_coeff_count; i++)
            {
                *(cipher0_ptr + i) = (*(buf0_ptr + i) << 1) + (*(buf1_ptr + i) << 1) +
                    *(buf1_ptr + i) + *(buf2_ptr + i) + *(buf3_ptr + i);
                *(cipher1_ptr + i) = (*(buf1_ptr + i) << 1) + (*(buf2_ptr + i) << 1) +
                    *(buf2_ptr + i) + *(buf3_ptr + i) + *(buf0_ptr + i);
                *(cipher2_ptr + i) = (*(buf2_ptr + i) << 1) + (*(buf3_ptr + i) << 1) +
                    *(buf3_ptr + i) + *(buf0_ptr + i) + *(buf1_ptr + i);
                *(cipher3_ptr + i) = (*(buf3_ptr + i) << 1) + (*(buf0_ptr + i) << 1) +
                    *(buf0_ptr + i) + *(buf1_ptr + i) + *(buf2_ptr + i);
            }
        }

        // Modulo operation
        for (int l = 0; l < BLOCKSIZE; l++)
        {
            for (int i = 0; i < coeff_mod_count; i++)
            {
                util::modulo_poly_coeffs_63(
                        cipher[l].data(j) + i * coeff_count,
                        coeff_count, coeff_modulus[i],
                        cipher[l].data(j) + i * coeff_count);
            }
        }
    }
}

int main()
{
    // Set Timer
    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_rks_gen, time_rks_ecd, \
                        time_ark, time_linear, time_cube, \
                        time_postproc_first, time_postproc_last, \
                        time_he_sub_first, time_he_sub_last, \
                        time_client_ckks, time_client_fv, time_client_enc;
    time_ark = chrono::microseconds(0);
    time_linear = chrono::microseconds(0);
    time_cube = chrono::microseconds(0);

    // Set random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> uniform_dist(0, MODULUS - 1);

    /* Sever-side Computation */
    // Set HE parameter
    sec_level_type sec_level = sec_level_type::tc128;
    EncryptionParameters parms(scheme_type::CKKS_FV);
    size_t poly_modulus_degree = POLY_MOD_DEG;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    auto coeff_modulus = CoeffModulus::BFVDefault(poly_modulus_degree, sec_level);
    for (int i = 0; i < INIT_MOD_DROP; i++)
    {
        coeff_modulus.pop_back();
    }
    parms.set_coeff_modulus(coeff_modulus);
    size_t plain_modulus = MODULUS;
    parms.set_plain_modulus(plain_modulus);

    auto context = SEALContext::Create(parms, true, sec_level);
    auto first_parms_id = context->first_parms_id();
    print_parameters(context);
    cout << endl;

    auto first_context_data = context->first_context_data();
    auto qualifiers = first_context_data->qualifiers();
    cout << "Batching enabled: " << boolalpha << qualifiers.using_batching << endl;

    // Init HE context
    KeyGenerator keygen(context);
    PublicKey public_key = keygen.public_key();
    SecretKey secret_key = keygen.secret_key();
    RelinKeys relin_keys = keygen.relin_keys();
    GaloisKeys gal_keys = keygen.galois_keys();
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    BatchEncoder batch_encoder(context);
    CKKSEncoder ckks_encoder(context);
    size_t slot_count = batch_encoder.slot_count();
    cout << "Plaintext slot count: " << slot_count << endl;

    // HERK parameters
    int rounds = ROUNDS;
    uint64_t nonce = 0x0123456789abcedf;
    uint64_t counter = 0;
    uint64_t key[BLOCKSIZE];
    cout << "Symmetric-key generation" << endl;
    cout << "  key: ";
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        key[i] = uniform_dist(gen);
        cout << key[i] << " ";
    }
    cout << endl;

    Hera cipher(key);
    cipher.init(nonce, counter);
    for (int i = 0; i < 20000; i++)
    {
        cipher.update(nonce, counter);
    }

    Plaintext p_key[BLOCKSIZE];
    Ciphertext c[BLOCKSIZE], c_key[BLOCKSIZE], c_temp;

    // Encrypt key
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        p_key[i] = make_plain_const(key[i], batch_encoder, slot_count);
        encryptor.encrypt(p_key[i], c_key[i]);
    }

    // Key randomizing source
    uint64_t temp_xof[(rounds + 1) * BLOCKSIZE];
    block_t xof = block_init(BLOCKSIZE * (rounds + 1) * poly_modulus_degree);
    size_t N_div_blocksize = poly_modulus_degree / BLOCKSIZE;
    time_start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < poly_modulus_degree; i++)
    {
        cipher.update(nonce, i);
        cipher.get_rand_vectors(temp_xof);
        for (int r = 0; r <= rounds; r++)
        {
            for (size_t j = 0; j < BLOCKSIZE; j++)
            {
                size_t offset = poly_modulus_degree * (r * BLOCKSIZE + j);
                size_t index = (i % N_div_blocksize) * BLOCKSIZE + (i / N_div_blocksize);
                xof[offset + index] = temp_xof[r * BLOCKSIZE + j];
            }
        }
    }
    time_end = chrono::high_resolution_clock::now();
    time_rks_gen = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

    Plaintext encoded_rand_vectors[(rounds + 1) * BLOCKSIZE];
    time_start = chrono::high_resolution_clock::now();
    for (int r = 0; r <= rounds; r++)
    {
        for (int j = 0; j < BLOCKSIZE; j++)
        {
            uint64_t *data = xof + poly_modulus_degree * (r * BLOCKSIZE + j);
            vector<uint64_t> data_vec(poly_modulus_degree);
            for (int i = 0; i < poly_modulus_degree; i++)
            {
                data_vec[i] = data[i];
            }
            batch_encoder.encode(data_vec, encoded_rand_vectors[r * BLOCKSIZE + j]);
        }
    }
    time_end = chrono::high_resolution_clock::now();
    time_rks_ecd = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

    // Set input constant ic = (1, 2, ..., 16) for evaluation input
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        c[i].resize(context, first_parms_id, 2);
        c[i].is_ntt_form() = false;
        c[i].parms_id() = first_parms_id;

        Plaintext p_ic;
        vector<uint64_t> data_ic(slot_count);
        for (size_t j = 0; j < slot_count; j++)
        {
            data_ic[j] = i + 1;
        }
        batch_encoder.encode(data_ic, p_ic);
        util::multiply_add_plain_with_scaling_variant(
            p_ic, *context->get_context_data(first_parms_id), c[i].data());
    }

    // Eval
    cout << "Eval..." << endl;
    Ciphertext c_buf[BLOCKSIZE];
    auto curr_coeff_modulus = context->get_context_data(first_parms_id)->parms().coeff_modulus();
    size_t coeff_mod_count = curr_coeff_modulus.size();
    util::Pointer<uint64_t> poly_buf[BLOCKSIZE];
    vector<uint64_t*> poly_buf_ptr;
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        poly_buf[i] = util::allocate_poly(poly_modulus_degree, coeff_mod_count, MemoryManager::GetPool());
        poly_buf_ptr.push_back(poly_buf[i].get());
    }
    for (int r = 0; r < rounds; r++)
    {
        cout << "  Round " << r + 1 << endl;
        // Add round key
        cout << "  | Add round key " << endl;
        time_start = chrono::high_resolution_clock::now();
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            evaluator.multiply_plain(c_key[i], encoded_rand_vectors[r * BLOCKSIZE + i], c_temp);
            evaluator.add_inplace(c[i], c_temp);
        }
        time_end = chrono::high_resolution_clock::now();
        time_ark += chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "  - Ciphertext noise budget: " << decryptor.invariant_noise_budget(c[0]) << " bits" << endl;

        // Linear layer
        cout << "  | Linear layer" << endl;
        time_start = chrono::high_resolution_clock::now();
        linear_layer(c, poly_buf_ptr, curr_coeff_modulus, poly_modulus_degree);
        time_end = chrono::high_resolution_clock::now();
        time_linear += chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "  - Ciphertext noise budget: " << decryptor.invariant_noise_budget(c[0]) << " bits" << endl;

        // Cube
        cout << "  | Cube " << endl;
        time_start = chrono::high_resolution_clock::now();
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            evaluator.square(c[i], c_temp);
            evaluator.relinearize_inplace(c_temp, relin_keys);
            evaluator.multiply_inplace(c[i], c_temp);
            evaluator.relinearize_inplace(c[i], relin_keys);
        }
        time_end = chrono::high_resolution_clock::now();
        time_cube += chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "  - Ciphertext noise budget: " << decryptor.invariant_noise_budget(c[0]) << " bits" << endl;
    }

    // Finalization
    cout << "  Finalization" << endl;
    cout << "  | MixColumns" << endl;
    time_start = chrono::high_resolution_clock::now();
    linear_layer(c, poly_buf_ptr, curr_coeff_modulus, poly_modulus_degree);
    time_end = chrono::high_resolution_clock::now();
    time_linear += chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "  - Ciphertext noise budget: " << decryptor.invariant_noise_budget(c[0]) << " bits" << endl;

    cout << "  | Add round key " << endl;
    time_start = chrono::high_resolution_clock::now();
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        evaluator.multiply_plain(c_key[i], encoded_rand_vectors[rounds * BLOCKSIZE + i], c_temp);
        evaluator.add_inplace(c[i], c_temp);
    }
    time_end = chrono::high_resolution_clock::now();
    time_ark += chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "  - Ciphertext noise budget: " << decryptor.invariant_noise_budget(c[0]) << " bits" << endl;

    // Post-process
    cout << "Post-process..." << endl;
    Ciphertext c_eval[BLOCKSIZE];
    time_start = chrono::high_resolution_clock::now();
    c_eval[0] = Ciphertext(c[0]);
    vector<uint64_t> mask(slot_count, 0ULL);
    for (int i = 0; i < slot_count / BLOCKSIZE; i++)
    {
        mask[BLOCKSIZE * i] = 1;
    }
    Plaintext plain_mask;
    batch_encoder.encode(mask, plain_mask);
    evaluator.multiply_plain_inplace(c_eval[0], plain_mask);
    for (int i = 1; i < BLOCKSIZE; i++)
    {
        evaluator.multiply_plain(c[i], plain_mask, c_temp);
        evaluator.rotate_rows_inplace(c_temp, slot_count / 2 - i, gal_keys);
        evaluator.add_inplace(c_eval[0], c_temp);
    }
    time_end = chrono::high_resolution_clock::now();
    time_postproc_first = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "  - Ciphertext noise budget: " << decryptor.invariant_noise_budget(c_eval[0]) << " bits" << endl;

    time_start = chrono::high_resolution_clock::now();
    for (int l = 1; l < BLOCKSIZE; l++)
    {
        for (int i = 0; i < slot_count / BLOCKSIZE; i++)
        {
            mask[BLOCKSIZE * i + l - 1] = 0;
            mask[BLOCKSIZE * i + l] = 1;
        }
        batch_encoder.encode(mask, plain_mask);
        evaluator.multiply_plain(c[0], plain_mask, c_eval[l]);
        evaluator.rotate_rows_inplace(c_eval[l], l, gal_keys);
        for (int i = 1; i < BLOCKSIZE; i++)
        {
            evaluator.multiply_plain(c[i], plain_mask, c_temp);
            if (l > i)
            {
                evaluator.rotate_rows_inplace(c_temp, l - i, gal_keys);
            }
            else if (l < i)
            {
                evaluator.rotate_rows_inplace(c_temp, slot_count / 2 + l - i, gal_keys);
            }
            evaluator.add_inplace(c_eval[l], c_temp);
        }
    }
    time_end = chrono::high_resolution_clock::now();
    time_postproc_last = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Finished" << endl;

    /* Client-side Computation */
    cout << "Client-side..." << endl;

    // CKKS encoding of client data
    cout << "CKKS encoding..." << endl;
    CKKS_Encoder client_ckks_encoder;
    client_ckks_encoder.init();

    size_t ckks_slot_count = poly_modulus_degree / 2;
    size_t max = 4;
    double precision_scale = static_cast<double>(MODULUS) / static_cast<double>(2 * max);

    complex<double> *client_data = (complex<double>*)malloc(sizeof(complex<double>) * BLOCKSIZE * ckks_slot_count);
    uint64_t client_encrypted_data[BLOCKSIZE * poly_modulus_degree];

    for (int l = 0; l < BLOCKSIZE; l++)
    {
        for (size_t i = 0; i < ckks_slot_count; i++)
        {
            double real = static_cast<double>(i % max + 1);
            complex<double> val = complex<double>(real, 0);
            client_data[l * ckks_slot_count + i] = val;
        }
    }

    time_start = chrono::high_resolution_clock::now();
    for (int l = 0; l < BLOCKSIZE; l++)
    {
        client_ckks_encoder.encode(client_data + l * ckks_slot_count, precision_scale,
                client_encrypted_data + l * poly_modulus_degree, MODULUS);
    }
    time_end = chrono::high_resolution_clock::now();
    time_client_ckks = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

    // Compute FV encoding of Hera keystream
    cout << "FV encoding..." << endl;
    FV_Encoder client_fv_encoder;
    client_fv_encoder.init();

    block_t keystream = block_init(BLOCKSIZE * poly_modulus_degree);
    block_t encoded_keystream = block_init(BLOCKSIZE * poly_modulus_degree);

    uint64_t temp_counter = 0;
    time_start = chrono::high_resolution_clock::now();
    for (int l = 0; l < BLOCKSIZE; l++)
    {
        uint64_t offset = l * poly_modulus_degree;
        for (size_t i = 0; i < poly_modulus_degree / BLOCKSIZE; i++)
        {
            cipher.update(nonce, temp_counter);
            temp_counter++;
            cipher.crypt(keystream + i * BLOCKSIZE + offset);
        }

        client_fv_encoder.encode(keystream + offset, encoded_keystream + offset);
    }
    time_end = chrono::high_resolution_clock::now();
    time_client_fv = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

    time_start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < BLOCKSIZE * poly_modulus_degree; i++)
    {
        client_encrypted_data[i] += encoded_keystream[i];
    }
    time_end = chrono::high_resolution_clock::now();
    time_client_enc = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Finished" << endl;

    /* Homomorphic Decryption */
    cout << "Homomorphic subtraction..." << endl;
    Plaintext p_client[BLOCKSIZE];
    Ciphertext c_client[BLOCKSIZE];

    time_start = chrono::high_resolution_clock::now();
    // Set plaintext
    p_client[0].resize(poly_modulus_degree);
    p_client[0].parms_id() = parms_id_zero;
    memcpy(p_client[0].data(), client_encrypted_data,
            sizeof(uint64_t) * poly_modulus_degree);

    // Trivial encryption
    auto parms_id = c_eval[0].parms_id();
    auto context_data = context->get_context_data(parms_id);
    c_client[0].resize(context, parms_id, 2);
    c_client[0].is_ntt_form() = false;
    c_client[0].parms_id() = parms_id;
    util::multiply_add_plain_with_scaling_variant(
            p_client[0], *context_data, c_client[0].data());

    // Subtract homomorphically evaluated keystream
    evaluator.sub_inplace(c_client[0], c_eval[0]);
    time_end = chrono::high_resolution_clock::now();
    time_he_sub_first = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "  - Ciphertext noise budget: " << decryptor.invariant_noise_budget(c_client[0]) << " bits" << endl;

    time_start = chrono::high_resolution_clock::now();
    for (int l = 1; l < BLOCKSIZE; l++)
    {
        // Set plaintext
        uint64_t offset = l * poly_modulus_degree;
        p_client[l].resize(poly_modulus_degree);
        p_client[l].parms_id() = parms_id_zero;
        memcpy(p_client[l].data(), client_encrypted_data + offset,
                sizeof(uint64_t) * poly_modulus_degree);

        // Trivial encryption
        c_client[l].resize(context, parms_id, 2);
        c_client[l].is_ntt_form() = false;
        c_client[l].parms_id() = parms_id;
        util::multiply_add_plain_with_scaling_variant(
                p_client[l], *context_data, c_client[l].data());

        // Subtract homomorphically evaluated keystream
        evaluator.sub_inplace(c_client[l], c_eval[l]);
    }
    time_end = chrono::high_resolution_clock::now();
    time_he_sub_last = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Finished" << endl;

    cout << "=============== Time [\u03bcs] ===============" << endl;
    cout << "0. RKS" << endl;
    cout << "| Time RKS gen    : " << time_rks_gen.count() << endl;
    cout << "| Time RKS ecd    : " << time_rks_ecd.count() << endl;
    cout << "1. HE Eval" << endl;
    cout << "| Time ARK        : " << time_ark.count() << endl;
    cout << "| Time Linear     : " << time_linear.count() << endl;
    cout << "| Time Cube       : " << time_cube.count() << endl;
    cout << "| First Post-proc : " << time_postproc_first.count() << endl;
    cout << "| Last Post-proc  : " << time_postproc_last.count() << endl;
    cout << "2. Client Side" << endl;
    cout << "| Time CKKS Ecd   : " << time_client_ckks.count() << endl;
    cout << "| Time FV Ecd     : " << time_client_fv.count() << endl;
    cout << "| Time Enc        : " << time_client_enc.count() << endl;
    cout << "3. Homomorphic Subtraction" << endl;
    cout << "| First HE Sub    : " << time_he_sub_first.count() << endl;
    cout << "| Last HE Sub     : " << time_he_sub_last.count() << endl;
    cout << "=========================================" << endl;

#if PRINT_RESULT
    Plaintext p_result[BLOCKSIZE];
    vector<vector<complex<double>>> decrypted_data(BLOCKSIZE);
    for (int l = 0; l < BLOCKSIZE; l++)
    {
        evaluator.transform_to_lowest_ckks_inplace(c_client[l]);
        decryptor.decrypt(c_client[l], p_result[l]);
        ckks_encoder.decode(p_result[l], decrypted_data[l]);
    }

    for (int l = 0; l < BLOCKSIZE; l++)
    {
        cout << "Result [" << l << "]" << endl;
#if PRINT_FULL
        for (size_t i = 0; i < ckks_slot_count; i++)
#else
        for (size_t i = 0; i < 64; i++)
#endif
        {
            double real = decrypted_data[l][i].real() / precision_scale;
            double imag = decrypted_data[l][i].imag() / precision_scale;
            cout << "  " << i << ": (" << real << ", " << imag << ")" << endl;
        }
    }
#endif

    free(client_data);
    return 0;
}

