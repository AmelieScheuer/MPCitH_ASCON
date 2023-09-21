/******************************************************************************
* Copyright (c) 2023 Thales
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* *
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
* *
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
-
******************************************************************************/



#include <cstdlib>
#include "banquet_instances.h"
#include "tape.h"


#ifndef GEN_ALEA_H
#define GEN_ALEA_H



std::pair<banquet_salt_t, std::vector<std::vector<uint8_t>>>
generate_salt_and_seeds(const banquet_instance_t &instance,
                        const banquet_keypair_t &keypair,
                        const uint8_t *message, size_t message_len) {
  // salt, seed_1, ..., seed_r = H(instance||sk||pk||m)
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update_uint16_le(&ctx, (uint16_t)instance.params);
  hash_update(&ctx, keypair.first.data(), keypair.first.size());
  hash_update(&ctx, keypair.second.data(), keypair.second.size());
  hash_update(&ctx, message, message_len);
  hash_final(&ctx);

  banquet_salt_t salt;
  hash_squeeze(&ctx, salt.data(), salt.size());
  std::vector<std::vector<uint8_t>> seeds;
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    std::vector<uint8_t> s(instance.seed_size);
    hash_squeeze(&ctx, s.data(), s.size());
    seeds.push_back(s);
  }
  return std::make_pair(salt, seeds);
}

RandomTape::RandomTape(const gsl::span<uint8_t> &seed,
                       const banquet_salt_t &salt, size_t rep_index,
                       size_t party_index) {
  hash_init(&ctx, seed.size() * 2);
  hash_update(&ctx, seed.data(), seed.size());
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, (uint16_t)rep_index);
  hash_update_uint16_le(&ctx, (uint16_t)party_index);
  hash_final(&ctx);
}

void RandomTape::squeeze_bytes(uint8_t *out, size_t len) {
  hash_squeeze(&ctx, out, len);
}

void RandomTapes::generate_tape(size_t repetition, size_t party,
                                const banquet_salt_t &salt,
                                const gsl::span<uint8_t> &seed) {
  hash_context ctx;
  hash_init(&ctx, seed.size() * 2);
  hash_update(&ctx, seed.data(), seed.size());
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, (uint16_t)repetition);
  hash_update_uint16_le(&ctx, (uint16_t)party);
  hash_final(&ctx);
  hash_squeeze(&ctx, random_tapes.get(repetition, party).data(),
               random_tape_size);
}

void RandomTapes::generate_4_tapes(size_t repetition, size_t start_party,
                                   const banquet_salt_t &salt,
                                   const gsl::span<uint8_t> &seed0,
                                   const gsl::span<uint8_t> &seed1,
                                   const gsl::span<uint8_t> &seed2,
                                   const gsl::span<uint8_t> &seed3) {

  hash_context_x4 ctx;
  hash_init_x4(&ctx, seed0.size() * 2);
  hash_update_x4_4(&ctx, seed0.data(), seed1.data(), seed2.data(), seed3.data(),
                   seed0.size());
  hash_update_x4_1(&ctx, salt.data(), salt.size());
  hash_update_x4_uint16_le(&ctx, (uint16_t)repetition);
  const uint16_t parties[4] = {
      (uint16_t)(start_party), (uint16_t)(start_party + 1),
      (uint16_t)(start_party + 2), (uint16_t)(start_party + 3)};
  hash_update_x4_uint16s_le(&ctx, parties);
  hash_final_x4(&ctx);
  hash_squeeze_x4_4(&ctx, random_tapes.get(repetition, parties[0]).data(),
                    random_tapes.get(repetition, parties[1]).data(),
                    random_tapes.get(repetition, parties[2]).data(),
                    random_tapes.get(repetition, parties[3]).data(),
                    random_tape_size);
}

gsl::span<uint8_t> RandomTapes::get_bytes(size_t repetition, size_t party,
                                          size_t start, size_t len) {
  auto tape = random_tapes.get(repetition, party);
  return tape.subspan(start, len);
}

#endif
