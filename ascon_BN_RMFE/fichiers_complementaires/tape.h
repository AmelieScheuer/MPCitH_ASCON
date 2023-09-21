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



#pragma once

extern "C" {
#include "kdf_shake.h"
}
#include "gsl-lite.hpp"
#include "types.h"
#include <cstdlib>

class RandomTape {
private:
  /* data */
  hash_context ctx;

public:
  RandomTape(const gsl::span<uint8_t> &seed, const banquet_salt_t &salt,
             size_t rep_index, size_t party_index);
  ~RandomTape() = default;

  void squeeze_bytes(uint8_t *out, size_t len);
};

class RandomTapes {
private:
  /* data */
  RepByteContainer random_tapes;
  size_t random_tape_size;

public:
  RandomTapes(size_t num_repetitions, size_t num_parties,
              size_t random_tape_size)
      : random_tapes(num_repetitions, num_parties, random_tape_size),
        random_tape_size(random_tape_size){};
  ~RandomTapes() = default;

  void generate_4_tapes(size_t repetition, size_t start_party,
                        const banquet_salt_t &salt,
                        const gsl::span<uint8_t> &seed0,
                        const gsl::span<uint8_t> &seed1,
                        const gsl::span<uint8_t> &seed2,
                        const gsl::span<uint8_t> &seed3);
  void generate_tape(size_t repetition, size_t party,
                     const banquet_salt_t &salt,
                     const gsl::span<uint8_t> &seed);
  gsl::span<uint8_t> get_bytes(size_t repetition, size_t party, size_t start,
                               size_t len);
};