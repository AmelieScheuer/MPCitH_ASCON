/*
 *  This file is part of the optimized implementation of the Picnic signature
 *  scheme.
 *  See the accompanying documentation for complete details.
 *  The code is provided under the MIT license:
 *
 * Copyright (c) 2019-2020 Sebastian Ramacher, AIT
 * Copyright (c) 2016-2020 Graz University of Technology
 * Copyright (c) 2017 Angela Promitzer

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the ""Software""), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef KDF_SHAKE_H
#define KDF_SHAKE_H

#include <stdint.h>

#include "macros.h"
#include "portable_endian.h"

#include "keccak/header.h"
// #include "keccak/KeccakHash.h"
// #include "keccak/KeccakHashtimes4.h"

typedef Keccak_HashInstance hash_context ATTR_ALIGNED(32);

/**
 * Initialize hash context based on the digest size used by Picnic. If the size
 * is 32 bytes, SHAKE128 is used, otherwise SHAKE256 is used.
 */
static inline void hash_init(hash_context *ctx, size_t digest_size) {
  if (digest_size == 32) {
    Keccak_HashInitialize_SHAKE128(ctx);
  } else {
    Keccak_HashInitialize_SHAKE256(ctx);
  }
}

static inline void hash_update(hash_context *ctx, const uint8_t *data,
                               size_t size) {
  Keccak_HashUpdate(ctx, data, size << 3);
}

static inline void hash_final(hash_context *ctx) {
  Keccak_HashFinal(ctx, NULL);
}

static inline void hash_squeeze(hash_context *ctx, uint8_t *buffer,
                                size_t buflen) {
  Keccak_HashSqueeze(ctx, buffer, buflen << 3);
}

#define hash_clear(ctx)

static inline void hash_update_uint16_le(hash_context *ctx, uint16_t data) {
  const uint16_t data_le = htole16(data);
  hash_update(ctx, (const uint8_t *)&data_le, sizeof(data_le));
}

static inline void hash_init_prefix(hash_context *ctx, size_t digest_size,
                                    const uint8_t prefix) {
  hash_init(ctx, digest_size);
  hash_update(ctx, &prefix, sizeof(prefix));
}

typedef hash_context kdf_shake_t;

#define kdf_shake_init(ctx, digest_size) hash_init((ctx), (digest_size))
#define kdf_shake_init_prefix(ctx, digest_size, prefix)                        \
  hash_init_prefix((ctx), (digest_size), (prefix))
#define kdf_shake_update_key(ctx, key, keylen)                                 \
  hash_update((ctx), (key), (keylen))
#define kdf_shake_update_key_uint16_le(ctx, key)                               \
  hash_update_uint16_le((ctx), (key))
#define kdf_shake_finalize_key(ctx) hash_final((ctx))
#define kdf_shake_get_randomness(ctx, dst, count)                              \
  hash_squeeze((ctx), (dst), (count))
#define kdf_shake_clear(ctx) hash_clear((ctx))

// 4x parallel hashing

/* Instances that work with 4 states in parallel. */
typedef Keccak_HashInstancetimes4 hash_context_x4 ATTR_ALIGNED(32);

static inline void hash_init_x4(hash_context_x4 *ctx, size_t digest_size) {
  if (digest_size == 32) {
    Keccak_HashInitializetimes4_SHAKE128(ctx);
  } else {
    Keccak_HashInitializetimes4_SHAKE256(ctx);
  }
}

static inline void hash_update_x4(hash_context_x4 *ctx, const uint8_t **data,
                                  size_t size) {
  Keccak_HashUpdatetimes4(ctx, data, size << 3);
}

static inline void hash_update_x4_4(hash_context_x4 *ctx, const uint8_t *data0,
                                    const uint8_t *data1, const uint8_t *data2,
                                    const uint8_t *data3, size_t size) {
  const uint8_t *data[4] = {data0, data1, data2, data3};
  hash_update_x4(ctx, data, size);
}

static inline void hash_update_x4_1(hash_context_x4 *ctx, const uint8_t *data,
                                    size_t size) {
  const uint8_t *tmp[4] = {data, data, data, data};
  hash_update_x4(ctx, tmp, size);
}

static inline void hash_init_prefix_x4(hash_context_x4 *ctx, size_t digest_size,
                                       const uint8_t prefix) {
  hash_init_x4(ctx, digest_size);
  hash_update_x4_1(ctx, &prefix, sizeof(prefix));
}

static inline void hash_final_x4(hash_context_x4 *ctx) {
  Keccak_HashFinaltimes4(ctx, NULL);
}

static inline void hash_squeeze_x4(hash_context_x4 *ctx, uint8_t **buffer,
                                   size_t buflen) {
  Keccak_HashSqueezetimes4(ctx, buffer, buflen << 3);
}

static inline void hash_squeeze_x4_4(hash_context_x4 *ctx, uint8_t *buffer0,
                                     uint8_t *buffer1, uint8_t *buffer2,
                                     uint8_t *buffer3, size_t buflen) {
  uint8_t *buffer[4] = {buffer0, buffer1, buffer2, buffer3};
  hash_squeeze_x4(ctx, buffer, buflen);
}

#define hash_clear_x4(ctx)
static inline void hash_update_x4_uint16_le(hash_context_x4 *ctx,
                                            uint16_t data) {
  const uint16_t data_le = htole16(data);
  hash_update_x4_1(ctx, (const uint8_t *)&data_le, sizeof(data_le));
}

static inline void hash_update_x4_uint16s_le(hash_context_x4 *ctx,
                                             const uint16_t data[4]) {
  const uint16_t data0_le = htole16(data[0]);
  const uint16_t data1_le = htole16(data[1]);
  const uint16_t data2_le = htole16(data[2]);
  const uint16_t data3_le = htole16(data[3]);
  hash_update_x4_4(ctx, (const uint8_t *)&data0_le, (const uint8_t *)&data1_le,
                   (const uint8_t *)&data2_le, (const uint8_t *)&data3_le,
                   sizeof(data[0]));
}

typedef hash_context_x4 kdf_shake_x4_t;

#define kdf_shake_x4_init(ctx, digest_size) hash_init_x4((ctx), (digest_size))
#define kdf_shake_x4_init_prefix(ctx, digest_size, prefix)                     \
  hash_init_prefix_x4((ctx), (digest_size), (prefix))
#define kdf_shake_x4_update_key(ctx, key, keylen)                              \
  hash_update_x4((ctx), (key), (keylen))
#define kdf_shake_x4_update_key_4(ctx, key0, key1, key2, key3, keylen)         \
  hash_update_x4_4((ctx), (key0), (key1), (key2), (key3), (keylen))
#define kdf_shake_x4_update_key_1(ctx, key, keylen)                            \
  hash_update_x4_1((ctx), (key), (keylen))
#define kdf_shake_x4_update_key_uint16_le(ctx, key)                            \
  hash_update_x4_uint16_le((ctx), (key))
#define kdf_shake_x4_update_key_uint16s_le(ctx, keys)                          \
  hash_update_x4_uint16s_le((ctx), (keys))
#define kdf_shake_x4_finalize_key(ctx) hash_final_x4((ctx))
#define kdf_shake_x4_get_randomness(ctx, dst, count)                           \
  hash_squeeze_x4((ctx), (dst), (count))
#define kdf_shake_x4_get_randomness_4(ctx, dst0, dst1, dst2, dst3, count)      \
  hash_squeeze_x4_4((ctx), (dst0), (dst1), (dst2), (dst3), (count))
#define kdf_shake_x4_clear(ctx) hash_clear_x4((ctx))
#endif