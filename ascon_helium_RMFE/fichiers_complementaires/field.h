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

#include "banquet_instances.h"
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <vector>
extern "C" {
#include <smmintrin.h>
#include <wmmintrin.h>
}

namespace field {
class GF2E;
}

field::GF2E dot_product(const std::vector<field::GF2E> &lhs,
                        const std::vector<field::GF2E> &rhs);

namespace field {
class GF2E {

  uint64_t data;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
  static std::function<uint64_t(__m128i)> reduce;
#pragma GCC diagnostic pop
  static size_t byte_size;
  static size_t lambda;
  static uint64_t modulus;

public:
  GF2E() : data(0){};
  GF2E(uint64_t data) : data(data) {}
  GF2E(const GF2E &other) = default;
  ~GF2E() = default;
  GF2E &operator=(const GF2E &other) = default;

  void clear() { data = 0; }
  void set_coeff(size_t idx) { data |= (1ULL << idx); }
  GF2E operator+(const GF2E &other) const;
  GF2E &operator+=(const GF2E &other);
  GF2E operator-(const GF2E &other) const;
  GF2E &operator-=(const GF2E &other);
  GF2E operator*(const GF2E &other) const;
  GF2E &operator*=(const GF2E &other);
  bool operator==(const GF2E &other) const;
  bool operator!=(const GF2E &other) const;

  GF2E inverse() const;

  void to_bytes(uint8_t *out) const;
  std::vector<uint8_t> to_bytes() const;
  void from_bytes(uint8_t *in);
  void afficher() const;
  static void init_extension_field(const banquet_instance_t &instance);
  friend GF2E(::dot_product)(const std::vector<field::GF2E> &lhs,
                             const std::vector<field::GF2E> &rhs);
};

const GF2E &lift_uint8_t(uint8_t value);

std::vector<GF2E> get_first_n_field_elements(size_t n);
std::vector<std::vector<GF2E>>
precompute_lagrange_polynomials(const std::vector<GF2E> &x_values);
std::vector<GF2E> interpolate_with_precomputation(
    const std::vector<std::vector<GF2E>> &precomputed_lagrange_polynomials,
    const std::vector<GF2E> &y_values);

std::vector<GF2E> build_from_roots(const std::vector<GF2E> &roots);
GF2E eval(const std::vector<GF2E> &poly, const GF2E &point);
} // namespace field

std::vector<field::GF2E> operator+(const std::vector<field::GF2E> &lhs,
                                   const std::vector<field::GF2E> &rhs);
std::vector<field::GF2E> &operator+=(std::vector<field::GF2E> &self,
                                     const std::vector<field::GF2E> &rhs);
std::vector<field::GF2E> operator*(const std::vector<field::GF2E> &lhs,
                                   const field::GF2E &rhs);
std::vector<field::GF2E> operator*(const field::GF2E &lhs,
                                   const std::vector<field::GF2E> &rhs);
std::vector<field::GF2E> operator*(const std::vector<field::GF2E> &lhs,
                                   const std::vector<field::GF2E> &rhs);