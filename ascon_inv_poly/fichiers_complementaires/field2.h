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



#include "field.h"

#include <array>
#include <cstring>
#include <iostream>
#include <stdexcept>



extern "C" {
#include "portable_endian.h"
}

namespace {
std::array<field::GF2E, 256> lifting_lut;

void init_lifting_lut(const field::GF2E &generator) {
  lifting_lut[0] = field::GF2E(0); // lut(0) = 0
  lifting_lut[1] = field::GF2E(1); // lut(1) = 1

  field::GF2E pow = generator;
  for (size_t bit = 1; bit < 8; bit++) {
    size_t start = (1ULL << bit);
    // copy last half of LUT and add current generator power
    for (size_t idx = 0; idx < start; idx++) {
      lifting_lut[start + idx] = lifting_lut[idx] + pow;
    }
    pow = pow * generator;
  }
}

inline __m128i clmul(uint64_t a, uint64_t b) {
  return _mm_clmulepi64_si128(_mm_set_epi64x(0, a), _mm_set_epi64x(0, b), 0);
}



uint64_t GF2_euclidean_div_quotient(uint64_t a, uint64_t b) {
  uint64_t quotient = 0;
  int diff = __builtin_clzl(b) - __builtin_clzl(a);
  while (diff >= 0 && a != 0) {
    quotient |= (1ULL << diff);
    a ^= (b << diff);
    diff = __builtin_clzl(b) - __builtin_clzl(a);
  }
  return quotient;
}

uint64_t mod_inverse(uint64_t a, uint64_t mod) {
  uint64_t t = 0;
  uint64_t new_t = 1;
  uint64_t r = mod;
  uint64_t new_r = a;
  uint64_t tmp;

  while (new_r != 0) {
    uint64_t quotient = GF2_euclidean_div_quotient(r, new_r);
    tmp = r;
    r = new_r;
    new_r = tmp ^ _mm_extract_epi64(clmul(quotient, new_r), 0);
    tmp = t;
    t = new_t;
    new_t = tmp ^ _mm_extract_epi64(clmul(quotient, new_t), 0);
  }

  return t;
}

} // namespace

namespace field {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
std::function<uint64_t(__m128i)> GF2E::reduce = nullptr;
#pragma GCC diagnostic pop
size_t GF2E::lambda = 0;
size_t GF2E::byte_size=0;
uint64_t GF2E::modulus = 0;

GF2E GF2E::operator+(const GF2E &other) const {
  return GF2E(this->data ^ other.data);
}
GF2E &GF2E::operator+=(const GF2E &other) {
  this->data ^= other.data;
  return *this;
}
GF2E GF2E::operator-(const GF2E &other) const {
  return GF2E(this->data ^ other.data);
}
GF2E &GF2E::operator-=(const GF2E &other) {
  this->data ^= other.data;
  return *this;
}
GF2E GF2E::operator*(const GF2E &other) const {
  return GF2E(reduce(clmul(this->data, other.data)));
}
GF2E &GF2E::operator*=(const GF2E &other) {
  this->data = reduce(clmul(this->data, other.data));
  return *this;
}
bool GF2E::operator==(const GF2E &other) const {
  return this->data == other.data;
}
bool GF2E::operator!=(const GF2E &other) const {
  return this->data != other.data;
}

GF2E GF2E::inverse() const { return GF2E(mod_inverse(this->data, modulus)); }

void GF2E::to_bytes(uint8_t *out) const {
  uint64_t be_data = htole64(data);
  memcpy(out, (uint8_t *)(&be_data), byte_size);
}
std::vector<uint8_t> GF2E::to_bytes() const {
  std::vector<uint8_t> buffer(byte_size);
  this->to_bytes(buffer.data());
  return buffer;
}

void GF2E::from_bytes(uint8_t *in) {
  data = 0;
  memcpy((uint8_t *)(&data), in, byte_size);
  data = le64toh(data);
}

void GF2E::afficher() const{
  printf("%0x",data);
}

const GF2E &lift_uint8_t(uint8_t value) { return lifting_lut[value]; }



uint64_t reduce_GF2E_k(__m128i in, size_t lambda, uint64_t modulus){
  // modulus = x^32 + x^7 + x^3 + x^2 + 1
  size_t lambda2=lambda*5;
  uint64_t P = modulus;
  uint64_t mu = P;
  uint64_t R = _mm_cvtsi128_si64(in);
  uint64_t T1 = _mm_cvtsi128_si64(clmul(R >> lambda2, mu));
  uint64_t T2 = _mm_cvtsi128_si64(clmul(T1 >> lambda2, P));

  uint64_t lower_mask=0;
  for(int i=0;i<lambda2;i++){
    lower_mask^=(1<<i);
  }
  return lower_mask & (R ^ T2);
}


void GF2E::init_extension_field(const banquet_instance_t &instance) {
  lambda=instance.lambda.lambda;
  byte_size=instance.lambda.byte_size;
  reduce= [](__m128i in){ return reduce_GF2E_k(in,lambda,modulus); };
  switch (lambda) {
  case 1: {
    // modulus = x^5 + x^2 + 1
    modulus =
        (1ULL << 5) | (1ULL << 2) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 2: {
    // modulus = x^10 + x^3 + 1
    modulus =
        (1ULL << 10) | (1ULL << 3) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 3: {
    // modulus = x^15 + x + 1
    modulus =
        (1ULL << 15) | (1ULL << 1) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 4: {
    // modulus = x^20 + x^3 + 1
    modulus =
        (1ULL << 20) | (1ULL << 3) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 5: {
    // modulus = x^25 + x^3 + 1
    modulus =
        (1ULL << 25) | (1ULL << 3) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 6: {
    // modulus = x^30 + x + 1
    modulus =
        (1ULL << 30) | (1ULL << 1) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 7: {
    // modulus = x^35 + x^2 + 1
    modulus =
        (1ULL << 35) | (1ULL << 2) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 8: {
    // modulus = x^40 + x^5 + x^4 + x^3 + 1
    modulus =
        (1ULL << 40) | (1ULL << 5) | (1ULL << 4) | (1ULL << 3) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 9: {
    // modulus = x^45 + x^4 + x^3 + x + 1
    modulus =
        (1ULL << 45) | (1ULL << 4) | (1ULL << 3) | (1ULL << 1) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 10: {
    // modulus = x^50 + x^4 + x^3 + x^2 + 1
    modulus =
        (1ULL << 50) | (1ULL << 4) | (1ULL << 3) | (1ULL << 2) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 11: {
    // modulus = x^55 + x^6 + x^2 + x + 1
    modulus =
        (1ULL << 55) | (1ULL << 6) | (1ULL << 2) | (1ULL << 1) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  case 12: {
    // modulus = x^60 + x + 1
    modulus =
        (1ULL << 60) | (1ULL << 1) | (1ULL << 0);
    init_lifting_lut(2);
  } break;
  default:
    throw std::runtime_error(
        "modulus for that specific lambda not implemented.");
  }
}


std::vector<GF2E> get_first_n_field_elements(size_t n) {
  std::vector<GF2E> result;
  result.reserve(n);
  for (size_t i = 0; i < n; i++) {
    result.push_back(lifting_lut[i]);
  }
  return result;
}

std::vector<std::vector<GF2E>>
precompute_lagrange_polynomials(const std::vector<GF2E> &x_values) {
  size_t m = x_values.size();
  std::vector<std::vector<GF2E>> precomputed_lagrange_polynomials;
  precomputed_lagrange_polynomials.reserve(m);

  std::vector<GF2E> x_except_k;
  GF2E denominator;
  for (size_t k = 0; k < m; k++) {
    denominator = GF2E(1);
    x_except_k.clear();
    x_except_k.reserve(m - 1);
    for (size_t j = 0; j < m; j++) {
      if (k != j) {
        denominator *= x_values[k] - x_values[j];
        x_except_k.push_back(x_values[j]);
      }
    }
    std::vector<GF2E> numerator = build_from_roots(x_except_k);

    numerator = numerator * denominator.inverse();
    precomputed_lagrange_polynomials.push_back(numerator);
  }

  return precomputed_lagrange_polynomials;
}

std::vector<GF2E> interpolate_with_precomputation(
    const std::vector<std::vector<GF2E>> &precomputed_lagrange_polynomials,
    const std::vector<GF2E> &y_values) {
  if (precomputed_lagrange_polynomials.size() != y_values.size() ||
      y_values.empty())
    throw std::runtime_error("invalid sizes for interpolation");

  std::vector<GF2E> res(precomputed_lagrange_polynomials[0].size());
  size_t m = y_values.size();
  for (size_t k = 0; k < m; k++) {
    res += precomputed_lagrange_polynomials[k] * y_values[k];
  }
  return res;
}

std::vector<GF2E> build_from_roots(const std::vector<GF2E> &roots) {
  size_t len = roots.size();

  std::vector<GF2E> poly(roots);
  poly.push_back(GF2E(0));

  GF2E tmp;
  for (size_t k = 1; k < len; k++) {
    tmp = poly[k];
    poly[k] = tmp + poly[k - 1];
    for (size_t i = k - 1; i >= 1; i--) {
      poly[i] = poly[i] * tmp + poly[i - 1];
    }
    poly[0] *= tmp;
  }
  poly[len] = GF2E(1);
  return poly;
}
// horner eval
GF2E eval(const std::vector<GF2E> &poly, const GF2E &point) {
  GF2E acc;
  long i;

  for (i = poly.size() - 1; i >= 0; i--) {
    acc *= point;
    acc += poly[i];
  }

  return acc;
}

} // namespace field

std::vector<field::GF2E> operator+(const std::vector<field::GF2E> &lhs,
                                   const std::vector<field::GF2E> &rhs) {
  if (lhs.size() != rhs.size())
    throw std::runtime_error("adding vectors of different sizes");

  std::vector<field::GF2E> result(lhs);
  for (size_t i = 0; i < lhs.size(); i++)
    result[i] += rhs[i];

  return result;
}

std::vector<field::GF2E> &operator+=(std::vector<field::GF2E> &lhs,
                                     const std::vector<field::GF2E> &rhs) {
  if (lhs.size() != rhs.size())
    throw std::runtime_error("adding vectors of different sizes");

  for (size_t i = 0; i < lhs.size(); i++)
    lhs[i] += rhs[i];

  return lhs;
}

// somewhat optimized inner product, only do one lazy reduction
field::GF2E dot_product(const std::vector<field::GF2E> &lhs,
                        const std::vector<field::GF2E> &rhs) {

  if (lhs.size() != rhs.size())
    throw std::runtime_error("adding vectors of different sizes");

  // field::GF2E result;
  // for (size_t i = 0; i < lhs.size(); i++)
  // result += lhs[i] * rhs[i];
  __m128i accum = _mm_setzero_si128();
  for (size_t i = 0; i < lhs.size(); i++)
    accum = _mm_xor_si128(accum, clmul(lhs[i].data, rhs[i].data));

  field::GF2E result(field::GF2E::reduce(accum));
  return result;
}

std::vector<field::GF2E> operator*(const std::vector<field::GF2E> &lhs,
                                   const field::GF2E &rhs) {
  std::vector<field::GF2E> result(lhs);
  for (size_t i = 0; i < lhs.size(); i++)
    result[i] *= rhs;

  return result;
}

std::vector<field::GF2E> operator*(const field::GF2E &lhs,
                                   const std::vector<field::GF2E> &rhs) {
  return rhs * lhs;
}

// naive polynomial multiplication
std::vector<field::GF2E> operator*(const std::vector<field::GF2E> &lhs,
                                   const std::vector<field::GF2E> &rhs) {

  std::vector<field::GF2E> result(lhs.size() + rhs.size() - 1);
  for (size_t i = 0; i < lhs.size(); i++)
    for (size_t j = 0; j < rhs.size(); j++)
      result[i + j] += lhs[i] * rhs[j];

  return result;
}

void afficher_poly(std::vector<field::GF2E> poly){
  for(size_t k=0;k<poly.size();k++){
    if(poly[k]!=0){
      poly[k].afficher();
      printf("X^%d + ",k);
    }
  }
}