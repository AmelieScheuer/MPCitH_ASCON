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



uint64_t reduce_GF2_16(__m128i in) {
  // modulus = x^16 + x^5 + x^3 + x + 1
  constexpr uint64_t lower_mask = 0xFFFFULL;
  uint64_t R_lower = _mm_cvtsi128_si64(in);
  uint64_t R_upper = R_lower >> 16;

  uint64_t T = R_upper;
  R_upper = R_upper ^ (T >> 11) ^ (T >> 13) ^ (T >> 15);
  R_lower = R_lower ^ (R_upper << 5) ^ (R_upper << 3) ^ (R_upper << 1) ^
            (R_upper << 0);
  return lower_mask & R_lower;
}

uint64_t reduce_GF2_24(__m128i in){
    // modulus = x^24 + x^16 + x^15 + x^14 + x^13 + x^10 + x^9 + x^7 + x^5 + x^3 + 1
    constexpr uint64_t lower_mask = 0xFFFFFFULL;
    uint64_t R_lower = _mm_cvtsi128_si64(in);
    uint64_t R_upper = R_lower >> 24;

    uint64_t T = R_upper;
    R_upper = R_upper ^ (T >>16) ^ (T >> 18) ^ (T >> 20) ^ (T >> 22) ^ (T >> 8) ^ (T >> 9) ^ (T >> 10) ^ (T >> 11) 
              ^ (T >> 14) ^ (T >> 15) ^ (T >> 17) ^ (T >> 19) ^ (T >> 21);
    R_lower = R_lower ^ (R_upper << 16) ^ (R_upper << 15) ^ (R_upper << 14) ^ (R_upper << 13) ^ (R_upper << 10) ^ (R_upper << 9) ^
              (R_upper << 7) ^ (R_upper << 5) ^ (R_upper << 3) ^ (R_upper << 0);
    return lower_mask & R_lower;
}

uint64_t reduce_GF2_32(__m128i in) {
  // modulus = x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1
  constexpr uint64_t lower_mask = 0xFFFFFFFFULL;
  uint64_t R_lower = _mm_cvtsi128_si64(in);
  uint64_t R_upper = R_lower >> 32;

  uint64_t T = R_upper;
  R_upper = R_upper ^ (T >> 17) ^ (T >> 23) ^ (T >> 25) ^ (T >> 28) ^ (T >> 29);
  R_lower = R_lower ^ (R_upper << 15) ^ (R_upper << 9) ^ (R_upper << 7) ^ (R_upper << 4) ^ (R_upper << 3) ^
            (R_upper << 0);
  return lower_mask & R_lower;
}

uint64_t reduce_GF2_40(__m128i in) {
  // modulus = x^40 + x^5 + x^4 + x^3 + 1
  constexpr uint64_t upper_mask = 0xFFFFULL;
  constexpr uint64_t lower_mask = 0xFFFFFFFFFFULL;
  uint64_t R_lower = _mm_extract_epi64(in, 0);
  uint64_t R_upper =
      ((_mm_extract_epi64(in, 1) & upper_mask) << 24) | (R_lower >> 40);

  uint64_t T = R_upper;
  R_upper = R_upper ^ (T >> 35) ^ (T >> 36) ^ (T >> 37);
  R_lower = R_lower ^ (R_upper << 5) ^ (R_upper << 4) ^ (R_upper << 3) ^
            (R_upper << 0);
  return lower_mask & R_lower;
}

uint64_t reduce_GF2_45(__m128i in) {
  // modulus = x^45 + x^20 + x^17 + x^15 + x^14 + x^12 + x^11 + x^6 + 1
  constexpr uint64_t upper_mask = 0x3FFFFFFULL; // mask for getting the relevant
                                                // bits of the upper 64 bit word
  constexpr uint64_t lower_mask =
      0x1FFFFFFFFFFFULL; // mask for getting the relevant bits of the lower 64
                         // bit word
  uint64_t R_lower = _mm_extract_epi64(in, 0);
  uint64_t R_upper = ((_mm_extract_epi64(in, 1) & upper_mask) << 19) |
                     (R_lower >> 45); // build the part to reduce from the upper
                                      // word + remaining from lower word
  // reduce
  uint64_t T = R_upper;
  R_upper = R_upper ^ (T >> 25) ^ (T >> 28) ^ (T >> 30) ^ (T >> 31) ^ (T >> 33) ^ (T >> 34) ^ (T >> 39);
  //R_upper = R_upper ^ (T >> 26);
  R_lower = R_lower ^ (R_upper << 20) ^ (R_upper << 17) ^ (R_upper << 15) ^ (R_upper << 14) ^ (R_upper << 12) ^
            (R_upper << 11) ^ (R_upper << 6) ^ (R_upper << 0);
  return lower_mask & R_lower;
}

uint64_t reduce_GF2_48(__m128i in) {
  // modulus = x^48 + x^25 + x^23 + x^17 + x^12 + x^11 + x^10 + x^8 + x^7 + x^3 + 1
  constexpr uint64_t upper_mask = 0xFFFFFFFFULL;
  constexpr uint64_t lower_mask = 0xFFFFFFFFFFFFULL;
  uint64_t R_lower = _mm_extract_epi64(in, 0);
  uint64_t R_upper =
      ((_mm_extract_epi64(in, 1) & upper_mask) << 16) | (R_lower >> 48);
  uint64_t T = R_upper;
  R_upper = R_upper ^ (T >> 46) ^ (T >> 23) ^ (T >> 25) ^ (T >> 31) ^ (T >> 36) ^ (T >> 37) ^ (T >> 38) ^ (T >> 40) ^ (T >> 41) ^ (T >> 45);
  R_lower = R_lower ^ (R_upper << 25) ^ (R_upper << 23) ^ (R_upper << 17) ^ (R_upper << 12) ^ (R_upper << 11) ^ (R_upper << 10) ^
            (R_upper << 8) ^ (R_upper << 7) ^ (R_upper << 3) ^ (R_upper << 0);
  return lower_mask & R_lower;
}

uint64_t reduce_GF2_64(__m128i in) {
  // modulus = x^64 + x^33 + x^30 + x^26 + x^25 + x^24 + x^23 + x^22 + x^21 + x^20 + x^18 + x^13 + x^12 + x^11 + x^10 + x^7 + x^5 + x^4 + x^2 + x + 1
  
  uint64_t R_lower = _mm_extract_epi64(in, 0);
  uint64_t R_upper = _mm_extract_epi64(in, 1);
  
  uint64_t T = R_upper;
  R_upper = R_upper ^ (T >> 62) ^ (T >> 31) ^ (T >> 34) ^ (T >> 38) ^ (T >> 39) ^ (T >> 40) ^ (T >> 41) ^ (T >> 42)
            ^ (T >> 43) ^ (T >> 44) ^ (T >> 46) ^ (T >> 51) ^ (T >> 52) ^ (T >> 53) ^ (T >> 54) ^ (T >> 57)
            ^ (T >> 59) ^ (T >> 60) ^ (T >> 62) ^ (T >> 63);
  //R_upper = R_upper ^ (T >> 26);
  R_lower = R_lower ^ (R_upper << 33) ^ (R_upper << 30) ^ (R_upper << 26) ^ (R_upper << 25) ^ (R_upper << 24) ^
            (R_upper << 23) ^ (R_upper << 22) ^ (R_upper << 21) ^ (R_upper << 20) ^ (R_upper << 18) ^ (R_upper << 13) ^
            (R_upper << 12) ^ (R_upper << 11) ^ (R_upper << 10) ^(R_upper << 7) ^ (R_upper << 5) ^ (R_upper << 4)
            ^ (R_upper << 2) ^ (R_upper << 1) ^ (R_upper << 0);
  return R_lower;
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
  uint64_t P = modulus;
  uint64_t mu = P;
  uint64_t R = _mm_cvtsi128_si64(in);
  uint64_t T1 = _mm_cvtsi128_si64(clmul(R >> lambda, mu));
  uint64_t T2 = _mm_cvtsi128_si64(clmul(T1 >> lambda, P));

  uint64_t lower_mask=0;
  for(int i=0;i<lambda;i++){
    lower_mask|=(1ULL<<i);
  }
  return lower_mask & (R ^ T2);
}


void GF2E::init_extension_field(const banquet_instance_t &instance) {
  lambda=instance.RMFE_fieldsize;
  if((lambda & 7) == 0){
    byte_size=lambda>>3;
  }
  else{
    byte_size=(lambda>>3)+1;
  }
  reduce= [](__m128i in){ return reduce_GF2E_k(in,lambda,modulus); };
  init_lifting_lut(2);
  switch (lambda) {
  case 2: {
    // modulus = x^2 + x + 1
    modulus =
        (1ULL << 2) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 3: {
    // modulus = x^3 + x + 1
    modulus =
        (1ULL << 3) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 4: {
    // modulus = x^4 + x + 1
    modulus =
        (1ULL << 4) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 5: {
    // modulus = x^5 + x^2 + 1
    modulus =
        (1ULL << 5) | (1ULL << 2) | (1ULL << 0);
  } break;
  case 6: {
    // modulus = x^6 + x + 1
    modulus =
        (1ULL << 6) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 7: {
    // modulus = x^7 + x + 1
    modulus =
        (1ULL << 7) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 8: {
    // modulus = x^8 + x^4 + x^3 + x + 1
    modulus =
        (1ULL << 8) | (1ULL << 4) | (1ULL << 3) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 9: {
    // modulus = x^9 + x^4 + 1
    modulus =
        (1ULL << 9) | (1ULL << 4) | (1ULL << 0);
  } break;
  case 10: {
    // modulus = x^10 + x^3 + 1
    modulus =
        (1ULL << 10) | (1ULL << 3) | (1ULL << 0);
  } break;
  case 11: {
    // modulus = x^11 + x^2 + 1
    modulus =
        (1ULL << 11) | (1ULL << 2) | (1ULL << 0);
  } break;
  case 12: {
    // modulus = x^12 + x^3 + 1
    modulus =
        (1ULL << 12) | (1ULL << 3) | (1ULL << 0);
  } break;
  case 13: {
    // modulus = x^13 + x^4 + x^3 + x + 1
    modulus =
        (1ULL << 13) | (1ULL << 4) | (1ULL << 3) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 14: {
    // modulus = x^14 + x^5 + 1
    modulus =
        (1ULL << 14) | (1ULL << 5) | (1ULL << 0);
  } break;
  case 15: {
    // modulus = x^15 + x + 1
    modulus =
        (1ULL << 15) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 16: {
    // modulus = x^16 + x^5 + x^3 + x + 1
    modulus =
        (1ULL << 16) | (1ULL << 5) | (1ULL << 3) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 17: {
    // modulus = x^17 + x^3 + 1
    modulus =
        (1ULL << 17) | (1ULL << 3) | (1ULL << 0);
  } break;
  case 18: {
    // modulus = x^18 + x^3 + 1
    modulus =
        (1ULL << 18) | (1ULL << 3) | (1ULL << 0);
  } break;
  case 19: {
    // modulus = x^19 + x^5 + x^2 + x + 1
    modulus =
        (1ULL << 19) | (1ULL << 5) | (1ULL << 2) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 20: {
    // modulus = x^20 + x^3 + 1
    modulus =
        (1ULL << 20) | (1ULL << 3) | (1ULL << 0);
  } break;
  case 21: {
    // modulus = x^21 + x^6 + x^5 + x^2 + 1
    modulus =
        (1ULL << 21) | (1ULL << 6) | (1ULL << 5) | (1ULL << 2) | (1ULL << 0);
  } break;
  case 22: {
    // modulus = x^22 + x + 1
    modulus =
        (1ULL << 22) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 23: {
    // modulus = x^23 + x^5 + 1
    modulus =
        (1ULL << 23) | (1ULL << 5) | (1ULL << 0);
  } break;
  case 24: {
    // modulus = x^24 + x^16 + x^15 + x^14 + x^13 + x^10 + x^9 + x^7 + x^5 + x^3 + 1
    modulus =
        (1ULL << 24) | (1ULL << 16) | (1ULL << 15) | (1ULL << 14) | (1ULL << 13) | (1ULL << 10) | (1ULL << 9) | (1ULL << 7) | (1ULL << 5) | (1ULL << 3) | (1ULL << 0);
    reduce=reduce_GF2_24;
  } break;
  case 25: {
    // modulus = x^25 + x^3 + 1
    modulus =
        (1ULL << 25) | (1ULL << 3) | (1ULL << 0);
  } break;
  case 26: {
    // modulus = x^26 + x^4 + x^3 + x + 1
    modulus =
        (1ULL << 26) | (1ULL << 4) | (1ULL << 3) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 27: {
    // modulus = x^27 + x^5 + x^2 + x + 1
    modulus =
        (1ULL << 26) | (1ULL << 5) | (1ULL << 2) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 28: {
    // modulus = x^28 + x + 1
    modulus =
        (1ULL << 28) | (1ULL << 1) | (1ULL << 0);
  } break;
  case 32: {
    // modulus = x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1
    modulus =
        (1ULL << 32) | (1ULL << 15) | (1ULL << 9) | (1ULL << 7) | (1ULL << 4) | (1ULL << 3) | (1ULL << 0);
    reduce=reduce_GF2_32;
  } break;
  case 40: {
    // modulus = x^40 + x^5 + x^4 + x^3 + 1
    modulus =
        (1ULL << 40) | (1ULL << 5) | (1ULL << 4) | (1ULL << 3) | (1ULL << 0);
  } break;
  case 45: {
    // modulus = x^45 + x^20 + x^17 + x^15 + x^14 + x^12 + x^11 + x^6 + 1
    modulus =
        (1ULL << 45) | (1ULL << 20) | (1ULL << 17) | (1ULL << 15) | (1ULL << 14) | (1ULL << 12) | (1ULL << 11) | (1ULL << 6) | (1ULL << 0);
    reduce=reduce_GF2_45;
  } break;
  case 48: {
    // modulus = x^48 + x^25 + x^23 + x^17 + x^12 + x^11 + x^10 + x^8 + x^7 + x^3 + 1
    modulus =
        (1ULL << 48) | (1ULL << 25) | (1ULL << 23) | (1ULL << 17) |  (1ULL << 12) | (1ULL << 11) | (1ULL << 10) | (1ULL << 8) | (1ULL << 7) | (1ULL << 3) |(1ULL << 0);
    reduce = reduce_GF2_48;
  } break;
  case 64: {
    // modulus = x^64 + x^33 + x^30 + x^26 + x^25 + x^24 + x^23 + x^22 + x^21 + x^20 + x^18 + x^13 + x^12 + x^11 + x^10 + x^7 + x^5 + x^4 + x^2 + x + 1
    modulus =
        (1ULL << 33) | (1ULL << 30) | (1ULL << 26) | (1ULL << 25) | (1ULL << 24) | (1ULL << 23) | 
        (1ULL << 22) | (1ULL << 21) | (1ULL << 20) | (1ULL << 18) | (1ULL << 13) | (1ULL << 12) | (1ULL << 11) | 
        (1ULL << 10) | (1ULL << 7) | (1ULL << 5) | (1ULL << 4) | (1ULL << 2) | (1ULL << 1) | (1ULL << 0);
    reduce = reduce_GF2_64;
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
  std::vector<uint8_t> vec;
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