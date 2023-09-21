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



#ifndef BANQUET_INSTANCES_H
#define BANQUET_INSTANCES_H

#include <cstdint>
#include <cstdlib>
#include <stdexcept>

/** Parameter set names */
enum banquet_params_t {
  PARAMETER_SET_INVALID = 0,
  Banquet_L1_Param1 = 1,
  Banquet_L1_Param2 = 2,
  Banquet_L1_Param3 = 3,
  Banquet_L1_Param4 = 4,
  Banquet_L1_Param5 = 5,
  Banquet_L1_Param6 = 6,
  Banquet_L1_Param7 = 7,
  Banquet_L1_Param8 = 8,
  Banquet_L1_Param9 = 9,
  Banquet_L1_Param10 = 10,
  Banquet_L3_Param1 = 11,
  Banquet_L3_Param2 = 12,
  Banquet_L3_Param3 = 13,
  Banquet_L3_Param4 = 14,
  Banquet_L3_Param5 = 15,
  Banquet_L3_Param6 = 16,
  Banquet_L3_Param7 = 17,
  Banquet_L5_Param1 = 18,
  Banquet_L5_Param2 = 19,
  Banquet_L5_Param3 = 20,
  Banquet_L5_Param4 = 21,
  Banquet_L5_Param5 = 22,
  Banquet_L5_Param6 = 23,
  PARAMETER_SET_MAX_INDEX = 24
};

struct banquet_instance_t {
  uint32_t digest_size;     /* bytes */
  uint32_t seed_size;       /* bytes */
  uint32_t num_rounds;      // T
  uint32_t num_MPC_parties; // N

  
  uint32_t d;     // d: number of check
  size_t deg_poly;
  size_t num_poly;
  size_t RMFE_fieldsize; // field size in RMFE
  size_t RMFE_vecsize; // size of vector in RMFE
  size_t RMFE_numbytes; // RMFE_vecsize/8

  banquet_params_t params;
};

#endif