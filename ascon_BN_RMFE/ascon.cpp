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



#include "ascon.h"


int main(){
    printf("IV :\n");
    afficher_type(IV);
    printf("\nclé secrète K :\n");
    ascon_type_t K={{0x0f,0x15,0x71,0xc9,
        0x47,0xd9,0xe8,0x59},
        {0x0c,0xb7,0xad,0xd6,
        0xaf,0x7f,0x67,0x98}};
    afficher_type(K);
    printf("\ncomplément N :\n");
    ascon_type_t N={{0xa3, 0x61, 0xcc, 0xdd,
        0x42, 0x54, 0xf5, 0x53},
        {0x1d, 0xbc, 0xcd, 0x0e, 
        0x9a, 0xf3, 0x52, 0x3c}};
    afficher_type(N);
    printf("\nLe message clair P :\n");
    ascon_type_t P={{0x02, 0x0c, 0xbd, 0xc8,
        0x4e, 0x1d, 0x83, 0x2d},
        {0x4e, 0x65, 0xd5, 0xc7,
        0x15, 0x1e, 0x57, 0x81}};
    afficher_type(P);
    printf("\n\n");
    
    ascon_type_t C=ascon_128(K,N,P);
    printf("ASCON-128(K,N,P) donne comme résultat :\n");
    afficher_type(C);
}