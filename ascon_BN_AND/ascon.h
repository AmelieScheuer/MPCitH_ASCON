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



#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <array> 
#include <vector>
#include <typeinfo>

/* Introduction des paramètres utilisés dans Ascon-128 */
#define key_bits 128
#define n_bits 128
#define r_bits 64
#define a_round 12
#define b_round 6


/* types spécifiques au chiffrement */
typedef std::array<std::array<uint8_t,8>,5> ascon_state_t;
typedef std::vector<std::array<uint8_t,8>> ascon_type_t;

#define IV ascon_type_t({{key_bits,r_bits,a_round,b_round,0,0,0,0}})





/* Opérations XOR et AND sur les states */ 
std::array<uint8_t,8> XOR(std::array<uint8_t,8> x1, std::array<uint8_t,8> x2){
    std::array<uint8_t,8> res;
    for(int i=0;i<8;i++){
        res[i]=(x1[i]^x2[i]);
    }
    return res;
}

void XOR_int(std::array<uint8_t,8> *x, std::array<uint8_t,8> x2){
    for(int i=0;i<8;i++){
        (*x)[i]^=x2[i];
    }
}

std::array<uint8_t,8> AND(std::array<uint8_t,8> x1, std::array<uint8_t,8> x2){
    std::array<uint8_t,8> res;
    for(int i=0;i<8;i++){
        res[i]=(x1[i]&x2[i]);
    }
    return res;
}


/* Round transformation */
void addition_constants(ascon_state_t *, int, int);
void Sbox(ascon_state_t *);
void linear_layer(ascon_state_t *);

void round_transformation(ascon_state_t *S, int num_round){
    for(int indice=0;indice<num_round;indice++){
        addition_constants(S,num_round,indice);
        Sbox(S);
        linear_layer(S);
    }
}


/* Addition of Constants */
const uint8_t constants[12]={0x4b,0x5a,0x69,0x78,0x87,0x96,0xa5,0xb4,0xc3,0xd2,0xe1,0xf0};

void addition_constants(ascon_state_t *S, int num_round, int indice){
    (*S)[3][7]^=constants[num_round-indice-1];
}


/* Substitution Layer (Sbox) */
const std::array<uint8_t,8> un={255,255,255,255,255,255,255,255};

void Sbox(ascon_state_t *S){
    XOR_int(&(*S)[0],(*S)[4]);
    XOR_int(&(*S)[2],(*S)[1]);
    XOR_int(&(*S)[4],(*S)[3]);

    std::array<uint8_t,8> tmp1=(*S)[0],tmp2=(*S)[1];
    XOR_int(&(*S)[0],AND((*S)[2],XOR(un,(*S)[1])));
    XOR_int(&(*S)[1],AND((*S)[3],XOR(un,(*S)[2])));
    XOR_int(&(*S)[2],AND((*S)[4],XOR(un,(*S)[3])));
    XOR_int(&(*S)[3],AND(tmp1,XOR(un,(*S)[4])));
    XOR_int(&(*S)[4],AND(tmp2,XOR(un,tmp1)));

    XOR_int(&(*S)[1],(*S)[0]);
    XOR_int(&(*S)[3],(*S)[2]);
    XOR_int(&(*S)[0],(*S)[4]);

    XOR_int(&(*S)[2],un);
}


/* Linear Diffusion Layer */
std::array<uint8_t,8> right_shift(std::array<uint8_t,8> x, int i){
    std::array<uint8_t,8> y;
    for(int j=0;j<8;j++){
        y[(j+(i>>3))&7]=x[j];
    }
    int i2=i&7;
    if(i2!=0){
        uint8_t decalage=0;
        for(int j=0;j<i2;j++){
            decalage<<=1;
            decalage^=1;
        }
        uint8_t tmp1=y[7]&decalage;
        uint8_t tmp2;
        for(int j=0;j<8;j++){
            tmp2=(tmp1<<(8-i2))^(y[j]>>i2);
            tmp1=y[j]&decalage;
            y[j]=tmp2;
        }
    }
    return y;
}

void linear_layer(ascon_state_t *S){
    XOR_int(&(*S)[0],XOR(right_shift((*S)[0],19),right_shift((*S)[0],28)));
    XOR_int(&(*S)[1],XOR(right_shift((*S)[1],61),right_shift((*S)[1],39)));
    XOR_int(&(*S)[2],XOR(right_shift((*S)[2],1),right_shift((*S)[2],6)));
    XOR_int(&(*S)[3],XOR(right_shift((*S)[3],10),right_shift((*S)[3],17)));
    XOR_int(&(*S)[4],XOR(right_shift((*S)[4],7),right_shift((*S)[4],41)));
}





/* Algorithme de chiffrement Ascon */
ascon_type_t ascon_128(ascon_type_t K, ascon_type_t N, ascon_type_t P){

    // Initialisation
    ascon_state_t S={IV[0],K[0],K[1],N[0],N[1]};
    round_transformation(&S,a_round);
    XOR_int(&S[3],K[0]);
    XOR_int(&S[4],K[1]);

    // Chiffrement de P
    ascon_type_t C;
    for(size_t t=0;t<P.size();t++){
        XOR_int(&S[0],P[t]);
        C.push_back(S[0]);
        if(t!=P.size()-1){
            round_transformation(&S,b_round);
        }
    }

    return C;
}
