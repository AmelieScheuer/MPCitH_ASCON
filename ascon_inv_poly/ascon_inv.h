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

#include "fichiers_complementaires/field2.h"
#include <math.h>

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



/* Round transformation */
void addition_constants(ascon_state_t *, int, int);
void Sbox(ascon_state_t *, const taille_lambda);
void linear_layer(ascon_state_t *);

void round_transformation(ascon_state_t *S, int num_round, const taille_lambda lambda){
    for(int indice=0;indice<num_round;indice++){
        addition_constants(S,num_round,indice);
        Sbox(S,lambda);
        linear_layer(S);
    }
}


/* Addition of Constants */
const uint8_t constants[12]={0x4b,0x5a,0x69,0x78,0x87,0x96,0xa5,0xb4,0xc3,0xd2,0xe1,0xf0};

void addition_constants(ascon_state_t *S, int num_round, int indice){
    (*S)[3][7]^=constants[num_round-indice-1];
}


/* Substitution Layer (Sbox) : inversion*/

void state_to_field(ascon_state_t *S, taille_lambda lambda, std::vector<field::GF2E> *out){
    uint8_t coeffs[8]={0x80,0x40,0x20,0x10,0x8,0x4,0x2,0x1};
    size_t ind_coeffs=0;
    size_t number_elt=lambda.number_elt;
    (*out) = std::vector<field::GF2E>(number_elt);

    size_t size_vec=lambda.byte_size;
    std::vector<uint8_t> vec(size_vec); 

    size_t compteur_vec=size_vec-1,compteur_S=0,compteur_taille=0;
    size_t lambda2=lambda.lambda;
    for(size_t i=0;i<number_elt;i++){
        if(i==number_elt-1 && 64%lambda.lambda!=0){
            size_vec=((64%lambda.lambda)*5>>3) + 1;
            lambda2=(64%lambda.lambda);
            for(size_t l=size_vec;l<lambda.byte_size;l++){
                vec[l]=0;
            }
            vec.resize(size_vec);
        }
        (*out)[i]=0;
        vec[size_vec-1]=0;
        compteur_vec=size_vec-1;
        compteur_taille=0;
        for(size_t j=0;j<lambda2;j++){
            for(size_t k=0;k<5;k++){
                vec[compteur_vec]<<=1;
                vec[compteur_vec]^=((*S)[k][compteur_S] & coeffs[ind_coeffs])>>(7-ind_coeffs);

                compteur_taille++;
                if(compteur_vec!=0 && compteur_taille==8 || compteur_vec==size_vec-1 && compteur_taille==((5*lambda2)&7)){
                    compteur_taille=0;
                    compteur_vec--;
                    vec[compteur_vec]=0;
                }
            }
            ind_coeffs++;
            if(ind_coeffs==8){
                ind_coeffs=0;
                compteur_S++;
            }
        }
        (*out)[i].from_bytes(vec.data());
    }    
}

void field_to_state(ascon_state_t *S, taille_lambda lambda, std::vector<field::GF2E> *in){
    uint8_t coeffs[8]={0x80,0x40,0x20,0x10,0x8,0x4,0x2,0x1};

    size_t colonne_S=0;
    size_t ligne_S=0;

    (*S)[0][0]=0; (*S)[1][0]=0; (*S)[2][0]=0; (*S)[3][0]=0; (*S)[4][0]=0; 
    size_t debut=8-((lambda.lambda*5) & 7);

    size_t ouv=0,compteur=0;
    for(size_t i=0;i<in->size();i++){
        std::vector<uint8_t> vec=(*in)[i].to_bytes();
        if(i==in->size()-1 && 64%lambda.lambda!=0){
            ouv=8-(((64%lambda.lambda)*5) & 7);
            vec.resize(((64%lambda.lambda)*5>>3) + 1);
        }
        else{
            ouv=debut;
        }
        for(int j=vec.size()-1;j>-1;j-=1){
            for(size_t k=ouv;k<8;k++){
                (*S)[ligne_S][colonne_S]<<=1;
                (*S)[ligne_S][colonne_S]^=(vec[j] & coeffs[k])>>(7-k);

                ligne_S++;
                if(ligne_S==5){
                    ligne_S=0;
                    compteur++;
                }
                if(compteur==8 && colonne_S<7){
                    compteur=0;
                    colonne_S++;
                    (*S)[0][colonne_S]=0; (*S)[1][colonne_S]=0; (*S)[2][colonne_S]=0; (*S)[3][colonne_S]=0; (*S)[4][colonne_S]=0; 
                }
            }
            ouv=0;
        }
    }
}

void Sbox(ascon_state_t *S, const taille_lambda lambda){
    std::vector<field::GF2E> tab_field;
    state_to_field(S,lambda,&tab_field);

    for(size_t i=0;i<tab_field.size();i++){
        if(tab_field[i]==field::GF2E(0)){
            printf("\nErreur de génération de clé.\n");
        }
        tab_field[i]=tab_field[i].inverse();
    }

    field_to_state(S,lambda,&tab_field);
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
ascon_type_t ascon_128(ascon_type_t K, ascon_type_t N, ascon_type_t P, const taille_lambda lambda){
    // Initialisation
    ascon_state_t S={IV[0],K[0],K[1],N[0],N[1]};
    round_transformation(&S,a_round,lambda);
    XOR_int(&S[3],K[0]);
    XOR_int(&S[4],K[1]);

    // Chiffrement de P
    ascon_type_t C;
    for(size_t t=0;t<P.size();t++){
        XOR_int(&S[0],P[t]);
        C.push_back(S[0]);
        if(t!=P.size()-1){
            round_transformation(&S,b_round,lambda);
        }
    }

    return C;
}


