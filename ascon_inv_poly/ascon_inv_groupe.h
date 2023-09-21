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



#include "ascon_inv.h"
#include "fichiers_complementaires/gen_alea.h"
#include <chrono>



std::vector<ascon_type_t> Share2(ascon_type_t t, 
            RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, 
            std::vector<std::vector<uint8_t>> *y_shares,
            std::vector<uint8_t> *delta_y){
    std::vector<ascon_type_t> y_shared(N);
    gsl::span<uint8_t> tape;
    ascon_type_t sum=t;
    int size=t.size();
    int size_y=(*y_shares)[0].size();
    for(int i=0;i<N;i++){
        y_shared[i]=t;
        tape=random_tapes.get_bytes(rep,i,*indice,8*size);
        for(int j=0;j<size;j++){
            for(int l=0;l<8;l++){
                y_shared[i][j][l]=tape[j*8+l];
                ((*y_shares)[i]).push_back(tape[j*8+l]);
                sum[j][l]^=tape[j*8+l];
            }
        }
    }

    

    for(int j=0;j<size;j++){
        for(int l=0;l<8;l++){
            delta_y->push_back(sum[j][l]);
            y_shared[0][j][l]^=sum[j][l];
            (*y_shares)[0][j*8+l+size_y]^=sum[j][l];
        }
    }

    *indice+=8*size;
    return y_shared;
}


std::vector<std::vector<field::GF2E>> Share(std::vector<field::GF2E> y, 
                                        RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, taille_lambda lambda,
                                        std::vector<std::vector<field::GF2E>> *Y_shares,
                                        std::vector<field::GF2E> *delta_Y){
    std::vector<std::vector<field::GF2E>> y_shared(N);
    gsl::span<uint8_t> tape;
    std::vector<field::GF2E> sum=y;
    size_t size_vec=y.size();
    size_t size_y=(*Y_shares)[0].size();
    std::vector<uint8_t> tmp_vec(lambda.byte_size);

    for(size_t i=0;i<N;i++){
        tape=random_tapes.get_bytes(rep,i,*indice,size_vec*lambda.byte_size);
        y_shared[i]=std::vector<field::GF2E>(size_vec);
        for(size_t j=0;j<size_vec;j++){
            for(size_t b=0;b<lambda.byte_size;b++){
                tmp_vec[b]=tape[j*lambda.byte_size+b];
            }
            y_shared[i][j].from_bytes(tmp_vec.data());
            (*Y_shares)[i].push_back(y_shared[i][j]);
            sum[j]+=y_shared[i][j];
        }
    }

    for(size_t j=0;j<size_vec;j++){
        delta_Y->push_back(sum[j]);
        y_shared[0][j]+=sum[j];
        (*Y_shares)[0][size_y+j]+=sum[j];
    }

    (*indice)+=size_vec*lambda.byte_size;
    return y_shared;
}


ascon_type_t Rec2(std::vector<ascon_type_t> x){
    ascon_type_t y(x[0].size());
    for(int i=0;i<x.size();i++){
        for(int j=0;j<y.size();j++){
            for(int l=0;l<8;l++){
                y[j][l]^=x[i][j][l];
            }
        }
    }
    return y;
}

ascon_state_t Rec(std::vector<ascon_state_t> x){
    ascon_state_t y=x[0];
    for(int i=1;i<x.size();i++){
        for(int j=0;j<5;j++){
            for(int l=0;l<8;l++){
                y[j][l]^=x[i][j][l];
            }
        }
    }
    return y;
}



void State_to_Field(std::vector<ascon_state_t> *S, taille_lambda lambda, 
                    std::vector<std::vector<field::GF2E>> *X_shares, std::vector<field::GF2E> *out_rec){
    uint8_t coeffs[8]={0x80,0x40,0x20,0x10,0x8,0x4,0x2,0x1};
    size_t number_elt=lambda.number_elt;
    size_t size_vec_deb=lambda.byte_size;

    field::GF2E valeur(0);
    (*out_rec)=std::vector<field::GF2E>(number_elt);

    for(size_t party=0;party<S->size();party++){
        size_t size_vec=size_vec_deb;
        std::vector<uint8_t> vec(size_vec); 

        size_t ind_coeffs=0;
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
            vec[size_vec-1]=0;
            compteur_vec=size_vec-1;
            compteur_taille=0;
            for(size_t j=0;j<lambda2;j++){
                for(size_t k=0;k<5;k++){
                    vec[compteur_vec]<<=1;
                    vec[compteur_vec]^=((*S)[party][k][compteur_S] & coeffs[ind_coeffs])>>(7-ind_coeffs);

                    compteur_taille++;
                    if(compteur_vec==size_vec-1 && compteur_taille==((5*lambda2)&7)  ||  compteur_vec!=0 && compteur_taille==8){
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
            valeur.from_bytes(vec.data());
            (*X_shares)[party].push_back(valeur);
            (*out_rec)[i]+=valeur;
        }    
    }
}

void Field_to_State(std::vector<ascon_state_t> *S, taille_lambda lambda, std::vector<std::vector<field::GF2E>> *in){
    uint8_t coeffs[8]={0x80,0x40,0x20,0x10,0x8,0x4,0x2,0x1};
    size_t debut=8-((lambda.lambda*5) & 7);
    size_t ouv,compteur;
    size_t ligne_S,colonne_S;

    size_t size_in=(*in)[0].size();

    for(size_t party=0;party<S->size();party++){
        colonne_S=0;
        ligne_S=0;

        (*S)[party][0][0]=0; (*S)[party][1][0]=0; (*S)[party][2][0]=0; (*S)[party][3][0]=0; (*S)[party][4][0]=0; 

        ouv=0; compteur=0;
        for(size_t i=0;i<size_in;i++){
            std::vector<uint8_t> vec=(*in)[party][i].to_bytes();
            if(i==(size_in-1) && 64%lambda.lambda!=0){
                ouv=8-(((64%lambda.lambda)*5) & 7);
                vec.resize(((64%lambda.lambda)*5>>3) + 1);
            }
            else{
                ouv=debut;
            }
            for(int j=vec.size()-1;j>-1;j-=1){
                for(size_t k=ouv;k<8;k++){
                    (*S)[party][ligne_S][colonne_S]<<=1;
                    (*S)[party][ligne_S][colonne_S]^=(vec[j] & coeffs[k])>>(7-k);

                    ligne_S++;
                    if(ligne_S==5){
                        ligne_S=0;
                        compteur++;
                    }
                    if(compteur==8 && colonne_S<7){
                        compteur=0;
                        colonne_S++;
                        (*S)[party][0][colonne_S]=0; (*S)[party][1][colonne_S]=0; (*S)[party][2][colonne_S]=0; (*S)[party][3][colonne_S]=0; (*S)[party][4][colonne_S]=0; 
                    }
                }
                ouv=0;
            }
        }
    }
}


void Sbox_groupe_prouveur(std::vector<ascon_state_t> *S,
                 RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, taille_lambda lambda,
                 std::vector<std::vector<field::GF2E>> *X_shares, 
                 std::vector<std::vector<field::GF2E>> *Y_shares, 
                 std::vector<field::GF2E> *delta_Y){
    std::vector<field::GF2E> tab_field;
    State_to_Field(S,lambda,X_shares,&tab_field);

    for(size_t i=0;i<tab_field.size();i++){
        tab_field[i]=tab_field[i].inverse();
    }

    std::vector<std::vector<field::GF2E>> tab_shared=Share(tab_field,random_tapes,N,rep,indice,lambda,Y_shares,delta_Y);
    Field_to_State(S,lambda,&tab_shared);
}




void round_transformation_groupe_prouveur(std::vector<ascon_state_t> *S, int num_round, 
                    RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, taille_lambda lambda,
                    std::vector<std::vector<field::GF2E>> *X_shares, 
                    std::vector<std::vector<field::GF2E>> *Y_shares, 
                    std::vector<field::GF2E> *delta_Y){
    for(int iter=0;iter<num_round;iter++){
        addition_constants(&(*S)[0],num_round,iter);
        Sbox_groupe_prouveur(S,random_tapes,N,rep,indice,lambda,X_shares,Y_shares,delta_Y);
        for(int i=0;i<N;i++){
            linear_layer(&(*S)[i]);
        }
    }
}



void ascon_128_groupe_prouveur(const std::vector<ascon_type_t> *key_in,
                        const ascon_type_t *plaintext_in, const ascon_type_t *nonce_in,
                        std::vector<ascon_type_t> *ciphertext_out,
                        std::vector<std::vector<field::GF2E>> *X_shares, 
                        std::vector<std::vector<field::GF2E>> *Y_shares, 
                        std::vector<field::GF2E> *delta_Y,
                        RandomTapes random_tapes, size_t rep, size_t *indice, taille_lambda lambda){
    // Initialisation
    size_t N=key_in->size();
    std::vector<ascon_state_t> S(N);
    for(int i=0;i<N;i++){
        if(i==0){
            S[i]={IV[0],(*key_in)[i][0],(*key_in)[i][1],(*nonce_in)[0],(*nonce_in)[1]};
        }
        else{
            S[i][1]=(*key_in)[i][0];
            S[i][2]=(*key_in)[i][1];
        }
        
    }

    round_transformation_groupe_prouveur(&S,a_round,random_tapes,N,rep,indice,lambda,X_shares,Y_shares,delta_Y);


    for(int i=0;i<N;i++){
        XOR_int(&S[i][3],(*key_in)[i][0]);
        XOR_int(&S[i][4],(*key_in)[i][1]);
    }

    for(size_t t=0;t<plaintext_in->size();t++){
        XOR_int(&S[0][0],(*plaintext_in)[t]);
        for(int i=0;i<N;i++){
            (*ciphertext_out)[i].push_back(S[i][0]);
        }
        if(t!=plaintext_in->size()-1){
            round_transformation_groupe_prouveur(&S,b_round,random_tapes,N,rep,indice,lambda,X_shares,Y_shares,delta_Y);
        }
    }
}



















std::vector<ascon_type_t> Share_verif2(size_t size, 
            RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, 
            std::vector<std::vector<uint8_t>> *y_shares,
            std::vector<uint8_t> delta_y, size_t *delta_ind){
    std::vector<ascon_type_t> y_shared(N);
    gsl::span<uint8_t> tape;
    int size_y=(*y_shares)[0].size();
    for(int i=0;i<N;i++){
        y_shared[i]=ascon_type_t(size);
        if(i<N_ind){
            tape=random_tapes.get_bytes(rep,i,*indice,8*size);
        }
        else{
            tape=random_tapes.get_bytes(rep,i+1,*indice,8*size);
        }
        for(int j=0;j<size;j++){
            for(int l=0;l<8;l++){
                y_shared[i][j][l]=tape[j*8+l];
                ((*y_shares)[i]).push_back(tape[j*8+l]);
            }
        }
    }

    if(N_ind!=0){
        for(int j=0;j<size;j++){
            for(int l=0;l<8;l++){
                y_shared[0][j][l]^=delta_y[*delta_ind];
                (*y_shares)[0][j*8+l+size_y]^=delta_y[*delta_ind];
                (*delta_ind)++;
            }
        }
    }
    

    *indice+=8*size;
    return y_shared;
}


std::vector<std::vector<field::GF2E>> Share_verif(RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, 
                                        taille_lambda lambda, std::vector<std::vector<field::GF2E>> *Y_shares,
                                        std::vector<field::GF2E> delta_Y, size_t *delta_ind){
    std::vector<std::vector<field::GF2E>> y_shared(N);
    gsl::span<uint8_t> tape;
    size_t size_y=(*Y_shares)[0].size();
    std::vector<uint8_t> tmp_vec(lambda.byte_size);

    for(size_t i=0;i<N;i++){
        if(i<N_ind){
            tape=random_tapes.get_bytes(rep,i,*indice,lambda.number_elt*lambda.byte_size);
        }
        else{
            tape=random_tapes.get_bytes(rep,i+1,*indice,lambda.number_elt*lambda.byte_size);
        }
        y_shared[i]=std::vector<field::GF2E>(lambda.number_elt);
        for(size_t j=0;j<lambda.number_elt;j++){
            for(int z=0;z<lambda.byte_size;z++){
              tmp_vec[z]=tape[j*lambda.byte_size+z];
            }
            y_shared[i][j].from_bytes(tmp_vec.data());
            (*Y_shares)[i].push_back(y_shared[i][j]);
        }
    }

    if(N_ind!=0){
        for(size_t j=0;j<lambda.number_elt;j++){
            y_shared[0][j]+=delta_Y[*delta_ind];
            (*Y_shares)[0][size_y+j]+=delta_Y[*delta_ind];
            (*delta_ind)++;
        }
    }

    (*indice)+=lambda.number_elt*lambda.byte_size;
    return y_shared;
}



void Sbox_groupe_verif(std::vector<ascon_state_t> *S,
                 RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, taille_lambda lambda,
                 std::vector<std::vector<field::GF2E>> *X_shares, 
                 std::vector<std::vector<field::GF2E>> *Y_shares, 
                 std::vector<field::GF2E> delta_Y, size_t *delta_ind){
    std::vector<field::GF2E> tab_field;
    State_to_Field(S,lambda,X_shares,&tab_field);

    std::vector<std::vector<field::GF2E>> tab_shared=Share_verif(random_tapes,N,N_ind,rep,indice,lambda,Y_shares,delta_Y,delta_ind);
    Field_to_State(S,lambda,&tab_shared);
}



void round_transformation_groupe_verif(std::vector<ascon_state_t> *S, int num_round, 
                    RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, taille_lambda lambda,
                    std::vector<std::vector<field::GF2E>> *X_shares, 
                    std::vector<std::vector<field::GF2E>> *Y_shares, 
                    std::vector<field::GF2E> delta_Y, size_t *delta_ind){
    for(int iter=0;iter<num_round;iter++){
        if(N_ind!=0){
            addition_constants(&(*S)[0],num_round,iter);
        }
        Sbox_groupe_verif(S,random_tapes,N,N_ind,rep,indice,lambda,X_shares,Y_shares,delta_Y,delta_ind);
        for(int i=0;i<N;i++){
            linear_layer(&(*S)[i]);
        }
    }
}



void ascon_128_groupe_verif(const std::vector<ascon_type_t> *key_in,
                        const ascon_type_t *plaintext_in, const ascon_type_t *nonce_in,
                        std::vector<ascon_type_t> *ciphertext_out,
                        std::vector<std::vector<field::GF2E>> *X_shares, 
                        std::vector<std::vector<field::GF2E>> *Y_shares, 
                        std::vector<field::GF2E> delta_Y,
                        RandomTapes random_tapes, size_t N_ind, size_t rep, size_t *indice, taille_lambda lambda){
    // Initialisation
    std::array<uint8_t,8> nul={0,0,0,0,0,0,0,0};
    size_t N=key_in->size();
    size_t delta_ind=0;
    std::vector<ascon_state_t> S(N);
    for(int i=0;i<N;i++){
        if(i==0 && N_ind!=0){
            S[i]={IV[0],(*key_in)[i][0],(*key_in)[i][1],(*nonce_in)[0],(*nonce_in)[1]};
        }
        else{
            S[i][1]=(*key_in)[i][0];
            S[i][2]=(*key_in)[i][1];
        }
        
    }

    round_transformation_groupe_verif(&S,a_round,random_tapes,N,N_ind,rep,indice,lambda,X_shares,Y_shares,delta_Y,&delta_ind);


    for(int i=0;i<N;i++){
        XOR_int(&S[i][3],(*key_in)[i][0]);
        XOR_int(&S[i][4],(*key_in)[i][1]);
    }

    for(size_t t=0;t<plaintext_in->size();t++){
        if(N_ind!=0){
            XOR_int(&S[0][0],(*plaintext_in)[t]);
        }
        for(int i=0;i<N+1;i++){
            if(i<N_ind){
                (*ciphertext_out)[i].push_back(S[i][0]);
            }
            if(i==N_ind){
                (*ciphertext_out)[i].push_back(nul);
            }
            if(i>N_ind){
                (*ciphertext_out)[i].push_back(S[i-1][0]);
            }
            
        }
        if(t!=plaintext_in->size()-1){
            round_transformation_groupe_verif(&S,b_round,random_tapes,N,N_ind,rep,indice,lambda,X_shares,Y_shares,delta_Y,&delta_ind);
        }
    }
}


