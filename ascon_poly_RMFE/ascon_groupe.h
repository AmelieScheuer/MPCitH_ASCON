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
#include "fichiers_complementaires/gen_alea.h"
#include "fichiers_complementaires/rmfe.h"
#include <chrono>

#define TRICHE 0


std::vector<std::vector<uint8_t>> triche_fin={{0x94, 0x3e, 0xc2, 0xa6, 0x61, 0x3b}, {0x3e, 0xfe, 0xac, 0xf5, 0xb6, 0x45},
                                                {0xd5, 0x68, 0x44, 0xd5, 0x06, 0x42}, {0xcb, 0x61, 0x8a, 0xb1, 0x3f, 0xa3}, 
                                                {0x31, 0x30, 0x3d, 0xcf, 0xb2, 0x9b}, {0x08, 0xc5, 0x7c, 0x91, 0xe5, 0x99}, 
                                                {0x22, 0xaa, 0xee, 0xdb, 0x2b, 0x04}, {0x14, 0x6a, 0xda, 0x5c, 0x9b, 0xe9}, 
                                                {0x74, 0x2f, 0xee, 0xe8, 0x6c, 0xe1}, {0xa3, 0x6c, 0xdb, 0xef, 0xdd, 0xa3}, 
                                                {0x2d, 0xd0, 0xe7, 0x10, 0x65, 0xf5}, {0x24, 0xa7, 0xe9, 0x15, 0xcf, 0xf9}, 
                                                {0x26, 0xd9, 0xb5, 0xa3, 0xfa, 0x17}, {0x7d, 0xf7, 0x17, 0xa4, 0xc6, 0xfe}, 
                                                {0x3a, 0x84, 0xd9, 0x6f, 0xe4, 0xdf}, {0x41, 0x88, 0x90, 0xe5, 0x0a, 0x12}, 
                                                {0x1c, 0x33, 0xee, 0xd7, 0x86, 0x15}, {0x40, 0xf8, 0xcd, 0x7d, 0x5b, 0x4c}, 
                                                {0x83, 0x24, 0x3a, 0x2a, 0xbb, 0xbc}, {0x78, 0xc4, 0x34, 0xde, 0x9f, 0x4d}};

std::vector<std::vector<uint8_t>> triche_fin_2={{0x3d, 0xc6, 0xce, 0xfb, 0xd6, 0xe6, }, {0xbd, 0x0d, 0xa7, 0xa8, 0x71, 0x74, }, 
                                                {0x97, 0x22, 0xb9, 0xbe, 0xed, 0x30, }, {0xf6, 0x5c, 0xf0, 0xbe, 0x66, 0xc5, }, 
                                                {0x17, 0xaa, 0x07, 0x3e, 0x17, 0x4a, }, {0x96, 0x35, 0x06, 0xce, 0x8b, 0x4f, }, 
                                                {0xc1, 0x1a, 0x75, 0x3a, 0x4d, 0x84, }, {0x60, 0x6a, 0xfd, 0x5f, 0x76, 0x14, }, 
                                                {0x82, 0x24, 0x7a, 0x34, 0xea, 0x9d, }, {0x50, 0xc8, 0x45, 0x15, 0x33, 0xfe, }, 
                                                {0x85, 0x2e, 0x74, 0x9f, 0x30, 0x70, }, {0xb5, 0xdf, 0x77, 0x3e, 0xe6, 0xbe, }, 
                                                {0x63, 0x0a, 0x08, 0xa0, 0xa0, 0xfe, }, {0x0a, 0x69, 0x41, 0xc3, 0x0a, 0x3a, }, 
                                                {0x1d, 0x6d, 0x7c, 0x6b, 0xda, 0x0c, }, {0x8e, 0x22, 0x60, 0x66, 0x5c, 0x30, }, 
                                                {0x7d, 0xb4, 0x04, 0x49, 0xe3, 0x3a, }, {0xbd, 0x11, 0xfe, 0x6a, 0x70, 0x39, }, 
                                                {0xa3, 0x98, 0xbf, 0x35, 0x6d, 0x9e, }, {0x11, 0x09, 0x04, 0xbc, 0x5e, 0xfe, }};



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


std::vector<std::vector<field::GF2E>> Share_RMFE(std::vector<field::GF2E> y, 
                                        RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, size_t lambda_size,
                                        std::vector<std::vector<field::GF2E>> *Y_shares,
                                        std::vector<field::GF2E> *delta_Y){
    std::vector<std::vector<field::GF2E>> y_shared(N);
    gsl::span<uint8_t> tape;
    std::vector<field::GF2E> sum=y;
    size_t size_vec=y.size();
    size_t size_y=(*Y_shares)[0].size();
    uint64_t tmp;

    for(size_t i=0;i<N;i++){
        tape=random_tapes.get_bytes(rep,i,*indice,size_vec*lambda_size);
        for(size_t j=0;j<size_vec;j++){
            tmp=0;
            for(int z=0;z<lambda_size;z++){
              tmp<<=8;
              tmp^=tape[j*lambda_size+z];
            }
            y_shared[i].push_back(field::GF2E(tmp));
            (*Y_shares)[i].push_back(field::GF2E(tmp));
            sum[j]+=field::GF2E(tmp);
        }
    }

    for(size_t j=0;j<size_vec;j++){
        delta_Y->push_back(sum[j]);
        y_shared[0][j]+=sum[j];
        (*Y_shares)[0][size_y+j]+=sum[j];
    }

    (*indice)+=size_vec*lambda_size;
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




std::vector<std::vector<uint8_t>> state_to_RMFEvec(ascon_state_t input, size_t vecsize){
    std::pair<size_t, size_t> vecpair={vecsize>>3,vecsize&7};
    size_t lambda_size=vecpair.first;
    size_t len_vec;
    
    if(vecpair.second==0){   
        len_vec=(size_t) 40./lambda_size;
        if(len_vec-40./lambda_size!=0){
            len_vec++;
        }
        std::vector<std::vector<uint8_t>> output(len_vec);
        size_t quotient=0,reste=0;
        for(size_t i=0;i<len_vec;i++){
            for(size_t j=0;j<lambda_size;j++){
                if(quotient<5){
                    output[i].push_back(input[quotient][reste]);
                    reste++;
                    if(reste>=8){
                        reste=0;
                        quotient++;
                    }
                }
                else{
                    output[i].push_back(0);
                }
            }
        }
        return output;
    }
    else{
        lambda_size++;
        len_vec=(size_t) 40./lambda_size;
        if(len_vec-40./lambda_size!=0){
            len_vec++;
        }
        throw std::runtime_error("not implemented for this specific RMFE parameters.");
    }
}

ascon_state_t RMFEvec_to_state(std::vector<std::vector<uint8_t>> input, size_t vecsize){
    size_t lambda_size=input[0].size();
    size_t len_vec=input.size();
    ascon_state_t output;

    if((vecsize&7)==0){
        size_t quotient=0,reste=0;
        for(size_t i=0;i<5;i++){
            for(size_t j=0;j<8;j++){
                output[i][j]=input[quotient][reste];
                reste++;
                if(reste==lambda_size){
                    reste=0;
                    quotient++;
                }
            }
        }
        return output;
    }
    else{
        throw std::runtime_error("not implemented for this specific RMFE parameters.");
    }
}


std::vector<field::GF2E> state_to_field(ascon_state_t input, banquet_instance_t &instance){
    size_t vecsize=instance.RMFE_vecsize;
    size_t lambda_size=instance.RMFE_numbytes;
    size_t len_vec;
    
    if((vecsize&7)==0){   
        len_vec=(size_t) 40./lambda_size;
        if(len_vec-40./lambda_size!=0){
            len_vec++;
        }
        std::vector<field::GF2E> output(len_vec);
        std::vector<uint8_t> tmp(lambda_size);
        size_t quotient=0,reste=0;
        for(size_t i=0;i<len_vec;i++){
            for(size_t j=0;j<lambda_size;j++){
                if(quotient<5){
                    tmp[j]=input[quotient][reste];
                    reste++;
                    if(reste>=8){
                        reste=0;
                        quotient++;
                    }
                }
                else{
                    tmp[j]=0;
                }
            }
            output[i]=RMFE::phi(tmp);
        }
        return output;
    }
    else{
        lambda_size++;
        len_vec=(size_t) 40./lambda_size;
        if(len_vec-40./lambda_size!=0){
            len_vec++;
        }
        throw std::runtime_error("not implemented for this specific RMFE parameters.");
    }
}

ascon_state_t field_to_state(std::vector<field::GF2E> input, banquet_instance_t &instance){
    size_t lambda_size=instance.RMFE_numbytes;
    size_t len_vec=input.size();
    ascon_state_t output;

    if((instance.RMFE_vecsize&7)==0){
        size_t quotient=0,reste=0;
        std::vector<uint8_t> tmp=RMFE::psi(input[quotient]);
        for(size_t i=0;i<5;i++){
            for(size_t j=0;j<8;j++){
                output[i][j]=tmp[reste];
                reste++;
                if(reste==lambda_size){
                    reste=0;
                    quotient++;
                    tmp=RMFE::psi(input[quotient]);
                }
            }
        }
        return output;
    }
    else{
        throw std::runtime_error("not implemented for this specific RMFE parameters.");
    }
}



void Sbox_groupe_prouveur(std::vector<ascon_state_t> *S,
                 RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, size_t vecsize, size_t lambda_size,
                 std::vector<std::vector<field::GF2E>> *X1_shares, 
                 std::vector<std::vector<field::GF2E>> *X2_shares, 
                 std::vector<std::vector<field::GF2E>> *Y_shares, 
                 std::vector<field::GF2E> *delta_Y){
    for(int i=0;i<N;i++){
        XOR_int(&(*S)[i][0],(*S)[i][4]);
        XOR_int(&(*S)[i][2],(*S)[i][1]);
        XOR_int(&(*S)[i][4],(*S)[i][3]);
    }
    
    std::vector<field::GF2E> S_mauvais;
    if(TRICHE!=0){
        if(*indice==1336 || *indice==2056){
            ascon_state_t T;
            T=Rec((*S));
            std::vector<std::vector<uint8_t>> tmp_S_mauvais;
            tmp_S_mauvais=state_to_RMFEvec(T,vecsize);
            for(size_t j=0;j<tmp_S_mauvais.size();j++){
                S_mauvais.push_back(RMFE::phi(tmp_S_mauvais[j]));
            }
        }
    }

    std::vector<field::GF2E> x1,x2,y;
    std::vector<std::vector<uint8_t>> tmp_x1,tmp_x2;
    field::GF2E f1,f2;
    for(size_t i=0;i<N;i++){
        if(i==0){
            tmp_x1=state_to_RMFEvec(XOR_state(ShiftRows((*S)[i],1),un_state),vecsize);
            x1=std::vector<field::GF2E>(tmp_x1.size());
            x2=std::vector<field::GF2E>(tmp_x1.size());
        }
        else{
            tmp_x1=state_to_RMFEvec(ShiftRows((*S)[i],1),vecsize);
        }
        tmp_x2=state_to_RMFEvec(ShiftRows((*S)[i],2),vecsize);
        for(size_t j=0;j<tmp_x1.size();j++){
            f1=RMFE::phi(tmp_x1[j]);
            f2=RMFE::phi(tmp_x2[j]);
            (*X1_shares)[i].push_back(f1);
            (*X2_shares)[i].push_back(f2);
            x1[j]+=f1;
            x2[j]+=f2;
        }
        
    }
    for(size_t j=0;j<x1.size();j++){
        if(TRICHE==0){
            y.push_back(x1[j]*x2[j]);
        }
        else{
            if(*indice==2056){
                field::GF2E res_tmp;
                res_tmp.from_bytes(triche_fin_2[j].data());
                res_tmp+=S_mauvais[j];
                y.push_back(res_tmp);
            }
            else{
                if(*indice == 1336){
                    field::GF2E res_tmp;
                    res_tmp.from_bytes(triche_fin[j].data());
                    res_tmp+=S_mauvais[j];
                    y.push_back(res_tmp);
                }
                else{
                    y.push_back(x1[j]*x2[j]);
                }
            }
        }
    }

    std::vector<std::vector<field::GF2E>> y_shares=Share_RMFE(y,random_tapes,N,rep,indice,lambda_size,Y_shares,delta_Y);
    

    std::vector<std::vector<uint8_t>> tmp_y(y_shares[0].size());
    for(int i=0;i<N;i++){
        for(size_t j=0;j<y_shares[0].size();j++){
            tmp_y[j]=RMFE::psi(y_shares[i][j]);
        }
        XOR_state_int(&(*S)[i],RMFEvec_to_state(tmp_y,vecsize));

        XOR_int(&(*S)[i][1],(*S)[i][0]);
        XOR_int(&(*S)[i][3],(*S)[i][2]);
        XOR_int(&(*S)[i][0],(*S)[i][4]);
    }
    XOR_int(&(*S)[0][2],un);
}




void round_transformation_groupe_prouveur(std::vector<ascon_state_t> *S, int num_round, 
                    RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, size_t vecsize, size_t lambda_size,
                    std::vector<std::vector<field::GF2E>> *X1_shares, 
                    std::vector<std::vector<field::GF2E>> *X2_shares, 
                    std::vector<std::vector<field::GF2E>> *Y_shares, 
                    std::vector<field::GF2E> *delta_Y){
    for(int iter=0;iter<num_round;iter++){
        addition_constants(&(*S)[0],num_round,iter);
        Sbox_groupe_prouveur(S,random_tapes,N,rep,indice,vecsize,lambda_size,X1_shares,X2_shares,Y_shares,delta_Y);
        for(int i=0;i<N;i++){
            linear_layer(&(*S)[i]);
        }
    }
}



void ascon_128_groupe_prouveur(const std::vector<ascon_type_t> *key_in,
                        const ascon_type_t *plaintext_in, const ascon_type_t *nonce_in,
                        std::vector<ascon_type_t> *ciphertext_out,
                        std::vector<std::vector<field::GF2E>> *X1_shares, 
                        std::vector<std::vector<field::GF2E>> *X2_shares, 
                        std::vector<std::vector<field::GF2E>> *Y_shares, 
                        std::vector<field::GF2E> *delta_Y,
                        RandomTapes random_tapes, size_t rep, size_t *indice, size_t vecsize, size_t lambda_size){
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

    round_transformation_groupe_prouveur(&S,a_round,random_tapes,N,rep,indice,vecsize,lambda_size,X1_shares,X2_shares,Y_shares,delta_Y);


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
            round_transformation_groupe_prouveur(&S,b_round,random_tapes,N,rep,indice,vecsize,lambda_size,X1_shares,X2_shares,Y_shares,delta_Y);
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


std::vector<std::vector<field::GF2E>> Share_verif_RMFE(RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, 
                                        size_t size_vec, size_t lambda_size, std::vector<std::vector<field::GF2E>> *Y_shares,
                                        std::vector<field::GF2E> delta_Y, size_t *delta_ind){
    std::vector<std::vector<field::GF2E>> y_shared(N);
    gsl::span<uint8_t> tape;
    size_t size_y=(*Y_shares)[0].size();
    uint64_t tmp;

    for(size_t i=0;i<N;i++){
        if(i<N_ind){
            tape=random_tapes.get_bytes(rep,i,*indice,size_vec*lambda_size);
        }
        else{
            tape=random_tapes.get_bytes(rep,i+1,*indice,size_vec*lambda_size);
        }
        for(size_t j=0;j<size_vec;j++){
            tmp=0;
            for(int z=0;z<lambda_size;z++){
              tmp<<=8;
              tmp^=tape[j*lambda_size+z];
            }
            y_shared[i].push_back(field::GF2E(tmp));
            (*Y_shares)[i].push_back(field::GF2E(tmp));
        }
    }

    if(N_ind!=0){
        for(size_t j=0;j<size_vec;j++){
            y_shared[0][j]+=delta_Y[*delta_ind];
            (*Y_shares)[0][size_y+j]+=delta_Y[*delta_ind];
            (*delta_ind)++;
        }
    }

    *indice+=size_vec*lambda_size;
    return y_shared;
}



void Sbox_groupe_verif(std::vector<ascon_state_t> *S,
                 RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, size_t vecsize, size_t lambda_size,
                 std::vector<std::vector<field::GF2E>> *X1_shares, 
                 std::vector<std::vector<field::GF2E>> *X2_shares, 
                 std::vector<std::vector<field::GF2E>> *Y_shares, 
                 std::vector<field::GF2E> delta_Y, size_t *delta_ind){
    for(int i=0;i<N;i++){
        XOR_int(&(*S)[i][0],(*S)[i][4]);
        XOR_int(&(*S)[i][2],(*S)[i][1]);
        XOR_int(&(*S)[i][4],(*S)[i][3]);
    }
    
    std::vector<std::vector<uint8_t>> tmp_x1,tmp_x2;
    field::GF2E f1,f2;
    for(size_t i=0;i<N;i++){
        if(i==0 && N_ind!=0){
            tmp_x1=state_to_RMFEvec(XOR_state(ShiftRows((*S)[i],1),un_state),vecsize);
        }
        else{
            tmp_x1=state_to_RMFEvec(ShiftRows((*S)[i],1),vecsize);
        }
        tmp_x2=state_to_RMFEvec(ShiftRows((*S)[i],2),vecsize);
        for(size_t j=0;j<tmp_x1.size();j++){
            f1=RMFE::phi(tmp_x1[j]);
            f2=RMFE::phi(tmp_x2[j]);
            (*X1_shares)[i].push_back(f1);
            (*X2_shares)[i].push_back(f2);
        }
        
    }
    size_t y_size=tmp_x1.size();

    std::vector<std::vector<field::GF2E>> y_shares=Share_verif_RMFE(random_tapes,N,N_ind,rep,indice,y_size,lambda_size,Y_shares,delta_Y,delta_ind);

    std::vector<std::vector<uint8_t>> tmp_y(y_size);
    for(int i=0;i<N;i++){
        for(size_t j=0;j<y_size;j++){
            tmp_y[j]=RMFE::psi(y_shares[i][j]);
        }
        XOR_state_int(&(*S)[i],RMFEvec_to_state(tmp_y,vecsize));

        XOR_int(&(*S)[i][1],(*S)[i][0]);
        XOR_int(&(*S)[i][3],(*S)[i][2]);
        XOR_int(&(*S)[i][0],(*S)[i][4]);
    }
    if(N_ind!=0){
        XOR_int(&(*S)[0][2],un);
    }
}



void round_transformation_groupe_verif(std::vector<ascon_state_t> *S, int num_round, 
                    RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, size_t vecsize, size_t lambda_size,
                    std::vector<std::vector<field::GF2E>> *X1_shares, 
                    std::vector<std::vector<field::GF2E>> *X2_shares, 
                    std::vector<std::vector<field::GF2E>> *Y_shares, 
                    std::vector<field::GF2E> delta_Y, size_t *delta_ind){
    for(int iter=0;iter<num_round;iter++){
        if(N_ind!=0){
            addition_constants(&(*S)[0],num_round,iter);
        }
        Sbox_groupe_verif(S,random_tapes,N,N_ind,rep,indice,vecsize,lambda_size,X1_shares,X2_shares,Y_shares,delta_Y,delta_ind);
        for(int i=0;i<N;i++){
            linear_layer(&(*S)[i]);
        }
    }
}



void ascon_128_groupe_verif(const std::vector<ascon_type_t> *key_in,
                        const ascon_type_t *plaintext_in, const ascon_type_t *nonce_in,
                        std::vector<ascon_type_t> *ciphertext_out,
                        std::vector<std::vector<field::GF2E>> *X1_shares, 
                        std::vector<std::vector<field::GF2E>> *X2_shares, 
                        std::vector<std::vector<field::GF2E>> *Y_shares, 
                        std::vector<field::GF2E> delta_Y,
                        RandomTapes random_tapes, size_t N_ind, size_t rep, size_t *indice, size_t vecsize, size_t lambda_size){
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

    round_transformation_groupe_verif(&S,a_round,random_tapes,N,N_ind,rep,indice,vecsize,lambda_size,X1_shares,X2_shares,Y_shares,delta_Y,&delta_ind);


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
            round_transformation_groupe_verif(&S,b_round,random_tapes,N,N_ind,rep,indice,vecsize,lambda_size,X1_shares,X2_shares,Y_shares,delta_Y,&delta_ind);
        }
    }
}


