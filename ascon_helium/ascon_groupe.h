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

std::vector<ascon_state_t> Share(ascon_state_t t, 
            RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, 
            std::vector<std::vector<uint8_t>> *y_shares,
            std::vector<uint8_t> *delta_y){
    std::vector<ascon_state_t> y_shared(N);
    gsl::span<uint8_t> tape;
    ascon_state_t sum=t;
    int size_y=(*y_shares)[0].size();
    for(int i=0;i<N;i++){
        tape=random_tapes.get_bytes(rep,i,*indice,40);
        for(int j=0;j<5;j++){
            for(int l=0;l<8;l++){
                y_shared[i][j][l]=tape[j*8+l];
                ((*y_shares)[i]).push_back(tape[j*8+l]);
                sum[j][l]^=tape[j*8+l];
            }
        }
    }
    for(int j=0;j<5;j++){
        for(int l=0;l<8;l++){
            delta_y->push_back(sum[j][l]);
            y_shared[0][j][l]^=sum[j][l];
            (*y_shares)[0][j*8+l+size_y]^=sum[j][l];
        }
    }
    *indice+=40;
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


void Sbox_groupe_prouveur(std::vector<ascon_state_t> *S,
                 RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, 
                 std::vector<std::vector<uint8_t>> *x1_shares, 
                 std::vector<std::vector<uint8_t>> *x2_shares, 
                 std::vector<std::vector<uint8_t>> *y_shares, 
                 std::vector<uint8_t> *delta_y){
    for(int i=0;i<N;i++){
        XOR_int(&(*S)[i][0],(*S)[i][4]);
        XOR_int(&(*S)[i][2],(*S)[i][1]);
        XOR_int(&(*S)[i][4],(*S)[i][3]);
    }
    uint8_t x1,x2;
    ascon_state_t y;
    for(int j=0;j<5;j++){
        for(int l=0;l<8;l++){
            x1=255; x2=0;
            for(int i=0;i<N;i++){
                x1^=(*S)[i][(j+1)%5][l];
                //Injecter erreure :
                    // if(rep==0 && i==2 && *indice>=50 && *indice<=90){
                    //     x1^=0;
                    // }
                    // else{
                    //     x1^=(*S)[i][(j+1)%5][l];
                    // }
                x2^=(*S)[i][(j+2)%5][l];
                if(i==0){
                    (*x1_shares)[i].push_back((*S)[i][(j+1)%5][l]^255);
                    (*x2_shares)[i].push_back((*S)[i][(j+2)%5][l]);
                }
                else{
                    (*x1_shares)[i].push_back((*S)[i][(j+1)%5][l]);
                    // Injecter erreure :
                        // if(rep==0 && i==2 && *indice>=50 && *indice<=90){
                        //     (*x1_shares)[i].push_back(0x00);
                        // }
                        // else{
                        //     (*x1_shares)[i].push_back((*S)[i][(j+1)%5][l]);
                        // }
                    (*x2_shares)[i].push_back((*S)[i][(j+2)%5][l]);
                }
            }
            y[j][l]=x1&x2;
        }
    }

    std::vector<ascon_state_t> y_shared=Share(y,random_tapes,N,rep,indice,y_shares,delta_y);


    for(int i=0;i<N;i++){
        for(int j=0;j<5;j++){
            XOR_int(&(*S)[i][j],y_shared[i][j]);
        }

        XOR_int(&(*S)[i][1],(*S)[i][0]);
        XOR_int(&(*S)[i][3],(*S)[i][2]);
        XOR_int(&(*S)[i][0],(*S)[i][4]);
    }
    XOR_int(&(*S)[0][2],un);
}




void round_transformation_groupe_prouveur(std::vector<ascon_state_t> *S, int num_round, 
                      RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, 
                      std::vector<std::vector<uint8_t>> *x1_shares,
                      std::vector<std::vector<uint8_t>> *x2_shares,
                      std::vector<std::vector<uint8_t>> *y_shares,
                      std::vector<uint8_t> *delta_y){
    for(int iter=0;iter<num_round;iter++){
        addition_constants(&(*S)[0],num_round,iter);
        Sbox_groupe_prouveur(S,random_tapes,N,rep,indice,x1_shares,x2_shares,y_shares,delta_y);
        for(int i=0;i<N;i++){
            linear_layer(&(*S)[i]);
        }
    }
}



void ascon_128_groupe_prouveur(const std::vector<ascon_type_t> *key_in,
                      const ascon_type_t *plaintext_in, const ascon_type_t *nonce_in,
                      std::vector<ascon_type_t> *ciphertext_out,
                      std::vector<std::vector<uint8_t>> *x1_shares,
                      std::vector<std::vector<uint8_t>> *x2_shares,
                      std::vector<std::vector<uint8_t>> *y_shares,
                      std::vector<uint8_t> *delta_y,
                      RandomTapes random_tapes, size_t rep, size_t *indice){
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

    round_transformation_groupe_prouveur(&S,a_round,random_tapes,N,rep,indice,x1_shares,x2_shares,y_shares,delta_y);
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
            round_transformation_groupe_prouveur(&S,b_round,random_tapes,N,rep,indice,x1_shares,x2_shares,y_shares,delta_y);
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



std::vector<ascon_state_t> Share_verif(RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, 
            std::vector<std::vector<uint8_t>> *y_shares, std::vector<uint8_t> delta_y, size_t *delta_ind){

    std::vector<ascon_state_t> y_shared(N);
    gsl::span<uint8_t> tape;
    int size_y=(*y_shares)[0].size();
    for(int i=0;i<N;i++){
        if(i<N_ind){
            tape=random_tapes.get_bytes(rep,i,*indice,40);
        }
        else{
            tape=random_tapes.get_bytes(rep,i+1,*indice,40);
        }
        for(int j=0;j<5;j++){
            for(int l=0;l<8;l++){
                y_shared[i][j][l]=tape[j*8+l];
                ((*y_shares)[i]).push_back(tape[j*8+l]);
            }
        }
    }

    if(N_ind!=0){
        for(int j=0;j<5;j++){
            for(int l=0;l<8;l++){
                y_shared[0][j][l]^=delta_y[*delta_ind];
                (*y_shares)[0][j*8+l+size_y]^=delta_y[*delta_ind];
                (*delta_ind)++;
            }
        }
    }
    
    *indice+=40;
    return y_shared;
}


void Sbox_groupe_verif(std::vector<ascon_state_t> *S,
                 RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, 
                 std::vector<std::vector<uint8_t>> *x1_shares, 
                 std::vector<std::vector<uint8_t>> *x2_shares, 
                 std::vector<std::vector<uint8_t>> *y_shares, 
                 std::vector<uint8_t> delta_y, size_t *delta_ind){

    for(int i=0;i<N;i++){
        XOR_int(&(*S)[i][0],(*S)[i][4]);
        XOR_int(&(*S)[i][2],(*S)[i][1]);
        XOR_int(&(*S)[i][4],(*S)[i][3]);
    }
    for(int j=0;j<5;j++){
        for(int l=0;l<8;l++){
            for(int i=0;i<N;i++){
                if(i==0 && N_ind!=0){
                    (*x1_shares)[i].push_back((*S)[i][(j+1)%5][l]^255);
                    (*x2_shares)[i].push_back((*S)[i][(j+2)%5][l]);
                }
                else{
                    (*x1_shares)[i].push_back((*S)[i][(j+1)%5][l]);
                    (*x2_shares)[i].push_back((*S)[i][(j+2)%5][l]);
                }
            }
        }
    }

    std::vector<ascon_state_t> y_shared=Share_verif(random_tapes,N,N_ind,rep,indice,y_shares,delta_y,delta_ind);
    

    for(int i=0;i<N;i++){
        for(int j=0;j<5;j++){
            XOR_int(&(*S)[i][j],y_shared[i][j]);
        }

        XOR_int(&(*S)[i][1],(*S)[i][0]);
        XOR_int(&(*S)[i][3],(*S)[i][2]);
        XOR_int(&(*S)[i][0],(*S)[i][4]);
    }
    if(N_ind!=0){
        XOR_int(&(*S)[0][2],un);
    }
}


void round_transformation_groupe_verif(std::vector<ascon_state_t> *S, int num_round, 
                      RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, 
                      std::vector<std::vector<uint8_t>> *x1_shares,
                      std::vector<std::vector<uint8_t>> *x2_shares,
                      std::vector<std::vector<uint8_t>> *y_shares,
                      std::vector<uint8_t> delta_y, size_t *delta_ind){
    for(int iter=0;iter<num_round;iter++){
        if(N_ind!=0){
            addition_constants(&(*S)[0],num_round,iter);
        }
        Sbox_groupe_verif(S,random_tapes,N,N_ind,rep,indice,x1_shares,x2_shares,y_shares,delta_y,delta_ind);
        for(int i=0;i<N;i++){
            linear_layer(&(*S)[i]);
        }
    }
}


void ascon_128_groupe_verif(const std::vector<ascon_type_t> *key_in,
                      const ascon_type_t *plaintext_in, const ascon_type_t *nonce_in,
                      std::vector<ascon_type_t> *ciphertext_out,
                      std::vector<std::vector<uint8_t>> *x1_shares,
                      std::vector<std::vector<uint8_t>> *x2_shares,
                      std::vector<std::vector<uint8_t>> *y_shares,
                      std::vector<uint8_t> delta_y,
                      RandomTapes random_tapes, size_t N_ind, size_t rep, size_t *indice){
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

    round_transformation_groupe_verif(&S,a_round,random_tapes,N,N_ind,rep,indice,x1_shares,x2_shares,y_shares,delta_y,&delta_ind);
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
            round_transformation_groupe_verif(&S,b_round,random_tapes,N,N_ind,rep,indice,x1_shares,x2_shares,y_shares,delta_y,&delta_ind);
        }
    }
}
