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



#include "ascon_groupe.h"
#include "fichiers_complementaires/tree2.h"


void afficher_hex(std::vector<uint8_t> vec){
    printf("{");
    for(int i=0;i<vec.size();i++){
        printf("%02x",vec[i]);
        if(i!=vec.size()-1){
            printf(", ");
        }
    }
    printf("}\n");
}

void afficher_seeds(SeedTree seed_tree, size_t N, size_t N_ind){
    printf("Les graines des participants sont :\n");
    for(size_t party=0;party<N;party++){
        if(party!=N_ind){
        for(size_t l=0;l<16;l++){
            printf("%02x, ",seed_tree.get_leaf(party).value()[l]);
        }
        printf("\n");
        }
    }
    printf("\n\n");
}

void afficher_salt(banquet_salt_t salt){
    printf("Le sel utilisé sera :\n");
    for(int i=0;i<salt.size();i++){
        printf("%02x",salt[i]);
        if(i!=salt.size()-1){
            printf(", ");
        }
    }
    printf("\n\n");
}



void afficher_type(ascon_type_t x){
    printf("{");
    for(size_t i=0;i<x.size();i++){
        printf("{");
        for(size_t j=0;j<8;j++){
            printf("%02x",x[i][j]);
            if(j!=7){
                printf(", ");
            }
        }
        printf("}");
        if(i!=x.size()-1){
            printf(",\n");
        }
    }
    printf("}\n");
}

void afficher_state(ascon_state_t x){
    printf("{");
    for(size_t i=0;i<5;i++){
        printf("{");
        for(size_t j=0;j<8;j++){
            printf("%02x",x[i][j]);
            if(j!=7){
                printf(", ");
            }
        }
        printf("}");
        if(i!=4){
            printf(",\n");
        }
    }
    printf("}\n");
}


void afficher_key_part(std::vector<ascon_type_t> key_part, size_t N, int rec){
    printf("\n Les parties de la clé données :\n\n");
    for(size_t party=0;party<N;party++){
        afficher_type(key_part[party]);
        printf("\n");
    }
    printf("\n");
    if(rec!=0){
        printf(" la clé recomposée :\n");
        afficher_type(Rec2(key_part));
        printf("\n\n");
    }   
}

void afficher_ciphertext_part(std::vector<ascon_type_t> ciphertext_part, size_t N, int rec){
    printf("\n Les parties de la sortie calculées :\n");
    for(size_t party=0;party<N;party++){
        afficher_type(ciphertext_part[party]);
        printf("\n");
    }
    printf("\n");
    if(rec!=0){
        printf(" la sortie recomposée :\n");
        afficher_type(Rec2(ciphertext_part));
        printf("\n\n");
    }
}

void afficher_in_out_part(std::vector<std::vector<uint8_t>> x1_shares,
                            std::vector<std::vector<uint8_t>> x2_shares,
                            std::vector<std::vector<uint8_t>> y_shares,
                            size_t N){
    printf("\n\n\n \t Les parties x1,x2 et y :\n\n");
    for(int i=0;i<N;i++){
        printf("x1(%d) :\n",i+1);
        for(int j=0;j<720;j++){
            printf("%02x, ",x1_shares[i][j]);
        }
        printf("\n\n");
    }
    for(int i=0;i<N;i++){
        printf("x2(%d) :\n",i+1);
        for(int j=0;j<720;j++){
            printf("%02x, ",x2_shares[i][j]);
        }
        printf("\n\n");
    }
    for(int i=0;i<N;i++){
        printf("y(%d) :\n",i+1);
        for(int j=0;j<720;j++){
            printf("%02x, ",y_shares[i][j]);
        }
        printf("\n\n");
    }
}

void afficher_in_out_rec(std::vector<std::vector<uint8_t>> x1_shares,
                            std::vector<std::vector<uint8_t>> x2_shares,
                            std::vector<std::vector<uint8_t>> y_shares,
                            size_t N){
    std::vector<uint8_t> sum1(720),sum2(720),sumy(720);
    for(int i=0;i<720;i++){
        for(int j=0;j<N;j++){
            sum1[i]^=x1_shares[j][i];
            sum2[i]^=x2_shares[j][i];
            sumy[i]^=y_shares[j][i];
        }
    }
    printf("\n\n x1 :\n");
    for(int i=0;i<720;i++){
        printf("%02x, ",sum1[i]);
    }
    printf("\n\n x2 :\n");
    for(int i=0;i<720;i++){
        printf("%02x, ",sum2[i]);
    }
    printf("\n\n y :\n");
    for(int i=0;i<720;i++){
        printf("%02x, ",sumy[i]);
    }
    printf("\n\n");
}

void afficher_alpha_v_part(std::vector<std::vector<uint8_t>> alpha_el_shares,
                           std::vector<std::vector<uint8_t>> v_e_shares, banquet_instance_t &instance){
    for(size_t party=0;party<instance.num_MPC_parties;party++){
        printf("Participant %d :\nalpha : ",party);
        for(size_t l=0;l<alpha_el_shares[party].size();l++){
            printf("%02x, ",alpha_el_shares[party][l]);
        }
        printf("\nv : ");
        for(size_t l=0;l<v_e_shares[party].size();l++){
            printf("%02x, ",v_e_shares[party][l]);
        }
        printf("\n\n");
    }
}

void afficher_a_c_part(std::vector<std::vector<uint8_t>> a_el_shares,
                       std::vector<std::vector<uint8_t>> c_e_shares, size_t N, size_t r){
    for(size_t party=0;party<N;party++){
        printf("Participant %d :\na : ",party);
        for(size_t l=0;l<r;l++){
            printf("%0x, ",a_el_shares[party][l]);
        }
        printf("\nc : ");
        for(size_t l=0;l<8;l++){
            printf("%0x, ",c_e_shares[party][l]);
        }
        printf("\n\n");
    }                    
}

void afficher_hash(std::vector<uint8_t> h_i, int i){
    printf("\nh_%d :\n",i);
    for(size_t l=0;l<h_i.size();l++){
        printf("%02x ",h_i[l]);
    }
    printf("\n\n");
}

void afficher_R_el(std::vector<uint8_t> R_el, size_t r){
    printf("On étend h_1 pour obtenir les R_el :\n");
    for(size_t l=0;l<r;l++){
        printf("%02x ",R_el[l]);
    }
    printf("\n\n");
}

void afficher_miss_part(uint16_t missing_party){
    printf("\nOn cache la branche %d.\n",missing_party);
}
