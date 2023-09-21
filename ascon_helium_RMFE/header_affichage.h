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

void afficher_in_out_part(std::vector<std::vector<field::GF2E>> X1_shares,
                            std::vector<std::vector<field::GF2E>> X2_shares,
                            std::vector<std::vector<field::GF2E>> Y_shares,
                            size_t N, size_t num_poly){
    printf("\n\n\n \t Les parties x1,x2 et y :\n\n");
    std::vector<uint8_t> vec;
    for(int i=0;i<N;i++){
        printf("x1(%d) :\n",i+1);
        for(int j=0;j<num_poly;j++){
            vec=X1_shares[i][j].to_bytes();
            for(size_t m=0;m<vec.size();m++){
                printf("%02x",vec[vec.size()-m-1]);
            }
            printf(", ");
        }
        printf("\n\n");
    }
    for(int i=0;i<N;i++){
        printf("x2(%d) :\n",i+1);
        for(int j=0;j<num_poly;j++){
            vec=X2_shares[i][j].to_bytes();
            for(size_t m=0;m<vec.size();m++){
                printf("%02x",vec[vec.size()-m-1]);
            }
            printf(", ");
        }
        printf("\n\n");
    }
    for(int i=0;i<N;i++){
        printf("y(%d) :\n",i+1);
        for(int j=0;j<num_poly;j++){
            vec=Y_shares[i][j].to_bytes();
            for(size_t m=0;m<vec.size();m++){
                printf("%02x",vec[vec.size()-m-1]);
            }
            printf(", ");
        }
        printf("\n\n");
    }
}

void afficher_in_out_rec(std::vector<std::vector<field::GF2E>> X1_shares,
                            std::vector<std::vector<field::GF2E>> X2_shares,
                            std::vector<std::vector<field::GF2E>> Y_shares,
                            banquet_instance_t &instance){
    std::vector<field::GF2E> sum1(instance.d),sum2(instance.d),sumy(instance.d);
        for(int i=0;i<instance.d;i++){
            sum1[i]=field::GF2E(0);
            sum2[i]=field::GF2E(0);
            sumy[i]=field::GF2E(0);
            for(int j=0;j<instance.num_MPC_parties;j++){
                sum1[i]+=X1_shares[j][i];
                sum2[i]+=X2_shares[j][i];
                sumy[i]+=Y_shares[j][i];
            }
        }
        std::vector<uint8_t> vec;
        printf("\n\n x1 :\n");
        for(int i=0;i<instance.d;i++){
            vec=sum1[i].to_bytes();
            for(size_t m=0;m<vec.size();m++){
                printf("%02x",vec[vec.size()-m-1]);
            }
            printf(", ");
        }
        printf("\n\n x2 :\n");
        for(int i=0;i<instance.d;i++){
            vec=sum2[i].to_bytes();
            for(size_t m=0;m<vec.size();m++){
                printf("%02x",vec[vec.size()-m-1]);
            }
            printf(", ");
        }
        printf("\n\n y :\n");
        for(int i=0;i<instance.d;i++){
            vec=sumy[i].to_bytes();
            for(size_t m=0;m<vec.size();m++){
                printf("%02x",vec[vec.size()-m-1]);
            }
            printf(", ");
        }
        printf("\n\n");
}

void afficher_alpha_v_part(std::vector<std::vector<field::GF2E>> alpha_el_shares,
                            std::vector<field::GF2E> v_e_shares, banquet_instance_t &instance){
    for(size_t party=0;party<instance.num_MPC_parties;party++){
        printf("Participant %d :\nalpha : ",party);
        for(size_t l=0;l<alpha_el_shares[party].size();l++){
            printf("%0x, ",alpha_el_shares[party][l]);
        }
        printf("\nv : %0x\n\n",v_e_shares[party]);
    }
}

void afficher_a_c_part(std::vector<std::vector<field::GF2E>> a_el_shares,
                       std::vector<field::GF2E> c_e_shares, size_t N, size_t r){
    for(size_t party=0;party<N;party++){
        printf("Participant %d :\na : ",party);
        for(size_t l=0;l<r;l++){
            printf("%0x, ",a_el_shares[party][l]);
        }
        printf("\nc : %0x\n\n",c_e_shares[party]);
    }                    
}

void afficher_eval_poly_part(std::vector<std::vector<field::GF2E>> X1_poly_eval,
                            std::vector<std::vector<field::GF2E>> X2_poly_eval,
                            std::vector<std::vector<field::GF2E>> Y_poly_eval,
                            std::vector<field::GF2E> R_el,
                            size_t N, size_t number){
    for(size_t party=0;party<N;party++){
        printf("X1(%d) numero %d = ",party,number); X1_poly_eval[party][number].afficher(); 
        printf("   R : "); R_el[number].afficher(); printf("\n\n");
        printf("X2(%d) numero %d = ",party,number); X2_poly_eval[party][number].afficher(); 
        printf("   R : "); R_el[number].afficher(); printf("\n\n");
        printf("Y(%d) numero %d = ",party,number); Y_poly_eval[party][number].afficher(); 
        printf("   R : "); R_el[number].afficher(); printf("\n\n");
    }
}

void afficher_poly_part(std::vector<std::vector<std::vector<field::GF2E>>> X1_poly_shares,
                        std::vector<std::vector<std::vector<field::GF2E>>> X2_poly_shares,
                        std::vector<std::vector<std::vector<field::GF2E>>> Y_poly_shares,
                        size_t N, size_t d, size_t number){
    for(size_t i=0;i<N;i++){
        printf("X1(%d) numero %d = ",i+1,number); afficher_poly(X1_poly_shares[i][number]); printf("\n\n");
        for(int j=0;j<d;j++){
        eval(X1_poly_shares[i][number],field::lift_uint8_t(j)).afficher(); printf(" ");
        }
        printf("\n\n");
        printf("X2(%d) numero %d = ",i+1,number); afficher_poly(X2_poly_shares[i][number]); printf("\n\n");
        for(int j=0;j<d;j++){
        eval(X2_poly_shares[i][number],field::lift_uint8_t(j)).afficher(); printf(" ");
        }
        printf("\n\n");
    }
    for(size_t i=0;i<N;i++){
        printf("Y(%d) numero %d = ",i+1,number); afficher_poly(Y_poly_shares[i][number]); printf("\n\n");
        for(int j=0;j<d;j++){
        eval(Y_poly_shares[i][number],field::lift_uint8_t(j)).afficher(); printf(" ");
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

void afficher_R_el(std::vector<field::GF2E> R_el, size_t r){
    printf("On étend h_1 pour obtenir les R_el :\n");
    for(size_t l=0;l<r;l++){
        R_el[l].afficher(); printf(" ");
    }
    printf("\n\n");
}


void afficher_miss_part(uint16_t missing_party){
    printf("\nOn cache la branche %d.\n",missing_party);
}


void afficher_eval_part(std::vector<std::vector<field::GF2E>> X1_poly_eval,
                        std::vector<std::vector<field::GF2E>> X2_poly_eval,
                        std::vector<std::vector<field::GF2E>> Y_poly_eval,
                        size_t N, size_t r){
    for(size_t party=0;party<N;party++){
        printf("\nX1 de la partie %d :\n",party);
        for(size_t l=0;l<r;l++){
            printf("%0x, ",X1_poly_eval[party][l]);
        }
        printf("\nX2 de la partie %d :\n",party);
        for(size_t l=0;l<r;l++){
            printf("%0x, ",X2_poly_eval[party][l]);
        }
        printf("\nY de la partie %d :\n",party);
        for(size_t l=0;l<r;l++){
            printf("%0x, ",Y_poly_eval[party][l]);
        }
        printf("\n");
    }
}