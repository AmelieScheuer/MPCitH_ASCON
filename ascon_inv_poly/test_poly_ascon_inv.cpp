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



#include "poly_ascon_inv.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <math.h>


#define Ecriture_temps 0


banquet_keypair_t key_gen(std::vector<uint8_t> key, std::vector<uint8_t> plaintext, std::vector<uint8_t> nonce, taille_lambda lambda){
    ascon_type_t Key(2),Plaintext(2),Nonce(2);
    for(size_t i=0;i<2;i++){
        for(size_t j=0;j<8;j++){
            Key[i][j]=key[(i<<3)+j];
            Plaintext[i][j]=plaintext[(i<<3)+j];
            Nonce[i][j]=nonce[(i<<3)+j];
        }
    }
    // printf("IV :\n"); afficher_type(IV);
    // printf("\nclé secrète K :\n"); afficher_type(Key);
    // printf("\ncomplément N :\n"); afficher_type(Nonce);
    // printf("\nLe message clair P :\n"); afficher_type(Plaintext); printf("\n\n");

    ascon_type_t Ciphertext=ascon_128(Key,Nonce,Plaintext,lambda);
    // printf("ASCON_inv(K,N,P) avec lambda = %d donne comme résultat :\n",lambda.lambda); afficher_type(Ciphertext); printf("\n\n");
    std::vector<uint8_t> ciphertext_expected(16);
    for(size_t i=0;i<2;i++){
        for(size_t j=0;j<8;j++){
            ciphertext_expected[(i<<3)+j]=Ciphertext[i][j];
        }
    }

    banquet_keypair_t keypair;
    keypair.first = key;
    keypair.second = plaintext;
    keypair.second.insert(keypair.second.end(), nonce.begin(), nonce.end());
    keypair.second.insert(keypair.second.end(), ciphertext_expected.begin(), ciphertext_expected.end());

    return keypair;
}


#define MESSAGE_DEFAULT "Exemple de texte."

int main(int argc, char *argv[]){
    std::ofstream monFlux("reccord_time/score_temps.txt", std::ios::app);

   std::string message = argc == 1 ? MESSAGE_DEFAULT : argv[1];
    banquet_params_t param{Banquet_L1_Param4};
    // On fait num_tours répétitions, avec num_MPC_parties participants
    // lambda*5 : la taile du corps
    // On construit num_poly de polynômes interpolés sur les entrées
    uint32_t num_tours,num_MPC_parties,num_verif;
    size_t vecsize,numbytes,deg_poly,num_poly;
    taille_lambda lambda;

    num_tours=28;
    num_MPC_parties=57;
    init_lambda(8,&lambda);
    deg_poly=12; num_poly=(size_t) ceil((float)lambda.non_lin_op_number/deg_poly);
    banquet_instance_t instance{32,16,num_tours,num_MPC_parties,lambda,deg_poly,num_poly,param};

    printf("On choisit de faire %d répétitions du MPC avec %d participants.\n",num_tours,num_MPC_parties);
    printf("version d'ASCON avec lambda = %d.\n",lambda.lambda);
    printf("nombre d'opérations non linéaires : %d\n",lambda.non_lin_op_number);
    printf("On fait la technique de vérification Banquet : les polynômes d'entrée seront de degré %d au nombre de %d.\n",deg_poly-1,num_poly);



    const std::vector<uint8_t> key = {0x0f,0x15,0x71,0xc9,
        0x47,0xd9,0xe8,0x59,
        0x0c,0xb7,0xad,0xd6,
        0xaf,0x7f,0x67,0x98};
    const std::vector<uint8_t> plaintext = {0x02, 0x0c, 0xbd, 0xc8, 0x4e, 0x1d,
                                            0x83, 0x2d, 0x4e, 0x65, 0xd5, 0xc7,
                                            0x15, 0x1e, 0x57, 0x81};
    const std::vector<uint8_t> nonce = {0xa3, 0x61, 0xcc, 0xdd, 0x42, 0x54, 0xf5, 
                                            0x53, 0x1d, 0xbc, 0xcd, 0x0e, 0x9a,
                                            0xf3, 0x52, 0x3c};  

    field::GF2E::init_extension_field(instance);
    banquet_keypair_t keypair=key_gen(key,plaintext,nonce,lambda);


    auto start=std::chrono::high_resolution_clock::now();
    banquet_signature_t sign=sign_poly_rmfe_ascon(instance,keypair,(const uint8_t *) message.c_str(), message.size());
    auto stop=std::chrono::high_resolution_clock::now();
    auto duration=std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start);
    printf("\n\nDurée de la signature : %.9f seconds.\n\n",duration.count()*1e-9);

    // écriture dans le fichier score_temps le temps de la signature
    if(Ecriture_temps!=0){
        if(monFlux){
            monFlux << duration.count()*1e-9 << "  ";
        }
        else{
            std::cerr << "ERREUR: Impossible d'ouvrir le fichier." << std::endl;
        }
    }

    size_t sign_size=0;
    
    sign_size= sign.salt.size() + sign.h_1.size() + sign.h_3.size();
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        for(size_t l=0;l<sign.proofs[rep].reveallist.first.size();l++){
            sign_size+=sign.proofs[rep].reveallist.first[l].size();
        }
        sign_size+=sign.proofs[rep].Commit.size();
        sign_size+=sign.proofs[rep].delta_key.size();
        for(size_t l=0;l<sign.proofs[rep].delta_Y.size();l++){
            sign_size+=sign.proofs[rep].delta_Y[l].to_bytes().size();
        }
        for(size_t l=0;l<sign.proofs[rep].delta_Z_poly_eval.size();l++){
            sign_size+=sign.proofs[rep].delta_Z_poly_eval[l].to_bytes().size();
        }
        for(size_t l=0;l<sign.proofs[rep].X_eval.size();l++){
            sign_size+=sign.proofs[rep].X_eval[l].to_bytes().size();
        }
        for(size_t l=0;l<sign.proofs[rep].Y_eval.size();l++){
            sign_size+=sign.proofs[rep].Y_eval[l].to_bytes().size();
        }
    }
    printf("Taille de la signature : %d octets.\n",sign_size);


    size_t sign_size2=0;
    sign_size2=(size_t)(ceil(log2(instance.num_MPC_parties)))*16+48+
                lambda.byte_size*(lambda.non_lin_op_number+instance.deg_poly-1+2*instance.num_poly);
    sign_size2*=instance.num_rounds;
    sign_size2+=96;

    if(sign_size2==sign_size){
        printf("Formule de la taille de la signature = 96 + tau * ( log2(N) * 16 + 48 + 5*lambda/8 * ( 18*64/lambda + deg + 2 * 18*64/(lambda*deg) ) ).\n\n");
    }
    else{
        printf("\n");
    }

    start=std::chrono::high_resolution_clock::now();
    bool verification=verif(instance,keypair.second,(const uint8_t *) message.c_str(), message.size(), sign);
    stop=std::chrono::high_resolution_clock::now();
    duration=std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start);
    printf("\n\nDurée de la vérification : %.9f seconds.\n\n",duration.count()*1e-9);

    // écriture dans le fichier score_temps le temps de la vérification
    if(Ecriture_temps!=0){    
        if(monFlux){
            monFlux << duration.count()*1e-9 << std::endl;
        }
        else{
            std::cout << "ERREUR: Impossible d'ouvrir le fichier." << std::endl;
        }
    }
}