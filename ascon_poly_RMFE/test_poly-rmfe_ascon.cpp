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



#include "poly-rmfe_ascon.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <math.h>


#define Ecriture_temps 0

int main(){
    std::ofstream monFlux("reccord_time/score_temps.txt", std::ios::app);

    const char *message = "Exemple de texte.";
    banquet_params_t param{Banquet_L1_Param4};
    // On fait num_tours répétitions, avec num_MPC_parties participants
    // On utilise RMFE avec les paramètres vecsize pour la taille de vecteurs, lambda pour la taile du corps
    // On construit num_poly de polynômes interpolés sur les entrées
    uint32_t num_tours,num_MPC_parties,num_verif;
    size_t lambda,vecsize,numbytes,deg_poly,num_poly;
    num_tours=27;
    num_MPC_parties=57;
    lambda=48;
    vecsize=16; num_verif=5760/vecsize; numbytes=vecsize/8;
    deg_poly=24; num_poly=num_verif/deg_poly;
    banquet_instance_t instance{32,16,num_tours,num_MPC_parties,num_verif,deg_poly,num_poly,lambda,vecsize,numbytes,param};

    printf("On choisit de faire %d répétitions du MPC avec %d participants.\n",num_tours,num_MPC_parties);
    printf("On utilise un (%d,%d)2-RMFE pour vérifier les AND.\n",vecsize,lambda);
    printf("On fait la technique de vérification Banquet couplée à du RMFE : les polynômes d'entrée seront de degré %d au nombre de %d.\n",deg_poly-1,num_poly);



    const std::vector<uint8_t> key = {0x0f,0x15,0x71,0xc9,
        0x47,0xd9,0xe8,0x59,
        0x0c,0xb7,0xad,0xd6,
        0xaf,0x7f,0x67,0x98};
    
    const std::vector<uint8_t> key_bad = {0x0f,0x15,0x71,0xc9,
        0x47,0xd9,0xe8,0x59,
        0x0c,0x17,0xad,0xd6,
        0xaf,0x7f,0x67,0x98};

    

    const std::vector<uint8_t> plaintext = {0x02, 0x0c, 0xbd, 0xc8, 0x4e, 0x1d,
                                            0x83, 0x2d, 0x4e, 0x65, 0xd5, 0xc7,
                                            0x15, 0x1e, 0x57, 0x81};
    const std::vector<uint8_t> nonce = {0xa3, 0x61, 0xcc, 0xdd, 0x42, 0x54, 0xf5, 
                                            0x53, 0x1d, 0xbc, 0xcd, 0x0e, 0x9a,
                                            0xf3, 0x52, 0x3c};                                      
    const std::vector<uint8_t> ciphertext_expected = {0x90, 0xf1, 0x90, 0x02, 0x88, 0xb4, 0xac, 0x53,
                                                      0xd2, 0x56, 0x5b, 0x96, 0x4f, 0x13, 0x3a, 0xff};

    banquet_keypair_t keypair;
    keypair.first = key;
    if(TRICHE!=0){
       keypair.first=key_bad;
    }
    keypair.second = plaintext;
    keypair.second.insert(keypair.second.end(), nonce.begin(),
                        nonce.end());
    keypair.second.insert(keypair.second.end(), ciphertext_expected.begin(),
                        ciphertext_expected.end());


    auto start=std::chrono::high_resolution_clock::now();
    banquet_signature_t sign=sign_poly_rmfe_ascon(instance,keypair,(const uint8_t *) message,strlen(message));
    auto stop=std::chrono::high_resolution_clock::now();
    auto duration=std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start);
    printf("\n\nDurée de la signature : %.9f seconds.\n\n",duration.count()*1e-9);


    // écriture dans le fichier score_temps le temps de la signature
    if(Ecriture_temps!=0){    
        if(monFlux){
            monFlux << duration.count()*1e-9 << "  ";
        }
        else{
            std::cout << "ERREUR: Impossible d'ouvrir le fichier." << std::endl;
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
        for(size_t l=0;l<sign.proofs[rep].delta_Y_poly.size();l++){
            sign_size+=sign.proofs[rep].delta_Y_poly[l].to_bytes().size();
        }
        for(size_t l=0;l<sign.proofs[rep].X1_eval.size();l++){
            sign_size+=sign.proofs[rep].X1_eval[l].to_bytes().size();
        }
        for(size_t l=0;l<sign.proofs[rep].X2_eval.size();l++){
            sign_size+=sign.proofs[rep].X2_eval[l].to_bytes().size();
        }
    }
    printf("Taille de la signature : %d octets.\n",sign_size);

    size_t lambda_size=instance.RMFE_fieldsize>>3;
    if((instance.RMFE_fieldsize&7)!=0){
        lambda_size++;
    } 

    size_t sign_size2=0;
    sign_size2=(size_t)(ceil(log2(instance.num_MPC_parties)))*16+48+lambda_size*(instance.d+instance.deg_poly-1+2*instance.num_poly);
    sign_size2*=instance.num_rounds;
    sign_size2+=96;
    if(sign_size2==sign_size){
        printf("Formule de la taille de la signature = 96 + tau * ( log2(N) * 16 + 48 + lambda/8 * ( 5760/vecsize + deg-1 + 2*5760/(deg*vecsize) ) ).\n\n");
    }
    else{
        printf("\n");
    }

    start=std::chrono::high_resolution_clock::now();
    bool verification=verif(instance,keypair.second,(const uint8_t *) message,strlen(message),sign);
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