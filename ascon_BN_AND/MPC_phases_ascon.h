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


#include <inttypes.h>
#include "header_affichage.h"
#include <cassert>
#include <chrono>

#define round_idx 0
#define number 20



ascon_type_t conversion_ascon(std::vector<uint8_t> vect, size_t indice, size_t len){
    size_t size=vect.size();
    if((size>>3)<indice+len){
        throw std::runtime_error("mauvaise dimension");
    }
    else{
        ascon_type_t res(len);
        for(int i=(indice<<3);i<((indice+len)<<3);i++){
            res[(i>>3)-indice][i&7]=vect[i];
        }
        return res;
    }
}



/*-------------- Commit des graines --------------*/

void commit_to_4_party_seeds(
    const banquet_instance_t &instance, const gsl::span<uint8_t> &seed0,
    const gsl::span<uint8_t> &seed1, const gsl::span<uint8_t> &seed2,
    const gsl::span<uint8_t> &seed3, const banquet_salt_t &salt, size_t rep_idx,
    size_t party_idx, gsl::span<uint8_t> com0, gsl::span<uint8_t> com1,
    gsl::span<uint8_t> com2, gsl::span<uint8_t> com3) {
  hash_context_x4 ctx;
  hash_init_x4(&ctx, instance.digest_size);
  hash_update_x4_1(&ctx, salt.data(), salt.size());
  hash_update_x4_uint16_le(&ctx, (uint16_t)rep_idx);
  const uint16_t party_idxs[4] = {
      (uint16_t)party_idx, (uint16_t)(party_idx + 1), (uint16_t)(party_idx + 2),
      (uint16_t)(party_idx + 3)};
  hash_update_x4_uint16s_le(&ctx, party_idxs);
  hash_update_x4_4(&ctx, seed0.data(), seed1.data(), seed2.data(), seed3.data(),
                   instance.seed_size);
  hash_final_x4(&ctx);

  hash_squeeze_x4_4(&ctx, com0.data(), com1.data(), com2.data(), com3.data(),
                    instance.digest_size);
}

void commit_to_party_seed(const banquet_instance_t &instance,
                          const gsl::span<uint8_t> &seed,
                          const banquet_salt_t &salt, size_t rep_idx,
                          size_t party_idx, gsl::span<uint8_t> commitment) {
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, (uint16_t)rep_idx);
  hash_update_uint16_le(&ctx, (uint16_t)party_idx);
  hash_update(&ctx, seed.data(), seed.size());
  hash_final(&ctx);

  hash_squeeze(&ctx, commitment.data(), commitment.size());
}




/*-----------------------------PHASES POUR SIGNATURE-----------------------------*/

/*------------------ PHASE 1 ------------------*/

std::vector<uint8_t> phase_1_MPC_compute(const uint8_t *msg, size_t msg_len, banquet_salt_t salt, const RepByteContainer &commitments,
                        banquet_instance_t &instance, RandomTapes random_tapes, banquet_keypair_t key,
                        size_t *indice, std::vector<std::vector<uint8_t>> *delta_key,
                        std::vector<std::vector<std::vector<uint8_t>>> *x1_shares,
                        std::vector<std::vector<std::vector<uint8_t>>> *x2_shares,
                        std::vector<std::vector<std::vector<uint8_t>>> *y_shares,
                        std::vector<std::vector<uint8_t>> *delta_y, 
                        std::vector<std::vector<std::vector<uint8_t>>> *a_el_shares,
                        std::vector<std::vector<std::vector<uint8_t>>> *c_e_shares,
                        std::vector<std::vector<uint8_t>> *delta_c){           
    
    ascon_type_t nonce_in,plaintext_in;
    plaintext_in=conversion_ascon(key.second,0,2);
    nonce_in=conversion_ascon(key.second,2,2);

    
    // Séparation de la clé secrète :
    *indice=0;
    size_t tmp;
    ascon_type_t K=conversion_ascon(key.first,0,2);
    std::vector<std::vector<uint8_t>> key_shares(instance.num_MPC_parties);
    (*delta_key)=std::vector<std::vector<uint8_t>>(instance.num_rounds);
    std::vector<std::vector<ascon_type_t>> key_in(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        tmp=*indice;
        key_in[rep]=Share2(K,random_tapes,instance.num_MPC_parties,rep,indice,&key_shares,&((*delta_key)[rep]));
        if(rep!=instance.num_rounds-1){
            *indice=tmp;
        }
    }

    // Affichage des partitions de clé et sa recomposition :
    // afficher_key_part(key_in[round_idx],instance.num_MPC_parties,1);

    // Calcul MPC d'Ascon :
    std::vector<std::vector<ascon_type_t>> ciphertext_out(instance.num_rounds);
    (*x1_shares) = std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);
    (*x2_shares) = std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);
    (*y_shares) = std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);
    (*delta_y)=std::vector<std::vector<uint8_t>>(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        ciphertext_out[rep]=std::vector<ascon_type_t>(instance.num_MPC_parties);
        (*x1_shares)[rep] = std::vector<std::vector<uint8_t>>(instance.num_MPC_parties);
        (*x2_shares)[rep] = std::vector<std::vector<uint8_t>>(instance.num_MPC_parties);
        (*y_shares)[rep] = std::vector<std::vector<uint8_t>>(instance.num_MPC_parties);
        tmp=*indice;
        ascon_128_groupe_prouveur(&(key_in[rep]),&plaintext_in,&nonce_in,&(ciphertext_out[rep]),&(*x1_shares)[rep],&(*x2_shares)[rep],&(*y_shares)[rep],
                                    &((*delta_y)[rep]),random_tapes,rep,indice);
        if(rep!=instance.num_rounds-1){
            *indice=tmp;
        }
    }

    // Affichage des partitions de la sortie et sa recomposition :
    // afficher_ciphertext_part(ciphertext_out[round_idx],instance.num_MPC_parties,1);


    // Affichage des parties des parties d'entrées et sortie des portes AND :
    // afficher_in_out_part((*x1_shares)[round_idx],(*x2_shares)[round_idx],(*y_shares)[round_idx],instance.num_MPC_parties);


    // Affichage des entrées et sorties recomposées :
    // afficher_in_out_rec((*x1_shares)[round_idx],(*x2_shares)[round_idx],(*y_shares)[round_idx],instance.num_MPC_parties);


    // Création des a_el et des c_e :
    gsl::span<uint8_t> tape;

    (*delta_c)=std::vector<std::vector<uint8_t>>(instance.num_rounds);
    std::vector<std::vector<uint8_t>> c_voulu(instance.num_rounds);


    std::vector<std::vector<uint8_t>> a_el(instance.num_rounds);
    (*a_el_shares)=std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);
    (*c_e_shares)=std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);

    for(size_t rep=0;rep<instance.num_rounds;rep++){
        (*a_el_shares)[rep]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties);
        (*c_e_shares)[rep]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties);
        (*delta_c)[rep]=std::vector<uint8_t>(8);
        a_el[rep]=std::vector<uint8_t>(720);
        for(size_t i=0;i<instance.num_MPC_parties;i++){
            (*a_el_shares)[rep][i]=std::vector<uint8_t>(720);
            (*c_e_shares)[rep][i]=std::vector<uint8_t>(8);
            tape=random_tapes.get_bytes(rep,i,*indice,728);
            for(size_t l=0;l<720;l++){
                (*a_el_shares)[rep][i][l]=tape[l];
                a_el[rep][l]^=(*a_el_shares)[rep][i][l];
            }
            for(size_t b=0;b<8;b++){
                (*c_e_shares)[rep][i][b]=tape[720+b];
                (*delta_c)[rep][b]^=(*c_e_shares)[rep][i][b];
            }
        }
    }

    std::vector<uint8_t> tmp_x2(8);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        c_voulu[rep]=std::vector<uint8_t>(8);
        for(size_t l=0;l<90;l++){
            for(size_t j=0;j<8;j++){
                tmp_x2[j]=0;
                for(size_t i=0;i<instance.num_MPC_parties;i++){
                    tmp_x2[j]^=(*x2_shares)[rep][i][(l<<3) | j];
                }
                c_voulu[rep][j]^=(a_el[rep][(l<<3) | j] & tmp_x2[j]);
            }
        }
        for(size_t j=0;j<8;j++){
            (*delta_c)[rep][j]^=c_voulu[rep][j];
            (*c_e_shares)[rep][0][j]^=(*delta_c)[rep][j];
        }
    }

    *indice+=728;


    
    // Faire le commit de toutes les données :
    hash_context ctx;
    hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_1);
    hash_update(&ctx, salt.data(), salt.size());
    hash_update(&ctx, key.second.data(), key.second.size());
    hash_update(&ctx, msg, msg_len);

    for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
        for (size_t party = 0; party < instance.num_MPC_parties; party++) {
            auto commitment = commitments.get(repetition, party);
            hash_update(&ctx, commitment.data(), commitment.size());
            auto output_broadcast = ciphertext_out[repetition][party];
            for(size_t j=0;j<output_broadcast.size();j++){
                hash_update(&ctx, output_broadcast[j].data(), output_broadcast[j].size());
            }
        }
        hash_update(&ctx, (*delta_key)[repetition].data(),
                    (*delta_key)[repetition].size());
        hash_update(&ctx, (*delta_c)[repetition].data(), 8);
        hash_update(&ctx, (*delta_y)[repetition].data(), 720);
    }
    hash_final(&ctx);

    std::vector<uint8_t> commitment(instance.digest_size);
    hash_squeeze(&ctx, commitment.data(), commitment.size());
    return commitment;
}





std::vector<std::vector<uint8_t>>
phase_1_expand(const banquet_instance_t &instance,
               const std::vector<uint8_t> &h_1) {
    hash_context ctx;
    hash_init(&ctx, instance.digest_size);
    hash_update(&ctx, h_1.data(), h_1.size());
    hash_final(&ctx);

    std::vector<std::vector<uint8_t>> eps_el(instance.num_rounds);
    for (size_t e = 0; e < instance.num_rounds; e++) {
        eps_el[e]=std::vector<uint8_t>(720);
        hash_squeeze(&ctx, eps_el[e].data(), 720);
    }
    return eps_el;
}





/*------------------ PHASE 2 ------------------*/


std::vector<uint8_t> phase_2_sacrifycing(const std::vector<uint8_t> &h_1, const banquet_salt_t &salt, std::vector<std::vector<std::vector<uint8_t>>> x1_shares,
                                        std::vector<std::vector<std::vector<uint8_t>>> x2_shares, std::vector<std::vector<std::vector<uint8_t>>> y_shares,
                                        std::vector<std::vector<std::vector<uint8_t>>> a_el_shares, 
                                        std::vector<std::vector<std::vector<uint8_t>>> c_e_shares, 
                                        banquet_instance_t &instance, std::vector<std::vector<uint8_t>> *alpha_el, 
                                        std::vector<std::vector<uint8_t>> epsilon_el){

    (*alpha_el)=std::vector<std::vector<uint8_t>>(instance.num_rounds);
    std::vector<std::vector<std::vector<uint8_t>>> alpha_el_shares(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        alpha_el_shares[rep]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties);
        (*alpha_el)[rep]=std::vector<uint8_t>(720);
        for(size_t party=0;party<instance.num_MPC_parties;party++){
            alpha_el_shares[rep][party]=std::vector<uint8_t>(720);            
            for(size_t l=0;l<720;l++){
                alpha_el_shares[rep][party][l]=(epsilon_el[rep][l] & x1_shares[rep][party][l]) ^ a_el_shares[rep][party][l];
                (*alpha_el)[rep][l]^=alpha_el_shares[rep][party][l];
            }
        }
    }

    std::vector<std::vector<std::vector<uint8_t>>> v_e_shares(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        v_e_shares[rep]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties);
        for(size_t party=0;party<instance.num_MPC_parties;party++){    
            v_e_shares[rep][party]=std::vector<uint8_t>(8);
            for(size_t j=0;j<8;j++){    
                v_e_shares[rep][party][j]=c_e_shares[rep][party][j];
                for(size_t l=0;l<90;l++){
                    v_e_shares[rep][party][j]^=((*alpha_el)[rep][(l<<3) | j] & x2_shares[rep][party][(l<<3) | j])^
                                                (epsilon_el[rep][(l<<3) | j] & y_shares[rep][party][(l<<3) | j]);
                }
            }
            
        }
    }


    // Affichage des alpha_el et v_e en parties :
    // afficher_alpha_v_part(alpha_el_shares[round_idx],v_e_shares[round_idx],instance);


    hash_context ctx;
    hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_3);
    hash_update(&ctx, salt.data(), salt.size());
    hash_update(&ctx, h_1.data(), h_1.size());
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        for(size_t party=0;party<instance.num_MPC_parties;party++){
            hash_update(&ctx,alpha_el_shares[rep][party].data(),720);
            hash_update(&ctx,v_e_shares[rep][party].data(),8);
        }
    }
    hash_final(&ctx);

    std::vector<uint8_t> commitment(instance.digest_size);
    hash_squeeze(&ctx, commitment.data(), commitment.size());
    return commitment;
}




std::vector<uint16_t> phase_2_expand(const banquet_instance_t &instance,
                                     const std::vector<uint8_t> &h_2) {
  assert(instance.num_MPC_parties < (1ULL << 16));
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, h_2.data(), h_2.size());
  hash_final(&ctx);
  size_t num_squeeze_bytes = instance.num_MPC_parties > 256 ? 2 : 1;

  std::vector<uint16_t> opened_parties;
  uint16_t mask = (1ULL << ceil_log2(instance.num_MPC_parties)) - 1;
  for (size_t e = 0; e < instance.num_rounds; e++) {
    uint16_t party;
    do {
      hash_squeeze(&ctx, (uint8_t *)&party, num_squeeze_bytes);
      party = le16toh(party);
      party = party & mask;
    } while (party >= instance.num_MPC_parties);
    opened_parties.push_back(party);
  }
  return opened_parties;
}














/*-----------------------------LES COMMITS h_i-----------------------------*/

std::vector<uint8_t> h_1_commit(banquet_instance_t &instance, std::vector<uint8_t> pk, const uint8_t *msg, size_t msg_len, 
                                banquet_signature_t sign, RepByteContainer party_seed_commitments, 
                                std::vector<std::vector<ascon_type_t>> ciphertext_out){
    hash_context ctx;
    hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_1);
    hash_update(&ctx, sign.salt.data(), sign.salt.size());
    hash_update(&ctx, pk.data(), pk.size());
    hash_update(&ctx, msg, msg_len);

    for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
        for (size_t party = 0; party < instance.num_MPC_parties; party++) {
            auto commitment = party_seed_commitments.get(repetition, party);
            hash_update(&ctx, commitment.data(), commitment.size());
            auto output_broadcast = ciphertext_out[repetition][party];
            for(size_t j=0;j<output_broadcast.size();j++){
                hash_update(&ctx, output_broadcast[j].data(), output_broadcast[j].size());
            }
        }
        hash_update(&ctx, sign.proofs[repetition].delta_key.data(),
                    sign.proofs[repetition].delta_key.size());
        hash_update(&ctx, sign.proofs[repetition].delta_c.data(), 8);
        hash_update(&ctx, sign.proofs[repetition].delta_y.data(), 720);
    }
    hash_final(&ctx);

    std::vector<uint8_t> h_1(instance.digest_size);
    hash_squeeze(&ctx, h_1.data(), h_1.size());
    return h_1;
}


std::vector<uint8_t> h_2_commit(banquet_instance_t &instance, banquet_signature_t sign,
                                std::vector<std::vector<std::vector<uint8_t>>> alpha_el_shares, 
                                std::vector<std::vector<std::vector<uint8_t>>> v_e_shares){
    hash_context ctx;
    hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_3);
    hash_update(&ctx, sign.salt.data(), sign.salt.size());
    hash_update(&ctx, sign.h_1.data(), sign.h_1.size());
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        for(size_t party=0;party<instance.num_MPC_parties;party++){
            hash_update(&ctx,alpha_el_shares[rep][party].data(),720);
            hash_update(&ctx,v_e_shares[rep][party].data(),8);
        }
    }
    hash_final(&ctx);

    std::vector<uint8_t> commitment(instance.digest_size);
    hash_squeeze(&ctx, commitment.data(), commitment.size());
    return commitment;
}












/*-----------------------------PHASES POUR VÉRIFICATION-----------------------------*/


void calcul_MPC_ascon(banquet_signature_t sign, banquet_instance_t &instance, std::vector<uint8_t> pk, 
                        RandomTapes random_tapes, std::vector<size_t> N_ind, size_t *indice,
                        std::vector<std::vector<ascon_type_t>> *ciphertext_out, 
                        std::vector<std::vector<std::vector<uint8_t>>> *x1_shares, std::vector<std::vector<std::vector<uint8_t>>> *x2_shares,
                        std::vector<std::vector<std::vector<uint8_t>>> *y_shares){
    ascon_type_t nonce_in,plaintext_in;
    plaintext_in=conversion_ascon(pk,0,2);
    nonce_in=conversion_ascon(pk,2,2);
    
    // Séparation de la clé secrète :
    size_t delta_ind=0;
    std::pair<size_t,size_t> tmp;
    std::vector<std::vector<uint8_t>> key_shares(instance.num_MPC_parties-1);
    std::vector<uint8_t> delta_key;
    std::vector<std::vector<ascon_type_t>> key_in(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        tmp.first=*indice;
        tmp.second=delta_ind;
        delta_key=sign.proofs[rep].delta_key;
        key_in[rep]=Share_verif2(2,random_tapes,instance.num_MPC_parties-1,N_ind[rep],rep,indice,&key_shares,delta_key,&delta_ind);
        if(rep!=instance.num_rounds-1){
            *indice=tmp.first;
            delta_ind=tmp.second;
        }
    }

    // Affichage des partitions de clés :
    // afficher_key_part(key_in[round_idx],instance.num_MPC_parties-1,0);

    // Calcul MPC d'Ascon :
    (*ciphertext_out)=std::vector<std::vector<ascon_type_t>>(instance.num_rounds);
    (*x1_shares)=std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);
    (*x2_shares)=std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);
    (*y_shares)=std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);
    std::vector<uint8_t> delta_y;
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        (*ciphertext_out)[rep]=std::vector<ascon_type_t>(instance.num_MPC_parties);
        (*x1_shares)[rep]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties-1);
        (*x2_shares)[rep]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties-1);
        (*y_shares)[rep]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties-1);
        tmp.first=*indice;
        delta_y=sign.proofs[rep].delta_y;
        ascon_128_groupe_verif(&(key_in[rep]),&plaintext_in,&nonce_in,&((*ciphertext_out)[rep]),&((*x1_shares)[rep]),&((*x2_shares)[rep]),&((*y_shares)[rep]),
                                    delta_y,random_tapes,N_ind[rep],rep,indice);
        (*ciphertext_out)[rep][N_ind[rep]]=Rec2({Rec2((*ciphertext_out)[rep]),conversion_ascon(pk,4,2)});
        if(rep!=instance.num_rounds-1){
            *indice=tmp.first;
        }
    }

    // Affichage des parties des parties d'entrées et sortie des portes AND :
    // afficher_in_out_part((*x1_shares)[round_idx],(*x2_shares)[round_idx],(*y_shares)[round_idx],instance.num_MPC_parties-1);



    // Affichage des partitions de la sortie :
    // afficher_ciphertext_part((*ciphertext_out)[round_idx],instance.num_MPC_parties,0);
}



void creation_a_c(RandomTapes random_tapes, size_t *indice, std::vector<std::vector<std::vector<uint8_t>>> *a_el_shares, 
                    std::vector<std::vector<std::vector<uint8_t>>> *c_e_shares, banquet_signature_t sign,
                    banquet_instance_t &instance, std::vector<size_t> N_ind){
    gsl::span<uint8_t> tape;

    (*a_el_shares)=std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);
    (*c_e_shares)=std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);

    for(size_t rep=0;rep<instance.num_rounds;rep++){
        (*a_el_shares)[rep]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties-1);
        (*c_e_shares)[rep]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties-1);
        for(size_t i=0;i<instance.num_MPC_parties-1;i++){
            if(i<N_ind[rep]){
                tape=random_tapes.get_bytes(rep,i,*indice,728);
            }
            else{
                tape=random_tapes.get_bytes(rep,i+1,*indice,728);
            }
            (*a_el_shares)[rep][i]=std::vector<uint8_t>(720);
            (*c_e_shares)[rep][i]=std::vector<uint8_t>(8);
            for(size_t l=0;l<720;l++){
                (*a_el_shares)[rep][i][l]=tape[l];
            }
            for(size_t b=0;b<8;b++){
                (*c_e_shares)[rep][i][b]=tape[720+b];
            }
        }
    }

    for(size_t rep=0;rep<instance.num_rounds;rep++){
        if(N_ind[rep]!=0){
            for(size_t b=0;b<8;b++){
                (*c_e_shares)[rep][0][b]^=sign.proofs[rep].delta_c[b];
            }
        }
    }
    
    *indice+=728;
}


void calcul_alpha_v_shares(std::vector<std::vector<std::vector<uint8_t>>> x1_shares, std::vector<std::vector<std::vector<uint8_t>>> x2_shares,
                            std::vector<std::vector<std::vector<uint8_t>>> y_shares, std::vector<std::vector<std::vector<uint8_t>>> a_el_shares, 
                            std::vector<std::vector<std::vector<uint8_t>>> c_e_shares, banquet_instance_t &instance, 
                            std::vector<std::vector<std::vector<uint8_t>>> *alpha_el_shares, std::vector<std::vector<std::vector<uint8_t>>> *v_e_shares,
                            std::vector<std::vector<uint8_t>> epsilon_el,
                            banquet_signature_t sign, std::vector<size_t> N_ind){

    int partynew;
    (*alpha_el_shares)=std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);
    (*v_e_shares)=std::vector<std::vector<std::vector<uint8_t>>>(instance.num_rounds);
    for(size_t repetition=0;repetition<instance.num_rounds;repetition++){
        (*alpha_el_shares)[repetition]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties);
        (*v_e_shares)[repetition]=std::vector<std::vector<uint8_t>>(instance.num_MPC_parties);
        (*alpha_el_shares)[repetition][N_ind[repetition]]=sign.proofs[repetition].alpha_l;
        (*v_e_shares)[repetition][N_ind[repetition]]=std::vector<uint8_t>(8);
        for(size_t party=0;party<instance.num_MPC_parties-1;party++){
            if(party<N_ind[repetition]){
                partynew=party;
            }
            else{
                partynew=party+1;
            }
            (*alpha_el_shares)[repetition][partynew]=std::vector<uint8_t>(720);
            (*v_e_shares)[repetition][partynew]=c_e_shares[repetition][party];
            for(size_t j=0;j<8;j++){
                for(size_t l=0;l<90;l++){
                    (*alpha_el_shares)[repetition][partynew][(l<<3) | j]=(epsilon_el[repetition][(l<<3) | j] & x1_shares[repetition][party][(l<<3) | j]) 
                                                                        ^ (a_el_shares[repetition][party][(l<<3) | j]);
                    (*v_e_shares)[repetition][partynew][j]^=(sign.proofs[repetition].alpha_l[(l<<3) | j] & x2_shares[repetition][party][(l<<3) | j]) 
                                                            ^ (epsilon_el[repetition][(l<<3) | j] & y_shares[repetition][party][(l<<3) | j]);
                    (*alpha_el_shares)[repetition][N_ind[repetition]][(l<<3) | j]^=(*alpha_el_shares)[repetition][partynew][(l<<3) | j];
                }
                (*v_e_shares)[repetition][N_ind[repetition]][j]^=(*v_e_shares)[repetition][partynew][j];
            }
        }
    }

    // Affichage des alpha_el et v_e en parties :
    // afficher_alpha_v_part((*alpha_el_shares)[round_idx],(*v_e_shares)[round_idx],instance);
}

