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



#include "header_affichage.h"
#include <cassert>
#include <chrono>

#define round_idx 0
#define number 20



std::vector<std::vector<field::GF2E>> Share_poly_eval(std::vector<field::GF2E> poly_evals,
                                                      RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, size_t byte_size,
                                                      std::vector<field::GF2E> *delta_Y_poly){
    std::vector<std::vector<field::GF2E>> res(N);
    gsl::span<uint8_t> tape;
    std::vector<field::GF2E> sum=poly_evals;
    uint64_t tmp;
    for(int i=0;i<N;i++){
        tape=random_tapes.get_bytes(rep,i,*indice,byte_size*poly_evals.size());
        res[i]=std::vector<field::GF2E>(poly_evals.size());
        for(int j=0;j<poly_evals.size();j++){
          tmp=0;
            for(int z=0;z<byte_size;z++){
              tmp<<=8;
              tmp^=tape[j*byte_size+z];
            }
            res[i][j]=field::lift_uint8_t(tmp);
            sum[j]+=res[i][j];
        }
    }

    for(int j=0;j<poly_evals.size();j++){
        delta_Y_poly->push_back(sum[j]);
        res[0][j]+=sum[j];
    }

    *indice+=byte_size*poly_evals.size();
    return res;
}


std::vector<field::GF2E> evals_d(std::vector<field::GF2E> Y_poly, size_t d){
  std::vector<field::GF2E> poly_evals(d-1);
  for(int i=0;i<d-1;i++){
    poly_evals[i]=field::eval(Y_poly,field::GF2E(d+i));
  }
  return poly_evals;
}



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

inline void hash_update_GF2E(hash_context *ctx,
                             size_t byte_size,
                             const field::GF2E &element) {
  // 8 bytes is enough for supported field sizes
  std::array<uint8_t, 8> buffer;
  element.to_bytes(buffer.data());
  hash_update(ctx, buffer.data(), byte_size);
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
                        std::vector<std::vector<field::GF2E>> *delta_Y,
                        std::vector<std::vector<field::GF2E>> *delta_Y_poly,
                        std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> *X1_poly_shares, 
                        std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> *X2_poly_shares, 
                        std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> *Y_poly_shares){
    size_t lambda_size=instance.RMFE_fieldsize>>3;
    if((instance.RMFE_fieldsize&7)!=0){
        lambda_size++;
    }                        
    
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
        tmp=(*indice);
        key_in[rep]=Share2(K,random_tapes,instance.num_MPC_parties,rep,indice,&key_shares,&((*delta_key)[rep]));
        if(rep!=instance.num_rounds-1){
            *indice=tmp;
        }
    }


    // Affichage des partitions de clé et sa recomposition :
    // afficher_key_part(key_in[round_idx],instance.num_MPC_parties,1);

    // Calcul MPC d'Ascon :
    std::vector<std::vector<ascon_type_t>> ciphertext_out(instance.num_rounds);
    std::vector<std::vector<std::vector<field::GF2E>>> X1_shares(instance.num_rounds);
    std::vector<std::vector<std::vector<field::GF2E>>> X2_shares(instance.num_rounds);
    std::vector<std::vector<std::vector<field::GF2E>>> Y_shares(instance.num_rounds);
    (*delta_Y)=std::vector<std::vector<field::GF2E>>(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        ciphertext_out[rep]=std::vector<ascon_type_t>(instance.num_MPC_parties);
        X1_shares[rep] = std::vector<std::vector<field::GF2E>>(instance.num_MPC_parties);
        X2_shares[rep] = std::vector<std::vector<field::GF2E>>(instance.num_MPC_parties);
        Y_shares[rep] = std::vector<std::vector<field::GF2E>>(instance.num_MPC_parties);
        tmp=*indice;
        ascon_128_groupe_prouveur(&(key_in[rep]),&plaintext_in,&nonce_in,&(ciphertext_out[rep]),&X1_shares[rep],&X2_shares[rep],&Y_shares[rep],
                                    &((*delta_Y)[rep]),random_tapes,rep,indice,instance.RMFE_vecsize,lambda_size);
        if(rep!=instance.num_rounds-1){
            *indice=tmp;
        }
    }

    // Affichage des partitions de la sortie et sa recomposition :
    // afficher_ciphertext_part(ciphertext_out[round_idx],instance.num_MPC_parties,1);


    // Affichage des parties d'entrées et sortie des portes AND :
    // afficher_in_out_part(X1_shares[round_idx],X2_shares[round_idx],Y_shares[round_idx],instance.num_MPC_parties,instance.d);


    // Affichage des entrées et sorties recomposées :
    // afficher_in_out_rec(X1_shares[round_idx],X2_shares[round_idx],Y_shares[round_idx],instance);
    

    (*delta_Y_poly)=std::vector<std::vector<field::GF2E>>(instance.num_rounds);
    
    size_t deg_poly=instance.deg_poly;
    size_t num_poly=instance.num_poly;

    std::vector<field::GF2E> x_values_for_interpolation_zero_to_d =
      field::get_first_n_field_elements(deg_poly);
    std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_d =
      field::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_d);
    std::vector<field::GF2E> x_values_for_interpolation_zero_to_2d =
      field::get_first_n_field_elements(2 * deg_poly - 1);
    std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_2d =
      field::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_2d);

    

    std::vector<std::vector<std::vector<field::GF2E>>> X1_poly(instance.num_rounds),X2_poly(instance.num_rounds);
    
    
    (*X1_poly_shares)=std::vector<std::vector<std::vector<std::vector<field::GF2E>>>>(instance.num_rounds);
    (*X2_poly_shares)=std::vector<std::vector<std::vector<std::vector<field::GF2E>>>>(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        (*X1_poly_shares)[rep]=std::vector<std::vector<std::vector<field::GF2E>>>(instance.num_MPC_parties);
        (*X2_poly_shares)[rep]=std::vector<std::vector<std::vector<field::GF2E>>>(instance.num_MPC_parties);
        X1_poly[rep]=std::vector<std::vector<field::GF2E>>(num_poly);
        X2_poly[rep]=std::vector<std::vector<field::GF2E>>(num_poly);
    for(size_t i=0;i<instance.num_MPC_parties;i++){
        (*X1_poly_shares)[rep][i]=std::vector<std::vector<field::GF2E>>(num_poly);
        (*X2_poly_shares)[rep][i]=std::vector<std::vector<field::GF2E>>(num_poly);
        for(size_t l=0;l<num_poly;l++){
        (*X1_poly_shares)[rep][i][l]=std::vector<field::GF2E>(deg_poly);
        (*X2_poly_shares)[rep][i][l]=std::vector<field::GF2E>(deg_poly);
        for (size_t k = 0; k < deg_poly; k++) {
            (*X1_poly_shares)[rep][i][l][k] = X1_shares[rep][i][l*deg_poly+k];
            (*X2_poly_shares)[rep][i][l][k] = X2_shares[rep][i][l*deg_poly+k];
            }
        (*X1_poly_shares)[rep][i][l] = field::interpolate_with_precomputation(
            precomputation_for_zero_to_d, (*X1_poly_shares)[rep][i][l]);
        (*X2_poly_shares)[rep][i][l] = field::interpolate_with_precomputation(
            precomputation_for_zero_to_d, (*X2_poly_shares)[rep][i][l]);
        if(i==0){
            X1_poly[rep][l]=std::vector<field::GF2E>(deg_poly);
            X2_poly[rep][l]=std::vector<field::GF2E>(deg_poly);
        }
        X1_poly[rep][l]+=(*X1_poly_shares)[rep][i][l];
        X2_poly[rep][l]+=(*X2_poly_shares)[rep][i][l];
        }
    }
    }


    
    std::vector<std::vector<std::vector<field::GF2E>>> Y_poly(instance.num_rounds);
    std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> Y_poly_eval_shares(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        Y_poly[rep]=std::vector<std::vector<field::GF2E>>(num_poly);
        Y_poly_eval_shares[rep]=std::vector<std::vector<std::vector<field::GF2E>>>(num_poly);
        tmp=(*indice);
        for(size_t l=0;l<num_poly;l++){
            Y_poly[rep][l]=X1_poly[rep][l]*X2_poly[rep][l];
            Y_poly_eval_shares[rep][l]=Share_poly_eval(evals_d(Y_poly[rep][l],deg_poly),random_tapes,instance.num_MPC_parties,rep,indice,lambda_size,&(*delta_Y_poly)[rep]);
        }        
        if(rep!=instance.num_rounds-1){
            *indice=tmp;
        }
    }
    
    
    (*Y_poly_shares)=std::vector<std::vector<std::vector<std::vector<field::GF2E>>>>(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        (*Y_poly_shares)[rep]=std::vector<std::vector<std::vector<field::GF2E>>>(instance.num_MPC_parties);
    for(size_t i=0;i<instance.num_MPC_parties;i++){
        (*Y_poly_shares)[rep][i]=std::vector<std::vector<field::GF2E>>(num_poly);
        for(size_t l=0;l<num_poly;l++){
            (*Y_poly_shares)[rep][i][l]=std::vector<field::GF2E>(2*deg_poly-1);
            for(size_t k=0;k<deg_poly;k++){
                (*Y_poly_shares)[rep][i][l][k]=Y_shares[rep][i][l*deg_poly+k];
            }
            for(size_t k2=deg_poly;k2<2*deg_poly-1;k2++){
                (*Y_poly_shares)[rep][i][l][k2]=Y_poly_eval_shares[rep][l][i][k2-deg_poly];
            }
            (*Y_poly_shares)[rep][i][l] = field::interpolate_with_precomputation(
                precomputation_for_zero_to_2d, (*Y_poly_shares)[rep][i][l]);
        }
    }
    }


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
        for(size_t k=0;k<instance.d;k++){
            hash_update_GF2E(&ctx, lambda_size, (*delta_Y)[repetition][k]);
        }
        for(size_t j=0;j<(deg_poly-1)*num_poly;j++){
            hash_update_GF2E(&ctx, lambda_size, (*delta_Y_poly)[repetition][j]);
        }
    }
    hash_final(&ctx);

    std::vector<uint8_t> commitment(instance.digest_size);
    hash_squeeze(&ctx, commitment.data(), commitment.size());
    return commitment;
}




std::vector<std::vector<field::GF2E>>
phase_1_expand(const banquet_instance_t &instance,
               const std::vector<uint8_t> &h_1) {
    size_t lambda_size=instance.RMFE_fieldsize>>3;
    if((instance.RMFE_fieldsize&7)!=0){
        lambda_size++;
    } 

    hash_context ctx;
    hash_init(&ctx, instance.digest_size);
    hash_update(&ctx, h_1.data(), h_1.size());
    hash_final(&ctx);

    std::vector<std::vector<field::GF2E>> R_el(instance.num_rounds);
    std::vector<uint8_t> lambda_sized_buffer(lambda_size);
    for (size_t e = 0; e < instance.num_rounds; e++) {
        R_el[e]=std::vector<field::GF2E>(instance.num_poly);
            for (size_t l = 0; l < instance.num_poly; l++) {
            hash_squeeze(&ctx, lambda_sized_buffer.data(),
                        lambda_size);
            R_el[e][l].from_bytes(lambda_sized_buffer.data());
            // on réduit epsilon sinon la multiplication se fait mal :
            if((instance.RMFE_fieldsize&7)!=0){
                R_el[e][l]*=field::GF2E(1);
            }
        }
    }
    return R_el;
}



/*------------------ PHASE 2 ------------------*/


std::vector<uint8_t> phase_2_masquage(const std::vector<uint8_t> &h_1, const banquet_salt_t &salt, RandomTapes random_tapes, size_t *indice, 
                                        std::vector<std::vector<std::vector<field::GF2E>>> X2_poly_eval, std::vector<std::vector<std::vector<field::GF2E>>> *a_el_shares, 
                                        std::vector<std::vector<field::GF2E>> *c_e_shares, std::vector<field::GF2E> *delta_c,  banquet_instance_t &instance, 
                                        size_t byte_size){
    gsl::span<uint8_t> tape;
    std::vector<uint8_t> tmp(byte_size);

    (*delta_c)=std::vector<field::GF2E>(instance.num_rounds);
    std::vector<field::GF2E> c_voulu(instance.num_rounds);


    

    std::vector<std::vector<field::GF2E>> a_el(instance.num_rounds);
    (*a_el_shares)=std::vector<std::vector<std::vector<field::GF2E>>>(instance.num_rounds);
    (*c_e_shares)=std::vector<std::vector<field::GF2E>>(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        a_el[rep]=std::vector<field::GF2E>(instance.num_poly);
        (*a_el_shares)[rep]=std::vector<std::vector<field::GF2E>>(instance.num_MPC_parties);
        (*c_e_shares)[rep]=std::vector<field::GF2E>(instance.num_MPC_parties);
        for(size_t i=0;i<instance.num_MPC_parties;i++){
            tape=random_tapes.get_bytes(rep,i,*indice,(instance.num_poly+1)*byte_size);
            (*a_el_shares)[rep][i]=std::vector<field::GF2E>(instance.num_poly);
            for(size_t l=0;l<instance.num_poly;l++){
                for(size_t b=0;b<byte_size;b++){
                    tmp[b]=tape[l*byte_size+b];
                }
                (*a_el_shares)[rep][i][l].from_bytes(tmp.data());
                a_el[rep][l]+=(*a_el_shares)[rep][i][l];
            }
            for(size_t b=0;b<byte_size;b++){
                tmp[b]=tape[instance.num_poly*byte_size+b];
            }
            (*c_e_shares)[rep][i].from_bytes(tmp.data());
            (*delta_c)[rep]+=(*c_e_shares)[rep][i];
        }
    }


    field::GF2E tmp_X2;
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        for(size_t l=0;l<instance.num_poly;l++){
            tmp_X2=field::lift_uint8_t(0);
            for(size_t i=0;i<instance.num_MPC_parties;i++){
                tmp_X2+=X2_poly_eval[rep][i][l];
            }
            c_voulu[rep]+=a_el[rep][l]*tmp_X2;
        }
        (*delta_c)[rep]+=c_voulu[rep];
        (*c_e_shares)[rep][0]+=(*delta_c)[rep];
    }
    
    *indice+=(instance.num_poly+1)*byte_size;



    hash_context ctx;
    hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_2);
    hash_update(&ctx, salt.data(), salt.size());
    hash_update(&ctx, h_1.data(), h_1.size());

    for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
        hash_update_GF2E(&ctx, byte_size, (*delta_c)[repetition]);
    }
    hash_final(&ctx);

    std::vector<uint8_t> commitment(instance.digest_size);
    hash_squeeze(&ctx, commitment.data(), commitment.size());
    return commitment;
}





std::vector<std::vector<field::GF2E>>
phase_2_expand(const banquet_instance_t &instance,
               const std::vector<uint8_t> &h_2) {
    size_t lambda_size=instance.RMFE_fieldsize>>3;
    if((instance.RMFE_fieldsize&7)!=0){
        lambda_size++;
    } 

    hash_context ctx;
    hash_init(&ctx, instance.digest_size);
    hash_update(&ctx, h_2.data(), h_2.size());
    hash_final(&ctx);

    std::vector<uint8_t> lambda_sized_buffer(lambda_size);
    std::vector<std::vector<field::GF2E>> epsilon_el(instance.num_rounds);

    for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
        epsilon_el[repetition]=std::vector<field::GF2E>(instance.num_poly);
        for(size_t l=0;l<instance.num_poly;l++){
            hash_squeeze(&ctx, lambda_sized_buffer.data(),
                        lambda_size);
            epsilon_el[repetition][l].from_bytes(lambda_sized_buffer.data());
        }
    }
    return epsilon_el;
}






/*------------------ PHASE 3 ------------------*/

std::vector<uint8_t> phase_3_sacrifycing(const std::vector<uint8_t> &h_2, const banquet_salt_t &salt, std::vector<std::vector<std::vector<field::GF2E>>> X1_poly_eval,
                                        std::vector<std::vector<std::vector<field::GF2E>>> X2_poly_eval, std::vector<std::vector<std::vector<field::GF2E>>> Y_poly_eval,
                                        std::vector<std::vector<std::vector<field::GF2E>>> a_el_shares, 
                                        std::vector<std::vector<field::GF2E>> c_e_shares, 
                                        banquet_instance_t &instance, size_t byte_size,
                                        std::vector<std::vector<field::GF2E>> *alpha_el, std::vector<std::vector<field::GF2E>> epsilon_el){

    (*alpha_el)=std::vector<std::vector<field::GF2E>>(instance.num_rounds);
    std::vector<std::vector<std::vector<field::GF2E>>> alpha_el_shares(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        alpha_el_shares[rep]=std::vector<std::vector<field::GF2E>>(instance.num_MPC_parties);
        (*alpha_el)[rep]=std::vector<field::GF2E>(instance.num_poly);
        for(size_t party=0;party<instance.num_MPC_parties;party++){
            alpha_el_shares[rep][party]=std::vector<field::GF2E>(instance.num_poly);
            for(size_t l=0;l<instance.num_poly;l++){
                alpha_el_shares[rep][party][l]=epsilon_el[rep][l]*X1_poly_eval[rep][party][l]+a_el_shares[rep][party][l];
                (*alpha_el)[rep][l]+=alpha_el_shares[rep][party][l];
            }
        }
    }


    std::vector<std::vector<field::GF2E>> v_e_shares(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        v_e_shares[rep]=std::vector<field::GF2E>(instance.num_MPC_parties);
        for(size_t party=0;party<instance.num_MPC_parties;party++){
            v_e_shares[rep][party]=c_e_shares[rep][party];
            for(size_t l=0;l<instance.num_poly;l++){
                v_e_shares[rep][party]+=(*alpha_el)[rep][l]*X2_poly_eval[rep][party][l]+epsilon_el[rep][l]*Y_poly_eval[rep][party][l];
            }
        }
    }

    // Affichage des alpha_el et v_e en parties :
    // afficher_alpha_v_part(alpha_el_shares[round_idx],v_e_shares[round_idx],instance);


    hash_context ctx;
    hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_3);
    hash_update(&ctx, salt.data(), salt.size());
    hash_update(&ctx, h_2.data(), h_2.size());
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        for(size_t party=0;party<instance.num_MPC_parties;party++){
            for(size_t l=0;l<instance.num_poly;l++){
                hash_update_GF2E(&ctx,byte_size,alpha_el_shares[rep][party][l]);
            }
            hash_update_GF2E(&ctx,byte_size,v_e_shares[rep][party]);
        }
    }
    hash_final(&ctx);

    std::vector<uint8_t> commitment(instance.digest_size);
    hash_squeeze(&ctx, commitment.data(), commitment.size());
    return commitment;
}




std::vector<uint16_t> phase_3_expand(const banquet_instance_t &instance,
                                     const std::vector<uint8_t> &h_3) {
  assert(instance.num_MPC_parties < (1ULL << 16));
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, h_3.data(), h_3.size());
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
                                std::vector<std::vector<ascon_type_t>> ciphertext_out, size_t lambda_size){
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
        for(size_t k=0;k<instance.d;k++){
        hash_update_GF2E(&ctx, lambda_size, sign.proofs[repetition].delta_Y[k]);
        }
        for(size_t k=0;k<sign.proofs[repetition].delta_Y_poly.size();k++){
        hash_update_GF2E(&ctx, lambda_size, sign.proofs[repetition].delta_Y_poly[k]);
        }
    }
    hash_final(&ctx);

    std::vector<uint8_t> h_1(instance.digest_size);
    hash_squeeze(&ctx, h_1.data(), h_1.size());
    return h_1;
}



std::vector<uint8_t> h_2_commit(banquet_instance_t &instance, banquet_signature_t sign, size_t byte_size){
    hash_context ctx;
    hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_2);
    hash_update(&ctx, sign.salt.data(), sign.salt.size());
    hash_update(&ctx, sign.h_1.data(), sign.h_1.size());
    for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
        hash_update_GF2E(&ctx, byte_size, sign.proofs[repetition].delta_c);
    }
    hash_final(&ctx);

    std::vector<uint8_t> h_2(instance.digest_size);
    hash_squeeze(&ctx, h_2.data(), h_2.size());
    return h_2;
}



std::vector<uint8_t> h_3_commit(banquet_instance_t &instance, banquet_signature_t sign, std::vector<uint8_t> h_2, size_t lambda_size,
                                std::vector<std::vector<std::vector<field::GF2E>>> alpha_el_shares, std::vector<std::vector<field::GF2E>> v_e_shares){
    hash_context ctx;
    hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_3);
    hash_update(&ctx, sign.salt.data(), sign.salt.size());
    hash_update(&ctx, h_2.data(), h_2.size());
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        for(size_t party=0;party<instance.num_MPC_parties;party++){
            for(size_t l=0;l<instance.num_poly;l++){
                hash_update_GF2E(&ctx,lambda_size,alpha_el_shares[rep][party][l]);
            }
            hash_update_GF2E(&ctx,lambda_size,v_e_shares[rep][party]);
        }
    }
    hash_final(&ctx);

    std::vector<uint8_t> commitment(instance.digest_size);
    hash_squeeze(&ctx, commitment.data(), commitment.size());
    return commitment;
}











// // /*-----------------------------PHASES POUR VÉRIFICATION-----------------------------*/


std::vector<std::vector<field::GF2E>> Share_poly_eval_verif(RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, 
                                                      size_t byte_size, size_t d,
                                                      std::vector<field::GF2E> delta_Y_poly, size_t *delta_ind){
    std::vector<std::vector<field::GF2E>> res(N);
    gsl::span<uint8_t> tape;
    uint64_t tmp;
    for(int i=0;i<N;i++){
        if(i<N_ind){
          tape=random_tapes.get_bytes(rep,i,*indice,byte_size*(d-1));
        }
        else{
          tape=random_tapes.get_bytes(rep,i+1,*indice,byte_size*(d-1));
        }
        res[i]=std::vector<field::GF2E>(d-1);
        tmp=0;
        for(int j=0;j<d-1;j++){
            for(int z=0;z<byte_size;z++){
              tmp<<=8;
              tmp^=tape[j*byte_size+z];
            }
            res[i][j]=field::lift_uint8_t(tmp);
        }
    }

    if(N_ind!=0){
      for(int j=0;j<d-1;j++){
        res[0][j]+=delta_Y_poly[*delta_ind];
        (*delta_ind)++;
      }
    }
    
    *indice+=byte_size*(d-1);
    return res;
}






void calcul_MPC_eval_poly(banquet_signature_t sign, banquet_instance_t &instance, std::vector<uint8_t> pk, 
                        RandomTapes random_tapes, std::vector<size_t> N_ind, size_t *indice,
                        std::vector<std::vector<ascon_type_t>> *ciphertext_out, std::vector<std::vector<field::GF2E>> R_el,
                        std::vector<std::vector<std::vector<field::GF2E>>> *X1_poly_eval,
                        std::vector<std::vector<std::vector<field::GF2E>>> *X2_poly_eval,
                        std::vector<std::vector<std::vector<field::GF2E>>> *Y_poly_eval){
    size_t lambda_size=instance.RMFE_fieldsize>>3;
    if((instance.RMFE_fieldsize&7)!=0){
        lambda_size++;
    } 
    
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
    size_t N=instance.num_MPC_parties-1;
    (*ciphertext_out)=std::vector<std::vector<ascon_type_t>>(instance.num_rounds);
    std::vector<std::vector<std::vector<field::GF2E>>> X1_shares(instance.num_rounds);
    std::vector<std::vector<std::vector<field::GF2E>>> X2_shares(instance.num_rounds);
    std::vector<std::vector<std::vector<field::GF2E>>> Y_shares(instance.num_rounds);
    std::vector<field::GF2E> delta_Y;
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        (*ciphertext_out)[rep]=std::vector<ascon_type_t>(instance.num_MPC_parties);
        X1_shares[rep]=std::vector<std::vector<field::GF2E>>(N);
        X2_shares[rep]=std::vector<std::vector<field::GF2E>>(N);
        Y_shares[rep]=std::vector<std::vector<field::GF2E>>(N);
        tmp.first=*indice;
        delta_Y=sign.proofs[rep].delta_Y;
        ascon_128_groupe_verif(&(key_in[rep]),&plaintext_in,&nonce_in,&((*ciphertext_out)[rep]),&(X1_shares[rep]),&(X2_shares[rep]),&(Y_shares[rep]),
                                    delta_Y,random_tapes,N_ind[rep],rep,indice,instance.RMFE_vecsize,lambda_size);
        (*ciphertext_out)[rep][N_ind[rep]]=Rec2({Rec2((*ciphertext_out)[rep]),conversion_ascon(pk,4,2)});
        if(rep!=instance.num_rounds-1){
            *indice=tmp.first;
        }
    }


    // Affichage des parties d'entrées et sortie des portes AND :
    // afficher_in_out_part(X1_shares[round_idx],X2_shares[round_idx],Y_shares[round_idx],instance.num_MPC_parties-1,instance.d);


    // Affichage des partitions de la sortie :
    // afficher_ciphertext_part((*ciphertext_out)[round_idx],instance.num_MPC_parties,0);


    size_t deg_poly=instance.deg_poly;
    size_t num_poly=instance.num_poly;

    std::vector<field::GF2E> x_values_for_interpolation_zero_to_d =
      field::get_first_n_field_elements(deg_poly);
    std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_d =
      field::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_d);
    std::vector<field::GF2E> x_values_for_interpolation_zero_to_2d =
      field::get_first_n_field_elements(2 * deg_poly - 1);
    std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_2d =
      field::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_2d);

    std::vector<field::GF2E> X1_poly_shares_tmp(deg_poly);
    std::vector<field::GF2E> X2_poly_shares_tmp(deg_poly);
    (*X1_poly_eval)=std::vector<std::vector<std::vector<field::GF2E>>>(instance.num_rounds);
    (*X2_poly_eval)=std::vector<std::vector<std::vector<field::GF2E>>>(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        (*X1_poly_eval)[rep]=std::vector<std::vector<field::GF2E>>(N);
        (*X2_poly_eval)[rep]=std::vector<std::vector<field::GF2E>>(N);
        for(size_t i=0;i<N;i++){
        (*X1_poly_eval)[rep][i]=std::vector<field::GF2E>(instance.num_poly);
        (*X2_poly_eval)[rep][i]=std::vector<field::GF2E>(instance.num_poly);
        for(size_t l=0;l<num_poly;l++){
        for (size_t k = 0; k < deg_poly; k++) {
            X1_poly_shares_tmp[k] = X1_shares[rep][i][l*deg_poly+k];
            X2_poly_shares_tmp[k] = X2_shares[rep][i][l*deg_poly+k];
        }
        X1_poly_shares_tmp = field::interpolate_with_precomputation(
            precomputation_for_zero_to_d, X1_poly_shares_tmp);
        X2_poly_shares_tmp = field::interpolate_with_precomputation(
            precomputation_for_zero_to_d, X2_poly_shares_tmp);
        (*X1_poly_eval)[rep][i][l]=field::eval(X1_poly_shares_tmp,R_el[rep][l]);
        (*X2_poly_eval)[rep][i][l]=field::eval(X2_poly_shares_tmp,R_el[rep][l]);
        }
    }
    }

    delta_ind=0;
    std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> Y_poly_eval_shares(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        tmp.first=*indice;
        tmp.second=delta_ind;
        Y_poly_eval_shares[rep]=std::vector<std::vector<std::vector<field::GF2E>>>(num_poly);
        for(size_t l=0;l<num_poly;l++){
            Y_poly_eval_shares[rep][l]=Share_poly_eval_verif(random_tapes,instance.num_MPC_parties-1,N_ind[rep],
                                    rep,indice,lambda_size,deg_poly,sign.proofs[rep].delta_Y_poly,&delta_ind);
        }
        if(rep!=instance.num_rounds-1){
            *indice=tmp.first;
            delta_ind=tmp.second;
        }
    }


    std::vector<field::GF2E> Y_poly_shares_tmp(2*deg_poly-1);
    (*Y_poly_eval)=std::vector<std::vector<std::vector<field::GF2E>>>(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        (*Y_poly_eval)[rep]=std::vector<std::vector<field::GF2E>>(N);
        for(size_t i=0;i<N;i++){
        (*Y_poly_eval)[rep][i]=std::vector<field::GF2E>(num_poly);
        for(size_t l=0;l<num_poly;l++){
            for(size_t k=0;k<deg_poly;k++){
                Y_poly_shares_tmp[k]=Y_shares[rep][i][l*deg_poly+k];
            }
            for(size_t k2=deg_poly;k2<2*deg_poly-1;k2++){
                Y_poly_shares_tmp[k2]=Y_poly_eval_shares[rep][l][i][k2-deg_poly];
            }
        Y_poly_shares_tmp = field::interpolate_with_precomputation(
            precomputation_for_zero_to_2d, Y_poly_shares_tmp);
        (*Y_poly_eval)[rep][i][l] = field::eval(Y_poly_shares_tmp,R_el[rep][l]);
        }
    }
    }
}


void creation_a_c(RandomTapes random_tapes, size_t *indice, std::vector<std::vector<std::vector<field::GF2E>>> *a_el_shares, 
                                        std::vector<std::vector<field::GF2E>> *c_e_shares, banquet_signature_t sign,
                                        banquet_instance_t &instance, std::vector<size_t> N_ind,size_t byte_size, size_t r){
    gsl::span<uint8_t> tape;
    std::vector<uint8_t> tmp(byte_size);

   (*a_el_shares)=std::vector<std::vector<std::vector<field::GF2E>>>(instance.num_rounds);
    (*c_e_shares)=std::vector<std::vector<field::GF2E>>(instance.num_rounds);

    for(size_t rep=0;rep<instance.num_rounds;rep++){
        (*a_el_shares)[rep]=std::vector<std::vector<field::GF2E>>(instance.num_MPC_parties-1);
        (*c_e_shares)[rep]=std::vector<field::GF2E>(instance.num_MPC_parties-1);
        for(size_t i=0;i<instance.num_MPC_parties-1;i++){
            if(i<N_ind[rep]){
                tape=random_tapes.get_bytes(rep,i,*indice,(r+1)*byte_size);
            }
            else{
                tape=random_tapes.get_bytes(rep,i+1,*indice,(r+1)*byte_size);
            }
            (*a_el_shares)[rep][i]=std::vector<field::GF2E>(r);
            for(size_t l=0;l<r;l++){
                for(size_t b=0;b<byte_size;b++){
                    tmp[b]=tape[l*byte_size+b];
                }
                (*a_el_shares)[rep][i][l].from_bytes(tmp.data());
            }
            for(size_t b=0;b<byte_size;b++){
                tmp[b]=tape[r*byte_size+b];
            }
            (*c_e_shares)[rep][i].from_bytes(tmp.data());
        }
    }

    for(size_t rep=0;rep<instance.num_rounds;rep++){
        if(N_ind[rep]!=0){
            (*c_e_shares)[rep][0]+=sign.proofs[rep].delta_c;
        }
    }
    
    *indice+=(r+1)*byte_size;
}


void calcul_alpha_v_shares(std::vector<std::vector<std::vector<field::GF2E>>> X1_poly_eval,
                                        std::vector<std::vector<std::vector<field::GF2E>>> X2_poly_eval, std::vector<std::vector<std::vector<field::GF2E>>> Y_poly_eval,
                                        std::vector<std::vector<std::vector<field::GF2E>>> a_el_shares, 
                                        std::vector<std::vector<field::GF2E>> c_e_shares, 
                                        banquet_instance_t &instance, size_t byte_size, size_t r,
                                        std::vector<std::vector<std::vector<field::GF2E>>> *alpha_el_shares, std::vector<std::vector<field::GF2E>> *v_e_shares,
                                        std::vector<std::vector<field::GF2E>> epsilon_el,
                                        banquet_signature_t sign, std::vector<size_t> N_ind){

    int partynew;
    (*alpha_el_shares)=std::vector<std::vector<std::vector<field::GF2E>>>(instance.num_rounds);
    (*v_e_shares)=std::vector<std::vector<field::GF2E>>(instance.num_rounds);
    for(size_t repetition=0;repetition<instance.num_rounds;repetition++){
        (*alpha_el_shares)[repetition]=std::vector<std::vector<field::GF2E>>(instance.num_MPC_parties);
        (*v_e_shares)[repetition]=std::vector<field::GF2E>(instance.num_MPC_parties);
        (*alpha_el_shares)[repetition][N_ind[repetition]]=sign.proofs[repetition].alpha_l;
        for(size_t party=0;party<instance.num_MPC_parties-1;party++){
            if(party<N_ind[repetition]){
                partynew=party;
            }
            else{
                partynew=party+1;
            }
            (*alpha_el_shares)[repetition][partynew]=std::vector<field::GF2E>(r);
            (*v_e_shares)[repetition][partynew]=c_e_shares[repetition][party];
            for(size_t l=0;l<r;l++){
                (*alpha_el_shares)[repetition][partynew][l]=epsilon_el[repetition][l]*X1_poly_eval[repetition][party][l]+a_el_shares[repetition][party][l];
                (*v_e_shares)[repetition][partynew]+=sign.proofs[repetition].alpha_l[l]*X2_poly_eval[repetition][party][l]+epsilon_el[repetition][l]*Y_poly_eval[repetition][party][l];
            }
            (*alpha_el_shares)[repetition][N_ind[repetition]]+=(*alpha_el_shares)[repetition][partynew];
            (*v_e_shares)[repetition][N_ind[repetition]]+=(*v_e_shares)[repetition][partynew];
        }
    }

    // Affichage des alpha_el et v_e en parties :
    // afficher_alpha_v_part((*alpha_el_shares)[round_idx],(*v_e_shares)[round_idx],instance);
}

