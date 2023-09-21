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



#include "MPC_phases_ascon.h"
#include <chrono>




banquet_signature_t sign_helium_ascon(banquet_instance_t &instance, const banquet_keypair_t &keypair, const uint8_t *msg, size_t msg_len){
    size_t lambda=instance.lambda;
    size_t byte_size;
    if((lambda & 7) == 0){
        byte_size=lambda>>3;
    }
    else{
        byte_size=(lambda>>3)+1;
    }
    int d=instance.d;
    int r=(int)ceil(5760/(float)d);

    

   // étape 1 : générer les graines et le sel
    auto [salt, master_seeds] =
      generate_salt_and_seeds(instance, keypair, msg, msg_len);

    std::vector<SeedTree> seed_trees;
    seed_trees.reserve(instance.num_rounds);

    size_t random_tape_size = 16+720+byte_size*(d*r+1);
    RandomTapes random_tapes(instance.num_rounds, instance.num_MPC_parties,
                            random_tape_size);

    RepByteContainer party_seed_commitments(
        instance.num_rounds, instance.num_MPC_parties, instance.digest_size);

    for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
        // generate seed tree for the N parties
        seed_trees.emplace_back(master_seeds[repetition], instance.num_MPC_parties,
                                salt, repetition);

        //commit to each party's seed;
        size_t party = 0;
        for (; party < (instance.num_MPC_parties >> 2) << 2; party += 4) {
            commit_to_4_party_seeds(
                instance, seed_trees[repetition].get_leaf(party).value(),
                seed_trees[repetition].get_leaf(party + 1).value(),
                seed_trees[repetition].get_leaf(party + 2).value(),
                seed_trees[repetition].get_leaf(party + 3).value(), salt,
                repetition, party, party_seed_commitments.get(repetition, party),
                party_seed_commitments.get(repetition, party + 1),
                party_seed_commitments.get(repetition, party + 2),
                party_seed_commitments.get(repetition, party + 3));
        }
        for (; party < instance.num_MPC_parties; party++) {
            commit_to_party_seed(
                instance, seed_trees[repetition].get_leaf(party).value(), salt,
                repetition, party, party_seed_commitments.get(repetition, party));
        }
        

        // create random tape for each party
        for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        random_tapes.generate_tape(
            repetition, party, salt,
            seed_trees[repetition].get_leaf(party).value());
        }
    }

    // Affichage du sel et des graines :
    // afficher_salt(salt);
    // afficher_seeds(seed_trees[round_idx],instance.num_MPC_parties,instance.num_MPC_parties);


    // On fait la phase 1 en calculant le commit h_1 puis en l'étendant :
    size_t indice=0;
    std::vector<std::vector<uint8_t>> delta_key,delta_y;
    std::vector<std::vector<field::GF2E>> delta_Y_poly_eval;
    std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> X1_poly_shares,X2_poly_shares,Y_poly_shares;

    std::vector<uint8_t> h_1=phase_1_MPC_compute(msg,msg_len,salt,party_seed_commitments,instance,random_tapes,
                            keypair,&indice,&delta_key,&delta_y,&delta_Y_poly_eval,&X1_poly_shares,&X2_poly_shares,&Y_poly_shares);
    //afficher_hash(h_1,1);

    std::vector<std::vector<field::GF2E>> R_el=phase_1_expand(instance,h_1);

    // Affichage des points R_el :
    // afficher_R_el(R_el[round_idx],r);


    // évaluation des polynômes en les R_el :
    std::vector<std::vector<std::vector<field::GF2E>>> X1_poly_eval(instance.num_rounds),X2_poly_eval(instance.num_rounds),Y_poly_eval(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        X1_poly_eval[rep]=std::vector<std::vector<field::GF2E>>(instance.num_MPC_parties);
        X2_poly_eval[rep]=std::vector<std::vector<field::GF2E>>(instance.num_MPC_parties);
        Y_poly_eval[rep]=std::vector<std::vector<field::GF2E>>(instance.num_MPC_parties);
        for(size_t party=0;party<instance.num_MPC_parties;party++){
            X1_poly_eval[rep][party]=std::vector<field::GF2E>(r);
            X2_poly_eval[rep][party]=std::vector<field::GF2E>(r);
            Y_poly_eval[rep][party]=std::vector<field::GF2E>(r);
            for(size_t l=0;l<r;l++){
                X1_poly_eval[rep][party][l]=eval(X1_poly_shares[rep][party][l],R_el[rep][l]);
                X2_poly_eval[rep][party][l]=eval(X2_poly_shares[rep][party][l],R_el[rep][l]);
                Y_poly_eval[rep][party][l]=eval(Y_poly_shares[rep][party][l],R_el[rep][l]);
            }
        }
    }

    // Affichage des évaluations en R_el des polynômes :
    // afficher_eval_poly_part(X1_poly_eval[round_idx],X1_poly_eval[round_idx],
    //                         X1_poly_eval[round_idx],R_el[round_idx],instance.num_MPC_parties,number);

    // On fait la phase 2 en calculant le commit h_2 puis en l'étendant :
    std::vector<std::vector<std::vector<field::GF2E>>> a_el_shares;
    std::vector<std::vector<field::GF2E>> c_e_shares;
    std::vector<field::GF2E> delta_c;
    std::vector<uint8_t> h_2=phase_2_masquage(h_1,salt,random_tapes,&indice,X2_poly_eval,&a_el_shares,&c_e_shares,&delta_c,instance,byte_size,r);
    //afficher_hash(h_2,2);

    // Affichage des a_el et c_e en parties :
    // afficher_a_c_part(a_el_shares[round_idx],c_e_shares[round_idx],instance.num_MPC_parties,r);

    std::vector<std::vector<field::GF2E>> epsilon_el=phase_2_expand(instance,h_2);

    // On fait la phase 3 en calculant le commit h_3 puis en l'étendant :
    std::vector<std::vector<field::GF2E>> alpha_el;
    std::vector<uint8_t> h_3=phase_3_sacrifycing(h_2,salt,X1_poly_eval,X2_poly_eval,Y_poly_eval,a_el_shares,c_e_shares,instance,byte_size,r,&alpha_el,epsilon_el);
    //afficher_hash(h_3,3);

    std::vector<uint16_t> missing_parties=phase_3_expand(instance,h_3);
    // afficher_miss_part(missing_parties[round_idx]);

    // Construction de la signature :
    std::vector<reveal_list_t> revealed_seeds;
    for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
        revealed_seeds.push_back(seed_trees[repetition].reveal_all_but(missing_parties[repetition]));
    }

    

    std::vector<banquet_repetition_proof_t> proofs;
    for(size_t repetition=0;repetition<instance.num_rounds;repetition++){
        std::vector<uint8_t> commitment(instance.digest_size);
        auto missing_commitment =
            party_seed_commitments.get(repetition, missing_parties[repetition]);
        std::copy(std::begin(missing_commitment), std::end(missing_commitment),
              std::begin(commitment));


        banquet_repetition_proof_t proof{
        revealed_seeds[repetition],
        commitment,
        delta_key[repetition],
        delta_y[repetition],
        delta_Y_poly_eval[repetition],
        alpha_el[repetition],
        delta_c[repetition]};

        proofs.push_back(proof);
    }

    banquet_signature_t signature{salt, h_1, h_3, proofs};
    return signature;
}








bool verif(banquet_instance_t &instance, std::vector<uint8_t> pk, const uint8_t *msg, size_t msg_len, banquet_signature_t sign){
    size_t lambda=instance.lambda;
    size_t byte_size;
    if((lambda & 7) == 0){
        byte_size=lambda>>3;
    }
    else{
        byte_size=(lambda>>3)+1;
    }
    int d=instance.d;
    int r=(int)ceil(5760/(float)d);


    


    std::vector<uint8_t> h_2=h_2_commit(instance,sign,byte_size);

    // Affichage des h_i :
    // afficher_hash(sign.h_1,1);
    // afficher_hash(h_2,2);
    // afficher_hash(sign.h_3,3);

    std::vector<std::vector<field::GF2E>> R_el=phase_1_expand(instance,sign.h_1);
    std::vector<std::vector<field::GF2E>> epsilon_el=phase_2_expand(instance,h_2);
    std::vector<uint16_t> missing_parties=phase_3_expand(instance,sign.h_3);

    // Affichage des points R_el :
    // afficher_R_el(R_el[round_idx],r);


    std::vector<size_t> N_ind(instance.num_rounds);
    for(size_t rep=0;rep<instance.num_rounds;rep++){
        N_ind[rep]=sign.proofs[rep].reveallist.second;
    }

    std::vector<SeedTree> seed_trees;
    seed_trees.reserve(instance.num_rounds);

    size_t random_tape_size = 16+720+byte_size*(d*r+1);
    RandomTapes random_tapes(instance.num_rounds, instance.num_MPC_parties,
                            random_tape_size);
    RepByteContainer party_seed_commitments(
        instance.num_rounds, instance.num_MPC_parties, instance.digest_size);


    // rebuild SeedTrees
    for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    const banquet_repetition_proof_t &proof = sign.proofs[repetition];
    // regenerate generate seed tree for the N parties (except the missing
    // one)
    if (missing_parties[repetition] != proof.reveallist.second)
      throw std::runtime_error(
          "modified signature between deserialization and verify");
    seed_trees.push_back(SeedTree(proof.reveallist, instance.num_MPC_parties,
                                  sign.salt, repetition));
    // commit to each party's seed, fill up missing one with data from proof
    {
      std::vector<uint8_t> dummy(instance.seed_size);
      size_t party = 0;
      for (; party < (instance.num_MPC_parties >> 2) << 2; party += 4) {
        auto seed0 = seed_trees[repetition].get_leaf(party).value_or(dummy);
        auto seed1 = seed_trees[repetition].get_leaf(party + 1).value_or(dummy);
        auto seed2 = seed_trees[repetition].get_leaf(party + 2).value_or(dummy);
        auto seed3 = seed_trees[repetition].get_leaf(party + 3).value_or(dummy);
        commit_to_4_party_seeds(
            instance, seed0, seed1, seed2, seed3, sign.salt, repetition,
            party, party_seed_commitments.get(repetition, party),
            party_seed_commitments.get(repetition, party + 1),
            party_seed_commitments.get(repetition, party + 2),
            party_seed_commitments.get(repetition, party + 3));
      }
      for (; party < instance.num_MPC_parties; party++) {
        if (party != missing_parties[repetition]) {
          commit_to_party_seed(instance,
                               seed_trees[repetition].get_leaf(party).value(),
                               sign.salt, repetition, party,
                               party_seed_commitments.get(repetition, party));
        }
      }
    }
    auto com =
        party_seed_commitments.get(repetition, missing_parties[repetition]);
    std::copy(std::begin(proof.Commit), std::end(proof.Commit), std::begin(com));

    // create random tape for each party, dummy one for missing party
    {
      size_t party = 0;
      std::vector<uint8_t> dummy(instance.seed_size);
      for (; party < (instance.num_MPC_parties >> 2) << 2; party += 4) {
        random_tapes.generate_4_tapes(
            repetition, party, sign.salt,
            seed_trees[repetition].get_leaf(party).value_or(dummy),
            seed_trees[repetition].get_leaf(party + 1).value_or(dummy),
            seed_trees[repetition].get_leaf(party + 2).value_or(dummy),
            seed_trees[repetition].get_leaf(party + 3).value_or(dummy));
      }
      for (; party < instance.num_MPC_parties; party++) {
        random_tapes.generate_tape(
            repetition, party, sign.salt,
            seed_trees[repetition].get_leaf(party).value_or(dummy));
      }
    }
    }

    // Affichage du sel et des graines :
    // afficher_salt(sign.salt);
    // afficher_seeds(seed_trees[round_idx],instance.num_MPC_parties,N_ind[round_idx]);


    // Vérification de h_1 :


    // Calcul MPC :
    size_t indice=0;
    std::vector<std::vector<ascon_type_t>> ciphertext_out;
    std::vector<std::vector<std::vector<uint8_t>>> x1_shares,x2_shares,y_shares;
    calcul_MPC_ascon(sign,instance,pk,random_tapes,N_ind,&indice,&ciphertext_out,&x1_shares,&x2_shares,&y_shares);

    std::vector<uint8_t> h_1=h_1_commit(instance,pk,msg,msg_len,sign,party_seed_commitments,ciphertext_out,byte_size);
    // afficher_hash(h_1,1);


    if(h_1!=sign.h_1){
        printf("\nh_1 n'est pas acceptée.\n");
        printf("\nLa signature n'est pas acceptée.\n");
        return 0;
    }
    else{
        printf("\nh_1 est acceptée.\n");
    }


    // Vérification de h_3 :
    field::GF2E::init_extension_field(instance);

    std::vector<std::vector<std::vector<field::GF2E>>> X1_poly_eval,X2_poly_eval,Y_poly_eval;
    interpolation_evaluation_poly(sign,instance,random_tapes,N_ind,&indice,d,r,lambda,byte_size,x1_shares,x2_shares,y_shares,R_el,&X1_poly_eval,&X2_poly_eval,&Y_poly_eval);


    // Affichage des évaluations :
    // afficher_eval_part(X1_poly_eval[round_idx],X2_poly_eval[round_idx],
    //                     Y_poly_eval[round_idx],instance.num_MPC_parties-1,r);

    
    std::vector<std::vector<std::vector<field::GF2E>>> a_el_shares;
    std::vector<std::vector<field::GF2E>> c_e_shares;
    creation_a_c(random_tapes,&indice,&a_el_shares,&c_e_shares,sign,instance,N_ind,byte_size,r);

    // Affichage des a_el et c_e en parties :
    // afficher_a_c_part(a_el_shares[round_idx],c_e_shares[round_idx],instance.num_MPC_parties-1,r);

    std::vector<std::vector<std::vector<field::GF2E>>> alpha_el_shares;
    std::vector<std::vector<field::GF2E>> v_e_shares;
    calcul_alpha_v_shares(X1_poly_eval,X2_poly_eval,Y_poly_eval,a_el_shares,c_e_shares,instance,byte_size,r,&alpha_el_shares,&v_e_shares,epsilon_el,sign,N_ind);

    // Affichage des alpha_el et v_e en parties :
    // afficher_alpha_v_part(alpha_el_shares[round_idx],v_e_shares[round_idx],instance);


    std::vector<uint8_t> h_3=h_3_commit(instance,sign,h_2,byte_size,r,alpha_el_shares,v_e_shares);
    // afficher_hash(h_3,3);


    if(h_3!=sign.h_3){
        printf("\nLa signature n'est pas acceptée.\n");
        return 0;
    }
    else{
        printf("\nLa signature est acceptée.\n");
        return 1;
    }
}
