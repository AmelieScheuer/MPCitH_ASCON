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



#include "fichiers_complementaires/field2.h"
#include "fichiers_complementaires/gen_alea.h"
#include "math.h"
#include <chrono>



std::vector<std::vector<field::GF2E>> Share_poly_eval(std::vector<field::GF2E> poly_evals,
                                                      RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, size_t byte_size, size_t lambda,
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
            res[i][j]=field::lift_uint8_t(tmp%(1<<lambda));
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

void Rec_polys(std::vector<std::vector<std::vector<field::GF2E>>> polys_shares, std::vector<std::vector<field::GF2E>> *polys){
  size_t N=polys_shares.size();
  size_t r=polys_shares[0].size();
  size_t d=polys_shares[0][0].size();
  (*polys)=std::vector<std::vector<field::GF2E>>(r);
  for(size_t l=0;l<r;l++){
    (*polys)[l]=std::vector<field::GF2E>(d);
    for(size_t i=0;i<N;i++){
      (*polys)[l]+=polys_shares[i][l];
    }
  }
}




std::vector<field::GF2E> evals_d(std::vector<field::GF2E> Y_poly, size_t d){
  std::vector<field::GF2E> poly_evals(d-1);
  for(int i=0;i<d-1;i++){
    poly_evals[i]=field::eval(Y_poly,field::lift_uint8_t(d+i));
  }
  return poly_evals;
}


void poly_shares(std::vector<std::vector<uint8_t>> x1_shares, std::vector<std::vector<uint8_t>> x2_shares,
                std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_d,
                size_t d, size_t r, size_t N, std::vector<std::vector<std::vector<field::GF2E>>> *X1_poly_shares,
                std::vector<std::vector<std::vector<field::GF2E>>> *X2_poly_shares){
  int quotient=0;
  int reste=0;
  (*X1_poly_shares)=std::vector<std::vector<std::vector<field::GF2E>>>(N);
  (*X2_poly_shares)=std::vector<std::vector<std::vector<field::GF2E>>>(N);
  for(size_t i=0;i<N;i++){
    (*X1_poly_shares)[i]=std::vector<std::vector<field::GF2E>>(r);
    (*X2_poly_shares)[i]=std::vector<std::vector<field::GF2E>>(r);
    reste=0;
    quotient=0;
    for(size_t l=0;l<r;l++){
      (*X1_poly_shares)[i][l]=std::vector<field::GF2E>(d);
      (*X2_poly_shares)[i][l]=std::vector<field::GF2E>(d);
      for (size_t k = 0; k < d; k++) {
        if(quotient<720){
          (*X1_poly_shares)[i][l][k] = field::lift_uint8_t((x1_shares[i][quotient]>>reste) & 1);
          (*X2_poly_shares)[i][l][k] = field::lift_uint8_t((x2_shares[i][quotient]>>reste) & 1);
        }
        reste++;
        if(reste==8){
          reste=0;
          quotient++;
        }
      }
      (*X1_poly_shares)[i][l] = field::interpolate_with_precomputation(
          precomputation_for_zero_to_d, (*X1_poly_shares)[i][l]);
      (*X2_poly_shares)[i][l] = field::interpolate_with_precomputation(
          precomputation_for_zero_to_d, (*X2_poly_shares)[i][l]);
    }
  }
}




void interpolation_poly(std::vector<std::vector<uint8_t>> x1_shares, std::vector<std::vector<uint8_t>> x2_shares, 
                        std::vector<std::vector<uint8_t>> y_shares, size_t d, size_t r, 
                        std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_d, std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_2d,
                        RandomTapes random_tapes, size_t N, size_t rep, size_t *indice, size_t lambda, size_t byte_size,
                        std::vector<field::GF2E> *delta_Y_poly,
                        std::vector<std::vector<std::vector<field::GF2E>>> *X1_poly_shares, std::vector<std::vector<std::vector<field::GF2E>>> *X2_poly_shares,
                        std::vector<std::vector<std::vector<field::GF2E>>> *Y_poly_shares){
  
  // 1ère étape : on interpole les parties des polynômes dans X1_poly_shares et X2_poly_shares
  auto start=std::chrono::high_resolution_clock::now();
  poly_shares(x1_shares,x2_shares,precomputation_for_zero_to_d,d,r,N,X1_poly_shares,X2_poly_shares);
  auto stop=std::chrono::high_resolution_clock::now();
  auto duration=std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start);
  // if(rep==0){
  //   printf("\nDurée du commit de l'étape 1 : %.9f seconds.\n",duration.count()*1e-9);
  // }


  // 2ème étape : on recompose les parties des polynômes pour avoir X1_poly et X2_poly
  start=std::chrono::high_resolution_clock::now();
  std::vector<std::vector<field::GF2E>> X1_poly,X2_poly;
  Rec_polys(*X1_poly_shares,&X1_poly);
  Rec_polys(*X2_poly_shares,&X2_poly);
  stop=std::chrono::high_resolution_clock::now();
  duration=std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start);
  // if(rep==0){
  //   printf("\nDurée du commit de l'étape 2 : %.9f seconds.\n",duration.count()*1e-9);
  // }


  // Affichage des polynômes X1,X2 reconstitués :
    // for(size_t l=0;l<2;l++){
    //   printf("X1 numero %d = ",l); afficher_poly(X1_poly[l]); printf("\n\n");
    //   for(int j=0;j<d;j++){
    //     eval(X1_poly[l],field::lift_uint8_t(j)).afficher(); printf(" ");
    //   }
    //   printf("\n\n");
    //   printf("X2 numero %d = ",l); afficher_poly(X2_poly[l]); printf("\n\n");
    //   for(int j=0;j<d;j++){
    //     eval(X2_poly[l],field::lift_uint8_t(j)).afficher(); printf(" ");
    //   }
    //   printf("\n\n");
    // }


  // 3ème étape : on calcule Y_poly=X1_poly*X2_poly 
  // puis on évalue Y_poly sur d-1 points que l'on partage ensuite avec la fonction Share
  start=std::chrono::high_resolution_clock::now();
  std::vector<std::vector<field::GF2E>> Y_poly(r);
  std::vector<std::vector<std::vector<field::GF2E>>> Y_poly_eval_shares(r);
  for(size_t k=0;k<r;k++){
    Y_poly[k]=X1_poly[k]*X2_poly[k];
    Y_poly_eval_shares[k]=Share_poly_eval(evals_d(Y_poly[k],d),random_tapes,N,rep,indice,byte_size,lambda,delta_Y_poly);
  }
  stop=std::chrono::high_resolution_clock::now();
  duration=std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start);
  // if(rep==0){
  //   printf("\nDurée du commit de l'étape 3 : %.9f seconds.\n",duration.count()*1e-9);
  // }


  // Affichage du polynôme produit Y :
    // for(size_t l=0;l<2;l++){
    //   printf("Y numero %d = ",l); afficher_poly(Y_poly[l]); printf("\n\n");
    //   for(int j=0;j<d;j++){
    //     eval(Y_poly[l],field::lift_uint8_t(j)).afficher(); printf(" ");
    //   }
    //   printf("\n\n");
    // }



  // 4ème étape : on interpole les parties des polynômes dans Y_poly_shares
  start=std::chrono::high_resolution_clock::now();
  *Y_poly_shares=std::vector<std::vector<std::vector<field::GF2E>>>(N);
  int quotient=0;
  int reste=0;
  for(size_t i=0;i<N;i++){
    (*Y_poly_shares)[i]=std::vector<std::vector<field::GF2E>>(r);
    reste=0;
    quotient=0;
    for(size_t l=0;l<r;l++){
      (*Y_poly_shares)[i][l]=std::vector<field::GF2E>(2*d-1);
      for(size_t k=0;k<d;k++){
        if(quotient<720){
          (*Y_poly_shares)[i][l][k]=field::lift_uint8_t((y_shares[i][quotient]>>reste) & 1);
        }
        reste++;
        if(reste==8){
          reste=0;
          quotient++;
        }
      }
      for(size_t k2=d;k2<2*d-1;k2++){
        (*Y_poly_shares)[i][l][k2]=Y_poly_eval_shares[l][i][k2-d];
      }
      (*Y_poly_shares)[i][l] = field::interpolate_with_precomputation(
          precomputation_for_zero_to_2d, (*Y_poly_shares)[i][l]);
    }
  }
  stop=std::chrono::high_resolution_clock::now();
  duration=std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start);
  // if(rep==0){
  //   printf("\nDurée du commit de l'étape 4 : %.9f seconds.\n",duration.count()*1e-9);
  // }
}

















std::vector<std::vector<field::GF2E>> Share_poly_eval_verif(RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, 
                                                      size_t byte_size, size_t lambda, size_t d,
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
            res[i][j]=field::lift_uint8_t(tmp%(1<<lambda));
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


void interpolation_poly_verif(std::vector<std::vector<uint8_t>> x1_shares, std::vector<std::vector<uint8_t>> x2_shares, 
                        std::vector<std::vector<uint8_t>> y_shares, size_t d, size_t r, 
                        std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_d, std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_2d,
                        RandomTapes random_tapes, size_t N, size_t N_ind, size_t rep, size_t *indice, size_t lambda, size_t byte_size,
                        std::vector<field::GF2E> delta_Y_poly, size_t *delta_ind,
                        std::vector<std::vector<std::vector<field::GF2E>>> *X1_poly_shares, std::vector<std::vector<std::vector<field::GF2E>>> *X2_poly_shares,
                        std::vector<std::vector<std::vector<field::GF2E>>> *Y_poly_shares){

  poly_shares(x1_shares,x2_shares,precomputation_for_zero_to_d,d,r,N,X1_poly_shares,X2_poly_shares);

  std::vector<std::vector<std::vector<field::GF2E>>> Y_poly_eval_shares(r);
  for(size_t k=0;k<r;k++){
    Y_poly_eval_shares[k]=Share_poly_eval_verif(random_tapes,N,N_ind,rep,indice,byte_size,lambda,d,delta_Y_poly,delta_ind);
  }

  *Y_poly_shares=std::vector<std::vector<std::vector<field::GF2E>>>(N);
  int quotient=0;
  int reste=0;
  for(size_t i=0;i<N;i++){
    (*Y_poly_shares)[i]=std::vector<std::vector<field::GF2E>>(r);
    reste=0;
    quotient=0;
    for(size_t l=0;l<r;l++){
      (*Y_poly_shares)[i][l]=std::vector<field::GF2E>(2*d-1);
      for(size_t k=0;k<d;k++){
        if(quotient<720){
          (*Y_poly_shares)[i][l][k]=field::lift_uint8_t((y_shares[i][quotient]>>reste) & 1);
        }
        reste++;
        if(reste==8){
          reste=0;
          quotient++;
        }
      }
      for(size_t k2=d;k2<2*d-1;k2++){
        (*Y_poly_shares)[i][l][k2]=Y_poly_eval_shares[l][i][k2-d];
      }
      (*Y_poly_shares)[i][l] = field::interpolate_with_precomputation(
          precomputation_for_zero_to_2d, (*Y_poly_shares)[i][l]);
    }
  }  
}

