//
// Created by oyliy on 2019/11/25.
//

#ifndef WONAL_TRSM_H
#define WONAL_TRSM_H


#include "basewo.h"

//square matrix minus method
void matrix_minus(float **mata,float **matb,int ar,int ac,int br, int bc,int n){
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mata[ar+i][ac + j] -= matb[br+i][bc+j];
        }
    }
}


void direct_trsm(float **mata, float **matb, float **matx,
                 int ar, int ac, int xr, int xc,
                 int n){
    //A is lower triangular matrix
    //A is m*k , B is m*n X is k*n ,wo find x
    //square n
#pragma omp parallel for num_threads(core_num) default(none) shared(mata,matb,matx,n,ar,ac,xr,xc)
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            int temp = 0;
            for (int k = 0; k < i; ++k) {
                temp += mata[ar+i][ac+k]*matx[xr+k][xc+j];
            }
            matx[xr+i][xc+j] = (matb[ar+i][xc+j] - temp)/mata[ar+i][ac+i];
            temp = 0;
        }
    }
}

//AX = B
void woco_trsm(float **mata, float **matb, float **matx, float** temp,
               int ar, int ac,int xr, int xc, int n){
    if (3*n*n < CACHE_SIZE){
        direct_trsm(mata,matb,matx,ar,ac,xr,xc,n);
    } else{
        //partitiion to 4
        int p = n/2;
        woco_trsm(mata,matb,matx,temp,ar,ac,xr,xc,p);
        woco_trsm(mata,matb,matx,temp,ar,ac,xr,xc+p,p);
        if (flag == 1){vex_wo(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
        else if(flag == 2){omp_vec_mat_mult(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
//        else if(flag == 3){matrix_multi(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
        matrix_minus(matb,temp,ar+p,xc,ar+p,xc,p);
        woco_trsm(mata,matb,matx,temp,ar+p,ac+p,xr+p,xc,p);
        if (flag ==1){vex_wo(mata,matx,temp,ar+p,ac,xr,xc+p,p,p,p);}
        else if (flag == 2){omp_vec_mat_mult(mata,matx,temp,ar+p,ac,xr,xc+p,p,p,p);}
//        else if(flag == 3){matrix_multi(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
        matrix_minus(matb,matx,ar+p,xc+p,ar+p,xc+p,p);
        woco_trsm(mata,matb,matx,temp,ar+p,ac+p,xr+p,xc+p,p);
    }
}

void direct_trsm_u(float **mata, float **matb, float **matx,
                 int ar, int ac, int xr, int xc,
                 int n){
    //A is upper triangular matrix
    //A is m*k , B is m*n X is k*n ,wo find x
    //square n
#pragma omp parallel for num_threads(core_num) default(none) shared(mata,matb,matx,n,ar,ac,xr,xc)
    for (int j = 0; j < n; ++j) {    
        for (int i = n - 1; i >= 0; --i) {
        int temp = 0;
            for (int k = i; k < n ; ++k) {
                temp += mata[ar+i][ac+k]*matx[xr+k][xc+j];
            }
            matx[xr+i][xc+j] = (matb[ar+i][xc+j] - temp)/mata[ar+i][ac+i];
            temp = 0;
        }
    }
}

//need test
void woco_trsm_u(float **mata, float **matb, float **matx, float** temp,
               int ar, int ac,int xr, int xc, int n){
    if (3*n*n < CACHE_SIZE){
        direct_trsm_u(mata, matb, matx, ar, ac, xr, xc, n);
    } else{
        //partitiion to 4
        int p = n/2;
        woco_trsm_u(mata, matb, matx, temp, ar + p, ac + p, xr + p, xc + p , p);
        woco_trsm_u(mata, matb, matx, temp, ar + p, ac + p, xr + p, xc,p);
        if (flag == 1){vex_wo(mata,matx,temp,ar,ac+p,xr+p,xc+p,p,p,p);}
        else if(flag == 2){omp_vec_mat_mult(mata,matx,temp,ar,ac+p,xr+p,xc+p,p,p,p);}
//        else if(flag == 3){matrix_multi(mata,matx,temp,ar,ac+p,xr+p,xc+p,p,p,p);}
        matrix_minus(matb,temp,ar,xc+p,ar,xc+p,p);
        woco_trsm_u(mata,matb,matx,temp,ar,ac,xr,xc+p,p);
        if (flag ==1){vex_wo(mata,matx,temp,ar,ac+p,xr+p,xc,p,p,p);}
        else if (flag == 2){omp_vec_mat_mult(mata,matx,temp,ar,ac+p,xr+p,xc,p,p,p);}
//        else if(flag == 3){matrix_multi(mata,matx,temp,ar,ac+p,xr+p,xc,p,p,p);}
        matrix_minus(matb,matx,ar,xc,ar,xc,p);
        woco_trsm_u(mata,matb,matx,temp,ar,ac,xr,xc,p);
    }
}

void trsm(char SIDE, char UPLO, char TRANSA, 
          float **mata, float **matb, float **matx,
          int ar, int ac, int xr, int xc, int n ){

    float **temp = new_matrix(n, n);
    init(temp,n);
    float **tran_a = new_matrix(n, n);
    init(tran_a,n);

    

    int error = 0; 

    if(SIDE == 'L' || SIDE == 'l'){
        if(UPLO == 'L' ||UPLO == 'l'){
            if(TRANSA == 'N' || TRANSA == 'n'){
                woco_trsm(mata, matb, matx, temp, ar, ac, xr, xc, n);
            }
            else if(TRANSA == 'T' || TRANSA == 't' || TRANSA == 'C' || TRANSA == 'c'){
                transport(mata, tran_a, n, n);
                woco_trsm_u(tran_a, matb, matx, temp, ar, ac, xr, xc, n);
            }
            else error = 1;
        }
        else if(UPLO == 'U' ||UPLO == 'u'){
            if(TRANSA == 'N' || TRANSA == 'n'){
                woco_trsm_u(mata, matb, matx, temp, ar, ac, xr, xc, n);
            }
            else if(TRANSA == 'T' || TRANSA == 't' || TRANSA == 'C' || TRANSA == 'c'){
                transport(mata, tran_a, n, n);
                woco_trsm(tran_a, matb, matx, temp, ar, ac, xr, xc, n);
            }
            else error = 1;
        }
        else error = 2;
    }
    else if(SIDE == 'R' || SIDE == 'r'){
		
        float **tran_x = new_matrix(n, n);
        init(tran_x,n);
        float **tran_b = new_matrix(n, n);
        init(tran_b,n);

        if(UPLO == 'L' ||UPLO == 'l'){
            if(TRANSA == 'N' || TRANSA == 'n'){
                transport(mata, tran_a, n, n);
                transport(matb, tran_b, n ,n);
                woco_trsm_u(tran_a, tran_b, tran_x, temp, ar, ac, xr, xc, n);
                transport(tran_x, matx, n, n);
            }
            else if(TRANSA == 'T' || TRANSA == 't' || TRANSA == 'C' || TRANSA == 'c'){
                transport(matb, tran_b, n ,n);
                woco_trsm(mata, tran_b, tran_x, temp, ar, ac, xr, xc, n);
                transport(tran_x, matx, n, n);
            }
            else error = 1;
        }
        else if(UPLO == 'U' ||UPLO == 'u'){
            if(TRANSA == 'N' || TRANSA == 'n'){
                transport(mata, tran_a, n, n);
                transport(matb, tran_b, n ,n);
                woco_trsm(tran_a, tran_b, tran_x, temp, ar, ac, xr, xc, n);
                transport(tran_x, matx, n, n);
            }
            else if(TRANSA == 'T' || TRANSA == 't' || TRANSA == 'C' || TRANSA == 'c'){
                transport(matb, tran_b, n ,n);
                woco_trsm_u(mata, tran_b, tran_x, temp, ar, ac, xr, xc, n);
                transport(tran_x, matx, n, n);
            }
            else error = 1;
        }
        else error = 2;

        free(tran_b);
        free(tran_x);
    }
    else error = 3;

    free(temp);
    free(tran_a);

    if(error == 1) printf("argument TRANSA error!");
    if(error == 2) printf("argument UPLO error!");
    if(error == 3) printf("argument SIDE error!");

}

#endif //WONAL_TRSM_H


