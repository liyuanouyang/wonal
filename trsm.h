//
// Created by oyliy on 2019/11/25.
//

#ifndef WONAL_TRSM_H
#define WONAL_TRSM_H


#include "basewo.h"

// matrix minus method
void matrix_minus(float **mata,float **matb,int ar,int ac,int br, int bc,int m,int n){
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            mata[ar+i][ac + j] -= matb[br+i][bc+j];
        }
    }
}


void direct_trsm(float **mata, float **matb, float **matx,
        int ar, int ac, int xr, int xc,
        int m, int n){
    //A is lower triangular matrix
    //A is m*m , B is m*n X is m*n ,wo find x
#pragma omp parallel for num_threads(core_num) default(none) shared(mata,matb,matx,m,n,ar,ac,xr,xc)
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            int temp = 0;
            for (int k = 0; k < i; ++k) {
                temp += mata[ar+i][ac+k]*matx[xr+k][xc+j];
            }
            matx[xr+i][xc+j] = (matb[ar+i][xc+j] - temp)/mata[ar+i][ac+i];
            temp = 0;
        }
    }
}


//need test
//A is m*m , B is m*n X is m*n ,wo find x
void woco_trsm(float **mata, float **matb, float **matx, float** temp,
        int ar, int ac,int xr, int xc, int m, int n){
    if (3*m*n < CACHE_SIZE){
        direct_trsm(mata,matb,matx,ar,ac,xr,xc,m,n);
    } else{
        //partitiion to 4
        int r = m / 2;
        int c = n / 2;
        woco_trsm(mata,matb,matx,temp,ar,ac,xr,xc,r,c);
        woco_trsm(mata,matb,matx,temp,ar,ac,xr,xc + c,r,c);
        vex_wo(mata,matx,temp,ar+r,ac,xr,xc,r,r,c);
        //        if (flag == 1){vex_wo(mata,matx,temp,ar+r,ac,xr,xc,r,r,c);}
        //        else if(flag == 2){omp_vec_mat_mult(mata,matx,temp,ar+r,ac,xr,xc,r,r,c);}
        //        else if(flag == 3){matrix_multi(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
        //       matrix_minus(matb,temp,ar+q,xc,ar+q,xc,p);
        matrix_minus(matb,temp,xr+r,xc,ar+r,xc,r,c);
        woco_trsm(mata,matb,matx,temp,ar+r,ac+r,xr+r,xc,r,c);
        vex_wo(mata,matx,temp,ar+r,ac,xr,xc+c,r,r,c);
        //        if (flag ==1){vex_wo(mata,matx,temp,ar+r,ac,xr,xc+c,r,r,c);}
        //        else if (flag == 2){omp_vec_mat_mult(mata,matx,temp,ar+r,ac,xr,xc+c,r,r,c);}
        //        else if(flag == 3){matrix_multi(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
        matrix_minus(matb,temp,xr+r,xc+c,ar+r,xc+c,r,c);
        woco_trsm(mata,matb,matx,temp,ar+r,ac+r,xr+r,xc+c,r,c);
    }
}


void direct_trsm_u(float **mata, float **matb, float **matx,
        int ar, int ac, int xr, int xc,
        int m, int n){
    //A is upper triangular matrix
    //A is m*m , B is m*n X is m*n ,wo find x
    //square n
#pragma omp parallel for num_threads(core_num) default(none) shared(mata,matb,matx,m,n,ar,ac,xr,xc)
    for (int j = 0; j < n; ++j) {
        for (int i = m - 1; i >= 0; --i) {
            int temp = 0;
            for (int k = i; k < m ; ++k) {
                temp += mata[ar+i][ac+k]*matx[xr+k][xc+j];
            }
            matx[xr+i][xc+j] = (matb[ar+i][xc+j] - temp)/mata[ar+i][ac+i];
            temp = 0;
        }
    }
}

//need test
void woco_trsm_u(float **mata, float **matb, float **matx, float** temp,
        int ar, int ac,int xr, int xc, int m, int n){
    if (3*m*n < CACHE_SIZE){
        direct_trsm_u(mata, matb, matx, ar, ac, xr, xc, m, n);
    } else{
        //partitiion to 4
        int r = m/2;
        int c = n/2;
        woco_trsm_u(mata, matb, matx, temp, ar + r, ac + r, xr + r, xc + c , r , c);
        woco_trsm_u(mata, matb, matx, temp, ar + r, ac + r, xr + r, xc, r, c);
        vex_wo(mata,matx,temp,ar,ac+r,xr+r,xc+c,r,r,c);
        //        if (flag == 1){vex_wo(mata,matx,temp,ar,ac+r,xr+r,xc+c,r,r,c);}
        //        else if(flag == 2){omp_vec_mat_mult(mata,matx,temp,ar,ac+r,xr+r,xc+c,r,r,c);}
        //        else if(flag == 3){matrix_multi(mata,matx,temp,ar,ac+r,xr+r,xc+c,r,r,c);}
        matrix_minus(matb,temp,xr,xc+c,ar,xc+c,r,c);
        woco_trsm_u(mata,matb,matx,temp,ar,ac,xr,xc+c,r,c);
        vex_wo(mata,matx,temp,ar,ac+r,xr+r,xc,r,r,c);
        //        if (flag ==1){vex_wo(mata,matx,temp,ar,ac+r,xr+r,xc,r,r,c);}
        //        else if (flag == 2){omp_vec_mat_mult(mata,matx,temp,ar,ac+r,xr+r,xc,r,r,c);}
        //        else if(flag == 3){matrix_multi(mata,matx,temp,ar,ac+r,xr+r,xc,r,r,c);}
        matrix_minus(matb,matx,xr,xc,ar,xc,r,c);
        woco_trsm_u(mata,matb,matx,temp,ar,ac,xr,xc,r,c);
    }
}

/*
DTRSM  solves one of the matrix equations

    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,

 where alpha is a scalar, X and B are m by n matrices, A is a unit, or
 non-unit,  upper or lower triangular matrix  and  op( A )  is one  of

    op( A ) = A   or   op( A ) = A**T.

 The matrix X is overwritten on B.
*/
void trsm(char SIDE, char UPLO, char TRANSA,
        float **mata, float **matb, float alpha,
        int ar, int ac, int xr, int xc, int m, int n ){
	float **matx = new_matrix(m, n);
    init(matx,m,n);
    int error = 0;
	mul(matb,m,n,alpha);
    if(SIDE == 'L' || SIDE == 'l'){

        float **temp = new_matrix(m, n);
        init(temp,m,n);
        float **tran_a = new_matrix(m, m);
        init(tran_a,m,m);

        if(UPLO == 'L' ||UPLO == 'l'){
            if(TRANSA == 'N' || TRANSA == 'n'){
                woco_trsm(mata, matb, matx, temp, ar, ac, xr, xc, m, n);
            }
            else if(TRANSA == 'T' || TRANSA == 't' || TRANSA == 'C' || TRANSA == 'c'){
                transport(mata, tran_a, m, m);
                woco_trsm_u(tran_a, matb, matx, temp, ar, ac, xr, xc, m, n);
            }
            else error = 1;
        }
        else if(UPLO == 'U' ||UPLO == 'u'){
            if(TRANSA == 'N' || TRANSA == 'n'){
                woco_trsm_u(mata, matb, matx, temp, ar, ac, xr, xc, m, n);
            }
            else if(TRANSA == 'T' || TRANSA == 't' || TRANSA == 'C' || TRANSA == 'c'){
                transport(mata, tran_a, m, m);
                woco_trsm(tran_a, matb, matx, temp, ar, ac, xr, xc, m, n);
            }
            else error = 1;
        }
        else error = 2;

        free(tran_a,m);
        free(temp,m);
    }
    else if(SIDE == 'R' || SIDE == 'r'){

        float **temp = new_matrix(n, m);
        init(temp,n,m);
        float **tran_a = new_matrix(n, n);
        init(tran_a,n,n);
        float **tran_x = new_matrix(n, m);
        init(tran_x,n,m);
        float **tran_b = new_matrix(n, m);
        init(tran_b,n,m);

        if(UPLO == 'L' ||UPLO == 'l'){
            if(TRANSA == 'N' || TRANSA == 'n'){
                transport(mata, tran_a, n, n);
                transport(matb, tran_b, m ,n);
                woco_trsm_u(tran_a, tran_b, tran_x, temp, ar, ac, xr, xc, n, m);
                transport(tran_x, matx, n, m);
            }
            else if(TRANSA == 'T' || TRANSA == 't' || TRANSA == 'C' || TRANSA == 'c'){
                transport(matb, tran_b, m ,n);
                woco_trsm(mata, tran_b, tran_x, temp, ar, ac, xr, xc, n, m);
                transport(tran_x, matx, n, m);
            }
            else error = 1;
        }
        else if(UPLO == 'U' ||UPLO == 'u'){
            if(TRANSA == 'N' || TRANSA == 'n'){
                transport(mata, tran_a, n, n);
                transport(matb, tran_b, m ,n);
                woco_trsm(tran_a, tran_b, tran_x, temp, ar, ac, xr, xc, n, m);
                transport(tran_x, matx, n, m);
            }
            else if(TRANSA == 'T' || TRANSA == 't' || TRANSA == 'C' || TRANSA == 'c'){
                transport(matb, tran_b, m ,n);
                woco_trsm_u(mata, tran_b, tran_x, temp, ar, ac, xr, xc, n, m);
                transport(tran_x, matx, n, m);
            }
            else error = 1;
        }
        else error = 2;

        free(tran_b,n);
        free(tran_x,n);
        free(tran_a,n);
        free(temp,n);
    }
    else error = 3;

    if(error == 1) printf("argument TRANSA error!");
    if(error == 2) printf("argument UPLO error!");
    if(error == 3) printf("argument SIDE error!");
	copy(matb,matx,m,n);
	free(matx,m);

}


#endif //WONAL_TRSM_H


