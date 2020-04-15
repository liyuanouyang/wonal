/*
 * @Author       : whe
 * @Date         : 2020-04-14 09:30:04
 * @LastEditTime : 2020-04-14 10:38:58
 * @Description  : POTRS
 */
#ifndef WONAL_POTRS_H
#define WONAL_POTRS_H

#include "basewo.h"
#include "trsm.h"

/*
DPOTRS solves a system of linear equations A*X = B with a symmetric
 positive definite matrix A using the Cholesky factorization
 A = U**T*U or A = L*L**T computed by POTRF.
 A is m * m, b is m * n.
 */
void potrs(float ** mata,float ** matb, float ** matx, int ar, int ac, int br ,int bc, int m, int n){
    float ** temp = new_matrix(m,n);
    init(temp,m,n);
    float ** tran_a = new_matrix(m,m);
    transport(mata,tran_a,m,m);
    trsm(mata,matb,temp,ar,ac,br,bc,m,n);
    trsm_u(tran_a,temp,matx,ar,ac,br,bc,m,n);
}

void test_potrs(){
    int m = 2048;
    int n = 2048;
    float ** mata = new_matrix(m,m);
    float ** matb = new_matrix(m,n);
    float ** matx = new_matrix(m,n);
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0; j < m ; ++j){
            if(i >= j) mata[i][j] = 1;
            else mata[i][j] = 0;
        }
    }
    for(int i = 0; i < m ; ++i){
        for (int j = 0; j < n; ++j){
            matb[i][j] = (i + 1)*(j + 1);
        }
    }
    init(matx,m,n);
    potrs(mata,matb,matx,0,0,0,0,m,n);
    display(matx,m,n);
    free(mata,m);
    free(matb,m);
    free(matx,m);
}

#endif