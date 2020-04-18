/*
 * @Author       : whe
 * @Date         : 2020-04-14 10:30:52
 * @LastEditTime : 2020-04-18 11:39:12
 * @Description  :
 */
#ifndef WONAL_GESV_H
#define WONAL_GESV_H

#include "basewo.h"
#include "LUfac.h"
#include "trsm.h"

/*
   GESV computes the solution to a real system of linear equations
   A * X = B,
   where A is an N-by-N matrix and X and B are N-by-NRHS matrices.

   The LU decomposition with partial pivoting and row interchanges is
   used to factor A as
   A = P * L * U,
   where P is a permutation matrix, L is unit lower triangular, and U is
   upper triangular.  The factored form of A is then used to solve the
   system of equations A * X = B
   */
// A is m * m, B is m * n.
void gesv(float ** mata, float ** matb, float ** matx, int ar, int ac, int br, int bc, int m, int n){
    float ** matp = new_matrix(m,m);
    float ** matl = new_matrix(m,m);
    float ** matu = new_matrix(m,m);
    float ** temp = new_matrix(m,n);
    init(matp,m);
    init(matl,m);
    init(matu,m);
    init(temp,m,n);
    lufac(mata,matl,matu,matp,ar,ac,m,m,m);
    trsm(matl,matb,matx,ar,ac,br,bc,m,n);
    trsm_u(matu,matx,temp,ar,ac,br,bc,m,n);
    init(matx,m,n);
    vex_wo(matp,temp,matx,ar,ac,br,bc,m,m,n);
    free(matp,m);
    free(matl,m);
    free(matu,m);
    free(temp,m);
}

void test_gesv(){
    int m = 2048;
    int n = 2048;
    float ** mata = new_matrix(m,m);
    float ** matb = new_matrix(m,n);
    float ** matx = new_matrix(m,n);
    init(matb,m,n);
    init(matx,m,n);
    for(int i = 0; i < m; i++){
        mata[i][i] = i + 1;
        for(int j = 0; j < m; j++){
            if(i > j) mata[i][j] = mata[j][j];
            if(i < j) mata[i][j] = mata[i][i];
        }
    }
    for(int i = 0; i < m ; ++i){
        for (int j = 0; j < n; ++j){
            matb[i][j] = (i + 1)*(j + 1);
        }
    }
    gesv(mata,matb,matx,0,0,0,0,m,n);
    display(matx,m,n);
}

#endif