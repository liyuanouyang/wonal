/*
 * @Author       : whe
 * @Date         : 2020-04-18 09:32:38
 * @LastEditTime : 2020-04-18 10:08:16
 * @Description  : 
 */
#ifndef WONAL_POSV_H
#define WONAL_POSV_H

#include "Cholesky_fac.h"

/*
POSV computes the solution to a real system of linear equations 
    A * X = B,
where A is an M-by-M symmetric positive definite matrix and 
    X and B are M-by-N matrices.

The Cholesky decomposition is used to factor A as  
    A = L * L**T.
The factored form of A is then used to solve the system of equations 
    A * X = B.
*/

void posv(float ** mata, float ** matb, float ** matx, int ar, int ac, int br, int bc, int m, int n){
    float ** matl = new_matrix(m,m);
    float ** temp = new_matrix(m,n);
    float ** tran_l = new_matrix(m,m);
    float ** tempc = new_matrix(m,m);
    float ** templ = new_matrix(m,m);
    init(matl,m,m);
    init(temp,m,n);
    init(tempc,m,m);
    init(templ,m,m);
    init(tran_l,m,m);
    woco_chol(mata,matl,tran_l,tempc,temp,templ,ar,ac,m);
    //direct_Cho(mata,matl,ar,ac,m);
    transport(matl,tran_l,m,m);
    init(temp,m,n);
    trsm(matl,matb,temp,ar,ac,br,bc,m,n);
    trsm_u(tran_l,temp,matx,ar,ac,br,bc,m,n);
    free(matl,m);
    free(temp,m);
    free(tran_l,m);
    free(tempc,m);
    free(templ,m);
}

void test_posv(){
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
    posv(mata,matb,matx,0,0,0,0,m,n);
    display(matx,m,n);
}

#endif