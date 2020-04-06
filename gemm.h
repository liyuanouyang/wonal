/*
 * @Author       : whe
 * @Date         : 2020-04-04 18:10:29
 * @LastEditTime : 2020-04-06 19:05:58
 * @Description  : 
 */


#ifndef WONAL_GEMM_H
#define WONAL_GEMM_H

#include "basewo.h"
#include <iostream>
using namespace std;
/*
   GEMM  performs one of the matrix-matrix operations

C := alpha* A * B  + beta*C,


*/
void gemm(float **mata, float **matb, float **matc,float alpha,float beta,
        int ar, int ac, int br, int bc, int m, int k, int n){

    float ** temp = new_matrix(ar + m, bc + n);
    init(temp,ar + m, bc + n);

    // A is m * k, B is k * n
    vex_wo(mata,matb,temp,ar,ac,br,bc,m,k,n);
    add(matc,temp,ar,bc,m,n,beta,alpha);
    free(temp,m);
}

void test_gemm(){
    printf("test for gemm:\n");
    int m = 2048;
    int n = 2048;
    int k = 2048;
    float alpha = 1;
    float beta = 1;
    float ** mata = new_matrix(m,k);
    float ** matb = new_matrix(k,n);
    float ** matc = new_matrix(m,n);
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0; j < m ; ++j){
            mata[i][j] = 1;
        }
    }
    for(int i = 0; i < m ; ++i){
        for (int j = 0; j < n; ++j){
            matb[i][j] = 2;
        }
    }
    init(matc,m,n);
    gemm(mata,matb,matc,alpha,beta,0,0,0,0,m,k,n);
    display(matc,m,n);
    free(mata,m);
    free(matb,k);
    free(matc,m);
}

#endif
