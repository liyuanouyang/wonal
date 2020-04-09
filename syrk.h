/*
 * @Author       : whe
 * @Date         : 2020-04-09 09:54:05
 * @LastEditTime : 2020-04-09 10:46:09
 * @Description  : 
 */
#ifndef WONAL_SYRK_H
#define WONAL_SYRK_H

#include "basewo.h"
#include "gemm.h"

/*
 SYRK  performs one of the symmetric rank k operations

C := alpha*A*A**T + beta*C

A is m * n, C is m * m
*/    
// need optimize
void syrk(float ** mata, float **matc, float alpha,float beta,int ar, int ac, int br, int bc, int m,int n){
    float ** tran_a = new_matrix(ac+n,ar+m);
    transport(mata,tran_a,ar+m,ac+n);
    gemm(mata,tran_a,matc,alpha,beta,ar,ac,br,bc,m,n,m);
    free(tran_a,ac+n);
}

void test_syrk(){
    printf("test for syrk:\n");
    int m = 2048;
    int n = 2048;
    float alpha = 1;
    float beta = 1;
    float ** mata = new_matrix(m,n);
    float ** matc = new_matrix(m,m);
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0; j < m ; ++j){
            mata[i][j] = 1;
        }
    }
    init(matc,m,n);
    syrk(mata,matc,alpha,beta,0,0,0,0,m,n);
    display(matc,m,n);
    free(mata,m);
    free(matc,m);
}

#endif