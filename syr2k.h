/*
 * @Author       : whe
 * @Date         : 2020-04-09 09:55:04
 * @LastEditTime : 2020-04-09 10:44:34
 * @Description  : 
 */
#ifndef WONAL_SYR2K_H
#define WONAL_SYR2K_H

#include "basewo.h"
#include "gemm.h"
/* 
DSYR2K  performs one of the symmetric rank 2k operations

    C := alpha*A*B**T + alpha*B*A**T + beta*C
    A is m * n, B is m * n, C is m * m
*/
// need optimize
void syr2k(float ** mata, float ** matb,float**matc,float alpha,float beta,int ar,int ac,int br,int bc, int m,int n){
    float ** tran_a = new_matrix(ac+n,ar+m);
    float ** tran_b = new_matrix(bc+n,br+m);
    transport(mata,tran_a,ar+m,ac+n);
    transport(matb,tran_b,br+m,bc+n);
    float ** temp1 = new_matrix(ar+m,br+m);
    float ** temp2 = new_matrix(br+m,ar+m);
    vex_wo(mata,tran_b,temp1,ar,ac,bc,br,m,n,m);
    vex_wo(matb,tran_a,temp2,br,bc,ac,ar,m,n,m);
    add(temp1,temp2,ar,br,m,m,alpha,alpha);
    add(matc,temp1,ar,br,m,m,1,beta);
    free(tran_a,ac+n);
    free(tran_b,bc+n);
    free(temp1,ar+m);
    free(temp2,br+m);
}

void test_syr2k(){
    printf("test for syr2k:\n");
    int m = 2048;
    int n = 2048;
    float alpha = 1;
    float beta = 1;
    float ** mata = new_matrix(m,n);
    float ** matb = new_matrix(m,n);
    float ** matc = new_matrix(m,m);
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0; j < n ; ++j){
            mata[i][j] = 1;
        }
    }
    for(int i = 0; i < m ; ++i){
        for (int j = 0; j < n; ++j){
            matb[i][j] = 2;
        }
    }
    init(matc,m,m);
    syr2k(mata,matb,matc,alpha,beta,0,0,0,0,m,n);
    display(matc,m,n);
    free(mata,m);
    free(matb,m);
    free(matc,m);
}


#endif