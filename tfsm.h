/*
 * @Author       : whe
 * @Date         : 2020-04-03 22:26:06
 * @LastEditTime : 2020-04-06 19:05:46
 * @Description  : 
 */

#ifndef WONAL_TFSM_H
#define WONAL_TFSM_H

#include "trsm.h"
#include "gemm.h"

/*
 Level 3 BLAS like routine for A in RFP Format.
 DTFSM  solves the matrix equation A *X = B  
 A is in Rectangular Full Packed (RFP) Format.
*/
//just for m is odd
void tfsm(float ** mata, float ** matb,float** matx,int ar, int ac, int br, int bc, int m, int n){
    int m1,m2;
    
    m2 = m / 2;
    m1 = m -m2;
    
    if(m1 == 1){
        trsm(mata,matb,matx,ar,ac,br,bc,m1,n);
    }
    else{
        trsm(mata,matb,matx,ar,ac,br,bc,m1,n);
        gemm(mata,matx,matb,-1,1,ar+m1,ac,br,bc,m2,m1,n);
        float** tran = new_matrix(m2,m2);
        transport(mata,tran,ar,ac+1,m2,m2);
        trsm(tran,matb,matx,0,0,br+m1,bc,m2,n);
        free(tran,m2);
    }        
}

void test_tfsm(){
    int m = 2049;
    int k = m / 2 + 1;   
    int n = 2048;
    float ** mata = new_matrix(m,k);
    float ** matb = new_matrix(m,n);
    float ** matx = new_matrix(m,n);
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0; j < k ; ++j){
            mata[i][j] = 1;
        }
    }
    for(int i = 0; i < m ; ++i){
        for (int j = 0; j < n; ++j){
            matb[i][j] = (i + 1)*(j + 1);
        }
    }
    init(matx,m,n);
    //trsm(mata,matb,matx,0,0,0,0,k,n);
    tfsm(mata,matb,matx,0,0,0,0,m,n);
    display(matx,m,n);
    free(mata,m);
    free(matb,m);
    free(matx,m);
}

#endif
