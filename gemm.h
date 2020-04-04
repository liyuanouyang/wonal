#ifndef WONAL_GEMM_H
#define WONAL_GEMM_H

#include "basewo.h"

/*
   GEMM  performs one of the matrix-matrix operations

C := alpha*op( A )*op( B ) + beta*C,

where  op( X ) is one of

op( X ) = X   or   op( X ) = X**T,
*/

void gemm(char TRANSA, char TRANSB,
        float **mata, float **matb, float **matc,float alpha, float beta,
        int ar, int ac, int br, int bc, int m, int k, int n){

    float ** temp = new_matrix(m,n);
    copy(temp,matc,m,n);

    // A is m * k, B is k * n
    if((TRANSA == 'N' || TRANSA == 'n') && (TRANSB == 'N' || TRANSB == 'n')){
        vex_wo(mata,matb,matc,ar,ac,br,bc,m,k,n);
    }
    //A is k * m, B is k * n
    else if((TRANSA == 'T' || TRANSA == 't' || TRANSA == 'C' || TRANSA == 'c') && (TRANSB == 'N' || TRANSB == 'n')){

        float ** tran_a = new_matrix(m,k);
        init(tran_a,m,k);
        transport(mata,tran_a,k,m);
        vex_wo(tran_a,matb,matc,ar,ac,br,bc,m,k,n);
        free(tran_a, m);
    }
    // A is m * k, B is n * k
    else if((TRANSB == 'T' || TRANSB == 't' || TRANSB == 'C' || TRANSB == 'c') && (TRANSA == 'N' || TRANSA == 'n')){

        float ** tran_b = new_matrix(k,n);
        init(tran_b,k,n);
        transport(matb,tran_b,n,k);
        vex_wo(mata,tran_b,matc,ar,ac,br,bc,m,n,k);
        free(tran_b,k);

    }
    //A is k * m, B is n * k
    else if((TRANSA == 'T' || TRANSA == 't' || TRANSA == 'C' || TRANSA == 'c') && (TRANSB == 'T' || TRANSB == 't' || TRANSB == 'C' || TRANSB == 'c')){
        float ** tran_a = new_matrix(m,k);
        float ** tran_b = new_matrix(k,n);
        init(tran_a,m,k);
        init(tran_b,k,n);
        transport(mata,tran_a,k,m);
        transport(matb,tran_b,n,k);
        vex_wo(tran_a,tran_b,matc,ar,ac,br,bc,m,k,n);
        free(tran_a,m);
        free(tran_b,n);
    }
    else{
        add(matc,temp,m,n,alpha,beta);
        free(temp,m);
        printf("argument error");
    }

    add(matc,temp,m,n,alpha,beta);
    free(temp,m);
}

#endif
