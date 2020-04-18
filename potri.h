/*
 * @Author       : whe
 * @Date         : 2020-04-14 10:22:11
 * @LastEditTime : 2020-04-18 12:09:11
 * @Description  : 
 */
#ifndef WONAL_POTRI_H
#define WONAL_POTRI_H

#include "basewo.h"
#include "inv.h"
#include "Cholesky_fac.h"
/*
POTRI computes the inverse of a real symmetric positive definite
 matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
 */
// just for A = L*L**T
void potri(float ** mata, float ** matx, int ar, int ac, int n){
    float ** matl = new_matrix(n,n);
    float ** temp = new_matrix(n,n);
    float ** tran_l = new_matrix(n,n);
    float ** tempc = new_matrix(n,n);
    float ** templ = new_matrix(n,n);
    float ** invl = new_matrix(n,n);
    float ** invu = new_matrix(n,n);
    init(matl,n,n);
    init(temp,n,n);
    init(tempc,n,n);
    init(templ,n,n);
    init(tran_l,n,n);
    init(invl,n,n);
    init(invu,n,n);
    woco_chol(mata,matl,tran_l,tempc,temp,templ,ar,ac,n);
    transport(matl,tran_l,n,n);
    woco_invl(matl,invl,temp,ar,ac,n);
    init(temp,n,n);
    woco_invu(tran_l,invu,temp,ar,ac,n);
    vex_wo(invu,invl,matx,ar,ac,ar,ac,n,n,n);
    free(invl,n);
    free(invu,n);
    free(tran_l,n);
    free(temp,n);
    free(matl,n);
    free(tempc,n);
    free(templ,n);
}

void test_potri(){
    int m = 2048;
    float ** mata = new_matrix(m,m);
    float ** matx = new_matrix(m,m);
    init(matx,m,m);
    for(int i = 0; i < m; i++){
        mata[i][i] = i + 1;
        for(int j = 0; j < m; j++){
            if(i > j) mata[i][j] = mata[j][j];
            if(i < j) mata[i][j] = mata[i][i];
        }
    }
    potri(mata,matx,0,0,m);
    display(matx,m,m);
    free(mata,m);
    free(matx,m);
}

#endif