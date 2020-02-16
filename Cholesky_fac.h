//
// Created by oyliy on 2019/11/25.
//

#ifndef WONAL_CHOLESKY_FAC_H
#define WONAL_CHOLESKY_FAC_H

#include <math.h>
#include "trsm.h"
#include "basewo.h"
//all ar,ac was be right,first,second position
void direct_Cho(
        float** mata, float** matl ,int ar , int ac, int m
        ){
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < i; ++j) {
            float plusres = 0.0;
            int ari = ar+i;
            int arj = ar+j;
            int acj = ac+j;
            for (int k = 0; k < j; ++k) {
                plusres += matl[ari][ac+k]*matl[arj][ac+k];
            }
            matl[ari][acj] = (mata[ari][acj] - plusres)/matl[arj][acj];
        }
        float multres = 0.0;
        for (int k = 0; k < i; ++k) {
            multres += matl[ar+i][ac+k]*matl[ar+i][ac+k];
        }
        matl[ar+i][ac+i] = mata[ar+i][ac+i] - multres;
    }
}

//void direct_Cho(
//        float** mata, float** matl ,int ar , int ac, int m
//){
//    for (int i = 0; i < m; ++i) {
//        for (int j = 0; j < i; ++j) {
//            float plusres = 0.0;
//            int ari = ar+i;
//            for (int k = 0; k < j; ++k) {
//                plusres += matl[i][k]*matl[j][k];
//            }
//            matl[i][j] = (mata[i][j] - plusres)/matl[j][j];
//        }
//        float multres = 0.0;
//        for (int k = 0; k < i; ++k) {
//            multres += matl[i][k]*matl[i][k];
//        }
//        matl[i][i] = mata[i][i] - multres;
//    }
//}

void woco_chol(float** mata, float** matl, float** matlt, float** tempc, float** temp, float** templ, int ar, int ac, int n){
    if (n*n < CACHE_SIZE){
        direct_Cho(mata,matl,ar,ac,n);
    } else{
        int m = n/2;
        int p = n-m;
        woco_chol(mata,matl,matlt,tempc,temp,templ,ar,ac,m);
        init(temp,ar,ac,m);
        woco_trsm(matl,mata,matlt,temp,ar,ac,ar,ac+m,p);
        transport(matl,templ,p,m);
//        float** tempc;
//        square_matrix(&tempc,2*p);
        init(tempc,ar,ac,n);
        if (flag == 1){vex_wo(matl,templ,tempc,ar+m,ac,ar,ac+m,p,m,p);}
        else if (flag == 2){omp_vec_mat_mult(matl,templ,tempc,ar+m,ac,ar,ac+m,p,m,p);}
        matrix_minus(mata,tempc,ar+m,ac+m,0,0,p);
        woco_chol(mata,matl,matlt,tempc,temp,templ,ar+m,ac+m,p);
    }
}



#endif //WONAL_CHOLESKY_FAC_H
