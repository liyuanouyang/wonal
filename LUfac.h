//
// Created by oyliy on 2019/11/25.
//
//only use for square matrix
#ifndef WONAL_LUFAC_H
#define WONAL_LUFAC_H

#include "basewo.h"
#include "trsm.h"
#include <math.h>

void matrix_copy(float** mata, float** matb, int m,int n){
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            matb[i][j] = mata[i][j];
        }
    }
}

void _swap(float** p, float** q){
    float temp = **p;
    **p = **q;
    **q = temp;
}

void _sum_k(int a, int b, float c, float* s){
    *s = 0;
    for (int k = a; k < b; ++k) {
        *s+=c;
    }
}

//only for cube
void matrix_pivot(float** mata, float** matp,int ar,int ac, int n){
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matp[ar+i][ac+j] = (i==j);
        }
    }
    for (int i = 0; i < n; ++i) {
        int max_j = i;
        for (int j = i; j < n; ++j) {
            if (fabs(mata[ar+j][ac+i]) > fabs(mata[ar+max_j][ac+i]))    max_j = j;
        }
        if (max_j != i){
            for (int k = 0; k < n; ++k) {
//                _swap(matp[i][k], &matp[max_j][k]);
                int ari = ar+i,ack = ac+k,armaxj = ar+max_j;
                float temp = matp[ari][ack];
                matp[ari][ack] = matp[armaxj][ack];
                matp[armaxj][ack] = temp;
            }
        }
    }
}
//void matrix_pivot(float** mata, float** matp, int n){
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            matp[i][j] = (i==j);
//        }
//    }
//    for (int i = 0; i < n; ++i) {
//        int max_j = i;
//        for (int j = i; j < n; ++j) {
//            if (fabs(mata[j][i]) > fabs(mata[max_j][i]))    max_j = j;
//        }
//        if (max_j != i){
//            for (int k = 0; k < n; ++k) {
////                _swap(matp[i][k], &matp[max_j][k]);
//                float temp = matp[i][k];
//                matp[i][k] = matp[max_j][k];
//                matp[max_j][k] = matp[i][k];
//            }
//        }
//    }
//}

void lufac(float** mata, float** matl, float** matu,float** matp,int ar, int ac,int n,int inputm,int inputn){

    init(matl,ar,ac,n);init(matu,ar,ac,n);
    matrix_pivot(mata,matp,ar,ac,n);
    float ** matAprime = NULL;
    if (flag == 1){
        matAprime = new float*[inputm];
        for (int i = 0; i < inputm; ++i) {
            matAprime[i] = new float[inputn];
        }
        vex_wo(mata,matp,matAprime,ar,ac,ar,ac,n,n,n);
    }
    else if (flag == 2){matAprime=omp_vec_mat_mult_return(mata,matp,ar,ac,ar,ac,n,n,n,inputm,inputn);}

    for (int i = 0; i < n; ++i) {
        matl[ar+i][ac+i] = 1;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            float s = 0.0;
            if (j<=i){
                for (int k = 0; k < j; ++k){s += matl[ar+j][ac+k]*matu[ar+k][ac+i];}
                matu[ar+j][ac+i] = matAprime[ar+j][ac+i] - s;
            }
            if (j>i){
                for (int k = 0; k < i; ++k) {
                    s+=matl[ar+j][ac+k]*matu[ar+k][ac+i];
                }
                matl[ar+j][ac+i] = (matAprime[ar+j][ac+i] - s)/matu[ar+i][ac+i];
            }
        }
    }

}
//void lufac(float** mata, float** matl, float** matu,float** matp,int ar, int ac,int n){
//
//    init(matl,ar,ac,n);init(matu,ar,ac,n);
//    matrix_pivot(mata,matp,n);
//
//    float ** matAprime = vex_wo_return(mata,matp,0,0,0,0,n,n,n);
//
//    for (int i = 0; i < n; ++i) {
//        matl[i][i] = 1;
//    }
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            float s = 0.0;
//            if (j<=i){
//                for (int k = 0; k < j; ++k){s += matl[j][k]*matu[k][i]}
//                matu[j][i] = matAprime[j][i] - s;
//            }
//            if (j>i){
//                for (int k = 0; k < i; ++k) {
//                    s+=matl[j][k]*matu[k][i];
//                }
//                matl[j][i] = (matAprime[j][i] - s)/matu[i][i];
//            }
//        }
//    }
//
//}

void woco_LUfac(float** mata,float** matl,float** matu,float**matp,float** trsm_temp,int ar,int ac,int n,int inputm,int inputn){
    //p!=q cant be factorize only odd can be resicuve
    int p = n/2;
    int q = n-p;
    if (n*n<CACHE_SIZE|| p!=q){
        lufac(mata,matl,matu,matp,ar,ac,n,inputm,inputn);
    }
    woco_LUfac(mata,matl,matu,matp,trsm_temp,ar,ac,p,inputm,inputn);
    woco_trsm(matl,mata,matu,trsm_temp,ar,ac,ar,ac+p,q,q);
    woco_trsm(matu,mata,matl,trsm_temp,ar,ac,ar+p,ac,q,q);
    if (flag ==1){
        float **temp = new float*[inputm];
        for (int i = 0; i < inputm; ++i) {
            temp[i] = new float[inputn];
        }
        vex_wo(matl,matu,temp,ar+p,ac,ar,ac+p,p,p,p);
        matrix_minus(mata,temp,ar,ac,0,0,p);
    } else if(flag == 2){
        matrix_minus(mata,omp_vec_mat_mult_return(matl,matu,ar+p,ac,ar,ac+p,p,p,p,inputm,inputn),ar,ac,0,0,p);
    }
    woco_LUfac(mata,matl,matu,matp,trsm_temp,ar+p,ac+p,q,inputm,inputn);
}

#endif //WONAL_LUFAC_H
