//
// Created by oyliy on 2019/11/25.
//

#ifndef WONAL_TRSM_H
#define WONAL_TRSM_H


#include "basewo.h"

//square matrix minus method
void matrix_minus(float **mata,float **matb,int ar,int ac,int br, int bc,int n){
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mata[ar+i][ac + j] -= matb[br+i][bc+j];
        }
    }
}


void direct_trsm(float **mata, float **matb, float **matx,
                 int ar, int ac, int xr, int xc,
                 int n){
    //A is m*k , B is m*n X is k*n ,wo find x
    //square n
#pragma omp parallel for num_threads(core_num) default(none) shared(mata,matb,matx,n,ar,ac,xr,xc)
    for (int i = 0; i < n; ++i) {
        int temp = 0;
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < i; ++k) {
                temp = mata[ar+i][ac+k]*matx[xr+k][xc+j];
            }
            matx[xr+i][xc+j] = (matb[ar+i][xc+j] - temp)/mata[ar+i][ac+i];
            temp = 0;
        }
    }
}

//AX = B
void woco_trsm(float **mata, float **matb, float **matx, float** temp,
               int ar, int ac,int xr, int xc, int n){
    if (3*n*n < CACHE_SIZE){
        direct_trsm(mata,matb,matx,ar,ac,xr,xc,n);
    } else{
        //partitiion to 4
        int p = n/2;
        woco_trsm(mata,matb,matx,temp,ar,ac,xr,xc,p);
        woco_trsm(mata,matb,matx,temp,ar,ac,xr,xc+p,p);
        if (flag == 1){vex_wo(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
        else if(flag == 2){omp_vec_mat_mult(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
//        else if(flag == 3){matrix_multi(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
        matrix_minus(matb,temp,ar+p,xc,ar+p,xc,p);
        woco_trsm(mata,matb,matx,temp,ar+p,ac+p,xr+p,xc,p);
        if (flag ==1){vex_wo(mata,matx,temp,ar+p,ac,xr,xc+p,p,p,p);}
        else if (flag == 2){omp_vec_mat_mult(mata,matx,temp,ar+p,ac,xr,xc+p,p,p,p);}
//        else if(flag == 3){matrix_multi(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
        matrix_minus(matb,matx,ar+p,xc+p,ar+p,xc+p,p);
        woco_trsm(mata,matb,matx,temp,ar+p,ac+p,xr+p,xc+p,p);
    }
}




#endif //WONAL_TRSM_H


