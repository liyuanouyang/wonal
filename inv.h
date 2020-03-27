//
// Created by whe on 2020/3/7
//

#ifndef WONAL_INV_H
#define WONAL_INV_H

#include<math.h>
#include "basewo.h"
#include "LUfac.h"

void matrix_sign_inv(float** mata, int ar, int ac, int n){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            mata[ar + i][ac + j] = 0 - mata[ar + i][ac + j];
        }
    }
}

// A is upper triangular matrix
void inv_u(float** mata, float** inv_mata,float** temp, int ar, int ac, int m){
    if(m == 1){
        inv_mata[ar][ac] = 1 / mata[ar][ac];
    }else{
        int p = m / 2;
        inv_u(mata, inv_mata,temp, ar, ac, p);
        inv_u(mata, inv_mata,temp, ar+p, ac +p, p);
        vex_wo(inv_mata, mata, temp, ar, ac, ar, ac + p, p, p ,p);
        matrix_sign_inv(temp, ar, ac + p, p);
        vex_wo(temp, inv_mata, inv_mata, ar, ac + p, ar + p, ac + p, p, p, p);
    }
}
// A is lower triangular matrix
void inv_l(float** mata, float** inv_mata,float** temp, int ar, int ac, int m){ 
    if(m == 1){
        inv_mata[ar][ac] = 1 / mata[ar][ac];
    }else{
        int p = m / 2;
        inv_l(mata, inv_mata, temp, ar, ac, p);
        inv_l(mata, inv_mata, temp, ar + p, ac + p, p);
        vex_wo(inv_mata, mata, temp, ar + p, ac + p, ar + p, ac, p, p, p);
        matrix_sign_inv(temp, ar + p, ac, p);
        vex_wo(temp, inv_mata, inv_mata, ar + p, ac, ar, ac, p, p, p);
    }
}
// A is m*m matrix, need test
void inv(
        float** mata, float** inv_mata, float** temp,int ar, int ac, int m
        ){
    float ** matl = NULL;
    float ** matu = NULL;
    float ** matp = NULL;
    float ** inv_matl = NULL;
    float ** inv_matu = NULL;
    square_matrix(&matl, m);
    square_matrix(&matu, m);
    square_matrix(&matp, m);
    square_matrix(&inv_matl,m);
    square_matrix(&inv_matu,m);
    init(matl,m);
    init(matu,m);
    init(matp,m);
    init(inv_matl,m);
    init(inv_matu,m);
    lufac(mata, matl, matu, matp, ar, ac, m, m ,m);
    inv_l(matl, inv_matl, temp, ar, ac, m);
    inv_u(matu, inv_matu, temp, ar, ac, m);
    vex_wo(inv_matl, inv_matu, inv_mata, ar, ac, ar, ac, m, m, m);
}

#endif //WONAL_INV_H
