//
// Created by whe on 2020/3/7
//

#ifndef WONAL_INV_H
#define WONAL_INV_H

#include<math.h>
#include "basewo.h"


void matrix_sign_inv(float** mata, int ar, int ac, int n){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            mata[ar + i][ac + j] = 0 - mata[ar + i][ac + j];
        }
    }
}

// A is upper triangular matrix
void inv(
        float** mata, float** inv_mata, float** temp,int ar, int ac, int m
        ){
    if(m == 1){
        if (0 == mata[ar][ac]) return ;//error, 不可求逆
        inv_mata[ar][ac] = 1 / mata[ar][ac];
    }else{
        int p = m / 2;
        inv(mata, inv_mata, temp, ar, ac, p);
        inv(mata, inv_mata, temp, ar + p, ac +p, p);

        vex_wo(mata, inv_mata, temp, ar, ac + p, ar + p, ac + p, p, p, p);
        matrix_sign_inv(temp, ar, ac + p, p);
        vex_wo(inv_mata, temp, inv_mata, ar, ac, ar, ac + p,  p, p, p);
    }
}



/*
// A is upper triangular matrix
void woco_inv(
    float** mata, float** inv_mata,float** temp, int ar, int ac, int m
    ){
    if(m * m < CACHE_SIZE){
        direct_inv(mata, inv_mata, temp, ar, ac, m);
    }else{
        int p = m / 2;
        woco_inv(mata, inv_mata, temp, ar, ac, p);
        woco_inv(mata, inv_mata, temp, ar + p, ac + p, p);
        
        vex_wo(mata, inv_mata, temp, ar, ac + p, ar + p, ac + p, p, p ,p);
        matrix_sign_inv(temp, ar, ac + p, p);
        

    }
}
*/
#endif //WONAL_INV_H