/*
 * @Author       : whe
 * @Date         : 2020-03-07 16:24:20
 * @LastEditTime : 2020-04-07 10:41:00
 * @Description  :
 */

#ifndef WONAL_INV_H
#define WONAL_INV_H

#include<math.h>
#include "basewo.h"
#include "LUfac.h"


// A is lower triangular matrix
void invl(float** mata,float** inv_mata,int ar,int ac,int n){
    for ( int i = 0; i < n; i++ )
    {
        if ( mata[ar + i][ac + i] == 0.0 )
        {
            cout << "*** Singular matrix ***\n";
            return ;
        }
        inv_mata[ar+i][ac+i] = 1.0 / mata[ar+i][ac+i];
        for ( int j = 0; j < i; j++ )
        {
            for ( int k = j; k < i; k++ ) inv_mata[ar+i][ac+j] += mata[ar+i][ac+k] * inv_mata[ac+k][ac+j];
            inv_mata[ar+i][ac+j] = -inv_mata[ar+i][ac+j] / mata[ar+i][ac+i];
        }
    }
}

void matrix_sign_inv(float** mata, int ar, int ac, int n){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            mata[ar + i][ac + j] = 0 - mata[ar + i][ac + j];
        }
    }
}

void woco_invl(float** mata, float** inv_mata,float** temp, int ar, int ac, int m){
    if(3 *m *m < CACHE_SIZE || m % 2 == 1){
        invl(mata,inv_mata,ar,ac,m);
    }else{
        int p = m / 2;
        woco_invl(mata, inv_mata, temp, ar, ac, p);
        woco_invl(mata, inv_mata, temp, ar + p, ac + p, p);
        init(temp,ar+m,ac+m);
        vex_wo(inv_mata, mata, temp, ar + p, ac + p, ar + p, ac, p, p, p);
        matrix_sign_inv(temp, ar + p, ac, p);
        vex_wo(temp, inv_mata, inv_mata, ar + p, ac, ar, ac, p, p, p);
    }
}

void test_invl(){
    int m = 2048;
    float ** mata = new_matrix(m,m);
    float ** matb = new_matrix(m,m);
    float ** temp = new_matrix(m,m);
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0; j < m ; ++j){
            if(i >= j) mata[i][j] = 1;
            else mata[i][j] = 0;
        }
    }
    init(matb,m,m);
    init(temp,m,m);
    woco_invl(mata,matb,temp,0,0,m);
    display(matb,m,m);
}

// A is upper triangular matrix

void invu(float** mata,float** inv_mata,int ar,int ac,int n){
    for ( int i = n - 1; i >= 0; i--)
    {
        if ( mata[ar + i][ac + i] == 0.0 )
        {
            cout << "*** Singular matrix ***\n";
            return ;
        }
        inv_mata[ar+i][ac+i] = 1.0 / mata[i][i];
        for ( int j = n - 1; j > i; j-- )
        {
            for ( int k = j; k > i; k-- ) inv_mata[ar+i][ac+j] += mata[ar+i][ac+k] * inv_mata[ar+k][ac+j];
            inv_mata[ar+i][ac+j] = -inv_mata[ar+i][ac+j] / mata[ar+i][ac+i];
        }
    }
}

void woco_invu(float** mata, float** inv_mata,float** temp, int ar, int ac, int m){
    if(m *m < CACHE_SIZE || m % 2 == 1){
        invu(mata,inv_mata,ar,ac,m);
    }else{
        int p = m / 2;
        woco_invu(mata, inv_mata,temp, ar, ac, p);
        woco_invu(mata, inv_mata,temp, ar+p, ac +p, p);
        init(temp,ar+m,ac+m);
        vex_wo(inv_mata, mata, temp, ar, ac, ar, ac + p, p, p ,p);
        matrix_sign_inv(temp, ar, ac + p, p);
        vex_wo(temp, inv_mata, inv_mata, ar, ac + p, ar + p, ac + p, p, p, p);
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
    float ** trsm_temp = NULL;
    square_matrix(&matl, m);
    square_matrix(&matu, m);
    square_matrix(&matp, m);
    square_matrix(&inv_matl,m);
    square_matrix(&inv_matu,m);
    square_matrix(&trsm_temp, m);
    init(matl,m);
    init(matu,m);
    init(matp,m);
    init(inv_matl,m);
    init(inv_matu,m);
    init(trsm_temp,m);
    woco_LUfac(mata, matl, matu, matp,trsm_temp, ar, ac, m, m ,m);
    woco_invl(matl, inv_matl, temp, ar, ac, m);
    init(temp, m);
    woco_invu(matu, inv_matu, temp, ar, ac, m);
    init(temp, m);
    vex_wo(inv_matu, inv_matl, temp, ar, ac, ar, ac, m, m, m);
    vex_wo(temp, matp, inv_mata, ar , ac, ar, ac, m, m, m);
    free(matl,m);
    free(matu,m);
    free(matp,m);
    free(inv_matu,m);
    free(inv_matl,m);
}

#endif //WONAL_INV_H

