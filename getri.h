/*
 * @Author       : whe
 * @Date         : 2020-04-18 11:32:42
 * @LastEditTime : 2020-04-18 12:03:30
 * @Description  : 
 */
#ifndef WONAL_GETRI_H
#define WONAL_GETRI_H

#include "basewo.h"
#include "inv.h"
#include "LUfac.h"
/*
GETRI computes the inverse of a real symmetric positive definite
 matrix A using the LU factorization AP = LU
 */
void getri(float ** mata, float ** inv_mata, int ar, int ac, int m){
    float ** matl = NULL;
    float ** matu = NULL;
    float ** matp = NULL;
    float ** inv_matl = NULL;
    float ** inv_matu = NULL;
    float ** trsm_temp = NULL;
    float ** temp = NULL;
    square_matrix(&matl, m);
    square_matrix(&matu, m);
    square_matrix(&matp, m);
    square_matrix(&inv_matl,m);
    square_matrix(&inv_matu,m);
    square_matrix(&trsm_temp, m);
    square_matrix(&temp,m);
    init(matl,m);
    init(matu,m);
    init(matp,m);
    init(inv_matl,m);
    init(inv_matu,m);
    init(trsm_temp,m);
    init(temp,m);
    woco_LUfac(mata, matl, matu, matp,trsm_temp, ar, ac, m, m ,m);
    woco_invl(matl, inv_matl, temp, ar, ac, m);
    init(temp, m);
    woco_invu(matu, inv_matu, temp, ar, ac, m);
    init(temp, m);
    vex_wo(inv_matu, inv_matl, temp, ar, ac, ar, ac, m, m, m);
    vex_wo(matp, temp, inv_mata, ar , ac, ar, ac, m, m, m);
    free(matl,m);
    free(matu,m);
    free(matp,m);
    free(inv_matu,m);
    free(inv_matl,m);
}

void test_getri(){
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
    getri(mata,matx,0,0,m);
    display(matx,m,m);
    free(mata,m);
    free(matx,m);
}
#endif