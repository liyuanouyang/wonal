/*
 * @Author       : whe
 * @Date         : 2020-04-07 20:18:01
 * @LastEditTime : 2020-04-07 20:53:25
 * @Description  : 
 */

#ifndef WONAL_TRMM_H
#define WONAL_TRMM_H

#include "basewo.h"

/*
A is lower triangular matrix
A * B = X, we find x.
A is m * m, B is m * n.
*/
void direct_trmm(float ** mata, float ** matb, float ** matx, int ar, int ac, int br, int bc, int m, int n){
    for(int i = 0; i < m; i++){
        for(int j = 0; j <= i; j++){
            for(int k = 0; k < n; k++)
                matx[ar+i][bc+k] += mata[ar+i][ac+j] * matb[br+j][bc+k];
        }
    }
}

void trmm(float **mata,float** matb,float ** matx,float ** temp,int ar,int ac,int br, int bc, int m, int n){
    if(3*m*n < CACHE_SIZE || m % 2 == 1 || n % 2 == 1){
        direct_trmm(mata,matb,matx,ar,ac,br,bc,m,n);
    }
    else{
        int r = m / 2;
        int c = n / 2;
        vex_wo(mata,matb,matx,ar,ac,br,bc,r,r,c);
        vex_wo(mata,matb,matx,ar,ac,br,bc+c,r,r,c);
        vex_wo(mata,matb,matx,ar+r,ac,br,bc,r,r,c);
        vex_wo(mata,matb,temp,ar+r,ac+r,br+r,bc,r,r,c);
        add(matx,temp,ar+r,bc,r,c,1,1);
        vex_wo(mata,matb,matx,ar+r,ac,br,bc+c,r,r,c);
        vex_wo(mata,matb,temp,ar+r,ac+r,br+r,bc+c,r,r,c);
        add(matx,temp,ar+r,bc+c,r,c,1,1);
    }
}
void test_trmm(){
    int m = 2048;
    int n = 2048;
    float ** mata = new_matrix(m,m);
    float ** matb = new_matrix(m,n);
    float ** matx = new_matrix(m,n);
    float ** temp = new_matrix(m,n);
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0; j < m ; ++j){
            if(i >= j) mata[i][j] = 1;
            else mata[i][j] = 0;
        }
    }
    for(int i = 0; i < m ; ++i){
        for (int j = 0; j < n; ++j){
            matb[i][j] = 2;
        }
    }
    init(matx,m,n);
    init(temp,m,n);
    direct_trmm(mata,matb,matx,0,0,0,0,m,n);
    //trmm(mata,matb,matx,temp,0,0,0,0,m,n);
    display(matx,m,n);
    free(mata,m);
    free(matb,m);
    free(matx,m);
    free(temp,m);
}
#endif