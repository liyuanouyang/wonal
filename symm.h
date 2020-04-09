/*
 * @Author       : whe
 * @Date         : 2020-04-08 11:51:26
 * @LastEditTime : 2020-04-08 19:42:22
 * @Description  : 
 */

#ifndef WONAL_SYMM_H
#define WONAL_SYMM_H

#include "basewo.h"

/*
A * B = X, we find x.
A is a symmetric matrix
A is m * m, B is m * n.
*/
// need optimize
void direct_symm(float ** mata, float ** matb, float ** matx, int ar, int ac, int br, int bc, int m, int n){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < m; j++){
            for(int k = 0; k < n; k++)
                matx[ar+i][bc+k] += mata[ar+i][ac+j] * matb[br+j][bc+k];
        }
    }
}

void symm(float **mata,float** matb,float ** matx,float ** temp,int ar,int ac,int br, int bc, int m, int n){
    if(3*m*n <  CACHE_SIZE|| m % 2 == 1 || n % 2 == 1){
        direct_symm(mata,matb,matx,ar,ac,br,bc,m,n);
    }
    else{
        int r = m / 2;
        int c = n / 2;
        direct_symm(mata,matb,matx,ar,ac,br,bc,r,c);
        vex_wo(mata,matb,temp,ar,ac+r,br+r,bc,r,r,c);
        

        direct_symm(mata,matb,matx,ar,ac,br,bc+c,r,c);
        vex_wo(mata,matb,temp,ar,ac+r,br+r,bc+c,r,r,c);

        direct_symm(mata,matb,matx,ar+r,ac+r,br+r,bc,r,c);
        vex_wo(mata,matb,temp,ar+r,ac,br,bc,r,r,c);

        direct_symm(mata,matb,matx,ar+r,ac+r,br+r,bc+c,r,c);
        vex_wo(mata,matb,matx,ar+r,ac,br,bc+c,r,r,c);

        add(matx,temp,ar,bc,m,n,1,1);
    }
}
void test_symm(){
    int m = 2048;
    int n = 2048;
    float ** mata = new_matrix(m,m);
    float ** matb = new_matrix(m,n);
    float ** matx = new_matrix(m,n);
    float ** temp = new_matrix(m,n);
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0; j < m ; ++j){
            if(i >= j) mata[i][j] = 1;
            else mata[i][j] = 1;
        }
    }
    for(int i = 0; i < m ; ++i){
        for (int j = 0; j < n; ++j){
            matb[i][j] = 2;
        }
    }
    init(matx,m,n);
    init(temp,m,n);
    //direct_symm(mata,matb,matx,0,0,0,0,m,n);
    symm(mata,matb,matx,temp,0,0,0,0,m,n);
    display(matx,m,n);
    free(mata,m);
    free(matb,m);
    free(matx,m);
    free(temp,m);
}
#endif