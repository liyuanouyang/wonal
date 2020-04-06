//
// Created by oyliy on 2019/11/25.
//

#ifndef WONAL_TRSM_H
#define WONAL_TRSM_H


#include "basewo.h"

// matrix minus method
void matrix_minus(float **mata,float **matb,int ar,int ac,int br, int bc,int m,int n){
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            mata[ar+i][ac + j] -= matb[br+i][bc+j];
        }
    }
}


void direct_trsm(float **mata, float **matb, float **matx,
        int ar, int ac, int xr, int xc,
        int m, int n){
    //A is lower triangular matrix
    //A is m*m , B is m*n X is m*n ,wo find x
#pragma omp parallel for num_threads(core_num) default(none) shared(mata,matb,matx,m,n,ar,ac,xr,xc)
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            int temp = 0;
            for (int k = 0; k < i; ++k) {
                if(ar + i >= ac + k)
                    temp += mata[ar+i][ac+k]*matx[xr+k][xc+j];
            }
            matx[xr+i][xc+j] = (matb[ar+i][xc+j] - temp)/mata[ar+i][ac+i];
            temp = 0;
        }
    }
}


//need test
//A is m*m , B is m*n X is m*n ,wo find x
void woco_trsm(float **mata, float **matb, float **matx, float** temp,
        int ar, int ac,int xr, int xc, int m, int n){
    if (3*m*n < CACHE_SIZE || m % 2 == 1 || n % 2 == 1){
        direct_trsm(mata,matb,matx,ar,ac,xr,xc,m,n);
    } else{
        //partitiion to 4
        int r = m / 2;
        int c = n / 2;
        woco_trsm(mata,matb,matx,temp,ar,ac,xr,xc,r,c);
        woco_trsm(mata,matb,matx,temp,ar,ac,xr,xc + c,r,c);
        vex_wo(mata,matx,temp,ar+r,ac,xr,xc,r,r,c);
        //        if (flag == 1){vex_wo(mata,matx,temp,ar+r,ac,xr,xc,r,r,c);}
        //        else if(flag == 2){omp_vec_mat_mult(mata,matx,temp,ar+r,ac,xr,xc,r,r,c);}
        //        else if(flag == 3){matrix_multi(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
        //       matrix_minus(matb,temp,ar+q,xc,ar+q,xc,p);
        matrix_minus(matb,temp,xr+r,xc,ar+r,xc,r,c);
        woco_trsm(mata,matb,matx,temp,ar+r,ac+r,xr+r,xc,r,c);
        vex_wo(mata,matx,temp,ar+r,ac,xr,xc+c,r,r,c);
        //        if (flag ==1){vex_wo(mata,matx,temp,ar+r,ac,xr,xc+c,r,r,c);}
        //        else if (flag == 2){omp_vec_mat_mult(mata,matx,temp,ar+r,ac,xr,xc+c,r,r,c);}
        //        else if(flag == 3){matrix_multi(mata,matx,temp,ar+p,ac,xr,xc,p,p,p);}
        matrix_minus(matb,temp,xr+r,xc+c,ar+r,xc+c,r,c);
        woco_trsm(mata,matb,matx,temp,ar+r,ac+r,xr+r,xc+c,r,c);
    }
}
 //A is m*m , B is m*n X is m*n 
void trsm(float ** mata, float **matb, float **matx,int ar,int ac,int br,int bc,int m,int n){
    float ** temp = new_matrix(ar + m, bc + n);
    init(temp,ar + m, bc + n);
    woco_trsm(mata,matb,matx,temp,ar,ac,br,bc,m,n);
}

void test_trsm(){
    int m = 2048;
    int n = 2048;
    float ** mata = new_matrix(m,m);
    float ** matb = new_matrix(m,n);
    float ** matx = new_matrix(m,n);
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0; j < m ; ++j){
            if(i >= j) mata[i][j] = 1;
            else mata[i][j] = 0;
        }
    }
    for(int i = 0; i < m ; ++i){
        for (int j = 0; j < n; ++j){
            matb[i][j] = (i + 1)*(j + 1);
        }
    }
    init(matx,m,n);
    trsm(mata,matb,matx,0,0,0,0,m,n);
    display(matx,m,n);
    free(mata,m);
    free(matb,m);
    free(matx,m);
}

/*
void direct_trsm_u(float **mata, float **matb, float **matx,
        int ar, int ac, int xr, int xc,
        int m, int n){
    //A is upper triangular matrix
    //A is m*m , B is m*n X is m*n ,wo find x
    //square n
#pragma omp parallel for num_threads(core_num) default(none) shared(mata,matb,matx,m,n,ar,ac,xr,xc)
    for (int j = 0; j < n; ++j) {
        for (int i = m - 1; i >= 0; --i) {
            int temp = 0;
            for (int k = i; k < m ; ++k) {
                temp += mata[ar+i][ac+k]*matx[xr+k][xc+j];
            }
            matx[xr+i][xc+j] = (matb[ar+i][xc+j] - temp)/mata[ar+i][ac+i];
            temp = 0;
        }
    }
}

//need test
void woco_trsm_u(float **mata, float **matb, float **matx, float** temp,
        int ar, int ac,int xr, int xc, int m, int n){
    if (3*m*n < CACHE_SIZE||m % 2 == 1 || n % 2 == 1){
        direct_trsm_u(mata, matb, matx, ar, ac, xr, xc, m, n);
    } else{
        //partitiion to 4
        int r = m/2;
        int c = n/2;
        woco_trsm_u(mata, matb, matx, temp, ar + r, ac + r, xr + r, xc + c , r , c);
        woco_trsm_u(mata, matb, matx, temp, ar + r, ac + r, xr + r, xc, r, c);
        vex_wo(mata,matx,temp,ar,ac+r,xr+r,xc+c,r,r,c);
        //        if (flag == 1){vex_wo(mata,matx,temp,ar,ac+r,xr+r,xc+c,r,r,c);}
        //        else if(flag == 2){omp_vec_mat_mult(mata,matx,temp,ar,ac+r,xr+r,xc+c,r,r,c);}
        //        else if(flag == 3){matrix_multi(mata,matx,temp,ar,ac+r,xr+r,xc+c,r,r,c);}
        matrix_minus(matb,temp,xr,xc+c,ar,xc+c,r,c);
        woco_trsm_u(mata,matb,matx,temp,ar,ac,xr,xc+c,r,c);
        vex_wo(mata,matx,temp,ar,ac+r,xr+r,xc,r,r,c);
        //        if (flag ==1){vex_wo(mata,matx,temp,ar,ac+r,xr+r,xc,r,r,c);}
        //        else if (flag == 2){omp_vec_mat_mult(mata,matx,temp,ar,ac+r,xr+r,xc,r,r,c);}
        //        else if(flag == 3){matrix_multi(mata,matx,temp,ar,ac+r,xr+r,xc,r,r,c);}
        matrix_minus(matb,matx,xr,xc,ar,xc,r,c);
        woco_trsm_u(mata,matb,matx,temp,ar,ac,xr,xc,r,c);
    }
}

*/

#endif //WONAL_TRSM_H


