//
// Created by oyliy on 2019/12/12.
//

#ifndef WONAL_QRFAC_H
#define WONAL_QRFAC_H
//#define mat float**

#include "basewo.h"
#include <math.h>


//tranpose on itself
void matrix_tranpose(float** mata,int m,int n){
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            float t = mata[i][j];
            mata[i][j] = mata[j][i];
            mata[j][i] = t;
        }
    }
}

float **matrix_minor(float** mata,int m,int n,int d){
    float** matm = new_matrix(m,n);
    for (int i = 0; i < d; ++i) {
        matm[i][i] = 1;
    }
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            matm[i][j] = mata[i][j];
        }
    }
    return matm;
}
//c = a+b*s
float* vmadd(float a[],float b[], float s, float c[], int n){
    for (int i = 0; i < n; ++i) {
        c[i] = a[i] + s*b[i];
    }
    return c;
}
// m = I-v v^T
float** vmul(float v[] , int n){
    float** matx = new_matrix(n,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matx[i][j] = -2*v[i]*v[j];
        }
    }
    for (int i = 0; i < n; ++i) {
        matx[i][i] +=1;
    }
    return matx;
}
/* ||x||*/
float vnorm(float x[], int n){
    float sum = 0;
    for (int i = 0; i < n; ++i) {
        sum+= x[i]*x[i];
    }
    return sqrt(sum);
}
//y = x/d
float* vdiv(float x[], float d, float y[] ,int n){
    for (int i = 0; i < n; ++i) {
        y[i] = x[i]/d;
    }
    return y;
}

//take a col of matm ,put in v
float* mcol(float **matm,int m, float *v,int c){
    for (int i = 0; i < m; ++i) {
        v[i] = matm[i][c];
    }
}
//m,n indicate sieze of matm
void householder(float** matm, float** *R,float** *Q,int m,int lenk, int n,int err, struct memkind *pmem_kind_unlimited){
    float** q[m];
    float **z = matm,**z1;
    for (int k = 0; k < n && k< m-1; ++k) {
        float e[m],x[m],a;
        z1 = matrix_minor(z,m,n,k);
        z = z1;
        mcol(z,m,x,k);
        a = vnorm(x,m);
        if (matm[k][k] > 0)    a = -a;
        for (int i = 0; i < m; ++i) {
            e[i] = (i==k)?1:0;
        }
        vmadd(x,e,a,e,m);
        vdiv(e,vnorm(e,m),e,m);
        q[k] = vmul(e,m);
        z1 = matrix_mul(q[k],z,m,lenk,n,err,&pmem_kind_unlimited);
        if (z != matm){
            free(z,m);
        }
        z = z1;
    }
    free(z,m);
    *Q = q[0];
    *R = matrix_mul(q[0],matm,m,lenk,n,err,&pmem_kind_unlimited);
    for (int i = 1; i < n&& i<m; ++i) {
        z1 = matrix_mul(q[i],*Q,m,lenk,n,err,&pmem_kind_unlimited);
        if (i>1)    free(*Q);
        *Q = z1;
        free(q[i],m);
    }
    free(q[0]);
    z = matrix_mul(*q,matm,m,lenk,n,err,&pmem_kind_unlimited);
    free(*R);
    *R = z;
    matrix_tranpose(*Q,m,lenk);
}


#endif //WONAL_QRFAC_H
