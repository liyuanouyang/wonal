//
// Created by oyliy on 2019/11/27.
//

//12.26 try to change b to row-main matrix
#ifndef WONAL_BASEWO_H
#define WONAL_BASEWO_H
#include <iostream>
#include <chrono>
#include <omp.h>
#include <math.h>
#include <memkind.h>
#include <immintrin.h>

#define CACHE_SIZE 7471104

static int core_num;
int flag = 3;
//int truu = 1;


int max(int x,int y){
    if (x>y){
        return x;
    } else{
        return y;
    }
}
//输出异常函数
static void print_err_message(int err)
{
    char error_message[MEMKIND_ERROR_MESSAGE_SIZE];
    memkind_error_message(err, error_message, MEMKIND_ERROR_MESSAGE_SIZE);
    fprintf(stderr, "%s\n", error_message);
}
//矩阵乘法，矩阵abc，a的起始位置(ar,ac)，b的起始，mkn表示相乘的大小，即从(ar,ac）开始大小为m,k的子矩阵和b的k,n子矩阵相乘
void matrix_multi(
        float **mata,float **matb,float **matc,
        int ar, int ac,int br,int bc, int m, int k ,int n){
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int l = 0;
            while (l<k){
                //原子化
#pragma omp atomic
                matc[ar + i][bc + j] += mata[ar + i][ac + l] * matb[br + l][bc + j];
                l++;
            }
        }
    }
}
//参数同上，加入了avx512向量化
void avx_mult(float **mata,float **matb,float **matc,
              int ar, int ac,int br,int bc, int m, int k ,int n){
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            __m512 sumx16 = _mm512_set1_ps(0.0);
            float sum = 0;
            int l = 0;
            while ( l < k) {
                __m512 a = _mm512_loadu_ps(&(mata[ar + i][ac + l]));
                __m512 b = _mm512_loadu_ps(&(matb[br + l][bc + j]));
                sumx16 = _mm512_fmadd_ps(a,b,sumx16);
                l+=16;
            }
            sum += _mm512_reduce_add_ps(sumx16);
#pragma omp atomic
            matc[ar + i][bc + j] += sum;
        }
    }
}


//bug : the last divisio doesnt need 16
//并行加向量化
void omp_vec_mat_mult(float **mata,float **matb,float **matc,
                      int ar, int ac,int br,int bc, int m, int k ,int n){
#pragma omp parallel for num_threads(core_num) default(none) shared(m,n,k,mata,matb,matc,ar,ac,br,bc)
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            __m512 sumx16 = _mm512_set1_ps(0.0);
            float sum = 0;
            int l = 0;
            while( l < k) {
                __m512 a = _mm512_loadu_ps(&(mata[ar + i][ac + l]));
                __m512 b = _mm512_loadu_ps(&(matb[br + l][bc + j]));
                sumx16 = _mm512_fmadd_ps(a,b,sumx16);
                l+=16;
            }
            sum += _mm512_reduce_add_ps(sumx16);
#pragma omp atomic
            matc[ar + i][bc + j] += sum;
        }
    }
}

//带返回值的函数
float** omp_vec_mat_mult_return(float **mata,float **matb,
                      int ar, int ac,int br,int bc, int m, int k ,int n ,int inputm,int inputn){
    float** matc = new float*[inputm];
    for (int i = 0; i < m; ++i) {
        matc[i] = new float[inputn];
    }
#pragma omp parallel for num_threads(core_num) default(none) shared(m,n,k,mata,matb,matc,ar,ac,br,bc)
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            __m512 sumx16 = _mm512_set1_ps(0.0);
            float sum = 0;
            int l = 0;
            while( l < k) {
                __m512 a = _mm512_loadu_ps(&(mata[ar + i][ac + l]));
                __m512 b = _mm512_loadu_ps(&(matb[br + l][bc + j]));
                sumx16 = _mm512_fmadd_ps(a,b,sumx16);
                l+=16;
            }
            sum += _mm512_reduce_add_ps(sumx16);
#pragma omp atomic
            matc[ar+i][bc+j] += sum;
        }
    }
    return matc;
}
//wo是我们论文中优化后的算法，是一个矩阵划分方式
void vex_wo(
        float **mata, float **matb, float **matc,
        int ar, int ac,int br, int bc, int m, int k, int n){
    if (m*n < CACHE_SIZE/15){
        //divide k
        int sec_len = k/core_num;
#pragma omp parallel for num_threads(core_num) default(none) shared(core_num,sec_len,mata,matb,matc,ar,ac,br,bc,m,n,k)
        for (int i = 0; i < core_num; ++i) {
            int start_point = i*sec_len;
            int sec = i != core_num - 1? sec_len: k- i*sec_len;
            avx_mult(mata,matb,matc,ar,ac+start_point,br+start_point,bc,m,sec,n);
        }
    }else if (m>n){
        //cut m
        int p = m/2;
        vex_wo(mata,matb,matc,ar,ac,br,bc,p,k,n);
        vex_wo(mata,matb,matc,ar+p,ac,br,bc,m-p,k,n);
    } else{
        //cut n
        int p = n/2;
        vex_wo(mata,matb,matc,ar,ac,br,bc,m,k,p);
        vex_wo(mata,matb,matc,ar,ac,br,bc+p,m,k,n-p);
    }
}


//void vex_wo(
//        float **mata, float **matb, float **matc,
//        int ar, int ac,int br, int bc, int m, int k, int n){
//    if (m*n < CACHE_SIZE/15){
//        //divide k
//        int sec_len = k/core_num;
////#pragma omp parallel for num_threads(core_num) default(none) shared(core_num,sec_len,mata,matb,matc,ar,ac,br,bc,m,n,k)
//        for (int i = 0; i < core_num; ++i) {
//            int start_point = i*sec_len;
//            int sec = i==core_num-1? k- i*sec_len : sec_len;
////            if (truu){
////                printf("%d,%d,%d\n",start_point,sec,i);
////            }
//            avx_mult(mata,matb,matc,ar,ac+start_point,br+start_point,bc,m,sec,n);
//        }
////        truu = 0;
//    }else if (m>n){
//        //cut m
//        int p = m/2;
//        vex_wo(mata,matb,matc,ar,ac,br,bc,p,k,n);
//        vex_wo(mata,matb,matc,ar+p,ac,br,bc,m-p,k,n);
//    } else{
//        //cut n
//        int p = n/2;
//        vex_wo(mata,matb,matc,ar,ac,br,bc,m,k,p);
//        vex_wo(mata,matb,matc,ar,ac,br,bc+p,m,k,n-p);
//    }
//}

//float** vex_wo_return(
//        float **mata, float **matb,
//        int ar, int ac,int br, int bc, int m, int k, int n){
//    float ** matc = new float*[m];
//    for (int i = 0; i < m; ++i) {
//        matc[i] = new float[n];
//    }
//    if (m*n < CACHE_SIZE/15){
//        //divide k
//        int sec_len = k/core_num;
//#pragma omp parallel for num_threads(core_num) default(none) shared(core_num,sec_len,mata,matb,matc,ar,ac,br,bc,m,n,k)
//        for (int i = 0; i < core_num; ++i) {
//            int start_point = i*sec_len;
//            int sec = i==core_num-1? k- i*sec_len : sec_len;
////            if (truu){
////                printf("%d,%d,%d\n",start_point,sec,i);
////            }
//            avx_mult(mata,matb,matc,ar,ac+start_point,br+start_point,bc,m,sec,n);
//        }
////        truu = 0;
//    }else if (m>n){
//        //cut m
//        int p = m/2;
//        vex_wo(mata,matb,matc,ar,ac,br,bc,p,k,n);
//        vex_wo(mata,matb,matc,ar+p,ac,br,bc,m-p,k,n);
//    } else{
//        //cut n
//        int p = n/2;
//        vex_wo(mata,matb,matc,ar,ac,br,bc,m,k,p);
//        vex_wo(mata,matb,matc,ar,ac,br,bc+p,m,k,n-p);
//    }
//    return matc;
//}

//转置
void transport(float** mata, float** matb, int m, int n){
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            matb[j][i] = mata[i][j];
        }
    }
}
//归零
void  init(float** matrix , int n){
    for (int i=0 ;i < n; i++){
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = 0;
        }
    }
}

void  init(float** matrix , int xr, int xc, int n){
    for (int i=0 ;i < n; i++){
        for (int j = 0; j < n; ++j) {
            matrix[xr+i][xc+j] = 0.0;
        }
    }
}

void  init(float** matrix ,int m, int n){
    for (int i=0 ;i < m; i++){
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = 0;
        }
    }
}

void  init(float** matrix , int xr, int xc, int m, int n){
    for (int i=0 ;i < m; i++){
        for (int j = 0; j < n; ++j) {
            matrix[xr+i][xc+j] = 0.0;
        }
    }
}

void free(float** matrix , int m){
    for (int i = 0;i<m;i++){
        free(matrix[i]);
    }
}

void square_matrix(float** *matrix,int n){
    *matrix = new float *[n];
    for (int i = 0; i < n; ++i) {
        (*matrix)[i] = new float[n];
    }
}

float** new_matrix(int m,int n){
    float **matrix = new float *[n];
    for (int i = 0; i < n; ++i) {
        matrix[i] = new float[n];
    }
    return matrix;
}


//&or*          need to pass by &
float** matrix_mul(float** mata,float** matb,int m,int k,int n,int err ,struct memkind **pmem_kind_unlimited){
    float** matrix_c = NULL;
    matrix_c = (float** )memkind_malloc(*pmem_kind_unlimited, m * sizeof(float*));
    for (int i = 0; i < m; ++i) {
        matrix_c[i] = (float *)memkind_malloc(*pmem_kind_unlimited, n * sizeof(float));
    }
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int l = 0;
            while (l<k){
#pragma omp atomic
                matrix_c[i][j] += mata[i][l] * matb[j][l];
                l++;
            }
        }
    }
}

//void write_optimal(
//        float **mata, float **matb, float **matc,
//        int ar, int ac,int br, int bc, int m, int k, int n){
//
//    if (m*k+k*n+m*n<CACHE_SIZE){
//        matrix_multi(mata,matb,matc,ar,ac,br,bc,m,k,n);
//    } else if (m*n < CACHE_SIZE){
//        //divide k
//        int p = k/2;
//        write_optimal(mata,matb,matc,ar,ac,br,bc,m,p,n);
//        write_optimal(mata,matb,matc,ar,ac+p,br+p,bc,m,k-p,n);
//    } else if (m>n){
//        //cut m
//        int p = m/2;
//        write_optimal(mata,matb,matc,ar,ac,br,bc,p,k,n);
//        write_optimal(mata,matb,matc,ar+p,ac,br,bc,m-p,k,n);
//    } else{
//        //cut n
//        int p = n/2;
//        write_optimal(mata,matb,matc,ar,ac,br,bc,m,k,p);
//        write_optimal(mata,matb,matc,ar,ac,br,bc+p,m,k,n-p);
//    }
//}



//float** write_optimal(float **mata, float **matb,
//                      int ar, int ac,int br, int bc, int m, int k, int n){
//    float **matc = new float *[m];
//    for (int i = 0; i < n; ++i) {
//        matc[i] = new float [n];
//    }
//    if (m*k+k*n+m*n<CACHE_SIZE){
//        matrix_multi(mata,matb,matc,ar,ac,br,bc,m,k,n);
//    } else if (m*n < CACHE_SIZE){
//        //divide k
//        int p = k/2;
//        write_optimal(mata,matb,matc,ar,ac,br,bc,m,p,n);
//        write_optimal(mata,matb,matc,ar,ac+p,br+p,bc,m,p,n);
//    } else if (m>n){
//        //cut m
//        int p = m/2;
//        write_optimal(mata,matb,matc,ar,ac,br,bc,p,k,n);
//        write_optimal(mata,matb,matc,ar+p,ac,br,bc,p,k,n);
//    } else{
//        //cut n
//        int p = n/2;
//        write_optimal(mata,matb,matc,ar,ac,br,bc,m,k,p);
//        write_optimal(mata,matb,matc,ar,ac,br,bc+p,m,k,p);
//    }
//}

#endif //WONAL_BASEWO_H
