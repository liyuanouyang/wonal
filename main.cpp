/* 测试代码用文件 */
#include "basewo.h"
//#include "trsm.h"
//#include "LUfac.h"
//#include "Cholesky_fac.h"


using namespace std;

//1024 2048 4096 6144 8192 10240 12288 14336 16384

//#define OMIGA 1
#define OMIGA 1

int INPUT_M = 2048;// size of matrix
int INPUT_N = 2048;
int INPUT_K = 2048;




//test mm
int main(int argc, char *argv[]){
    int i = 0;
//openmp 配置
//    core_num = omp_get_max_threads();
    core_num = atoi(argv[1]);
    cout<<core_num<<endl;

//memkind配置
    struct memkind *pmem_kind_unlimited = NULL;
    int err = 0;

    err = memkind_create_pmem("/mnt/persist-memory/pmem10/", 0, &pmem_kind_unlimited);
    if (err) {
        print_err_message(err);
        return 1;
    }
//初始化矩阵
    float **matrix_a = NULL;
    float **matrix_c = NULL;
    float **matrix_b = NULL;

    matrix_a = (float **)memkind_malloc(pmem_kind_unlimited, 512*(INPUT_M * sizeof(float*)/512+1));
    for (int i = 0; i < INPUT_M; i++) {
        matrix_a[i] = (float *)memkind_malloc(pmem_kind_unlimited, 512*(INPUT_K * sizeof(float)/512 +1));
    }
    matrix_b = (float **)memkind_malloc(pmem_kind_unlimited, 512*(INPUT_K * sizeof(float*)/512+1));
    for (int i = 0; i < INPUT_K; ++i) {
        matrix_b[i] = (float *)memkind_malloc(pmem_kind_unlimited, 512*(INPUT_N * sizeof(float)/512+1));
    }

    matrix_c = (float **)memkind_malloc(pmem_kind_unlimited, 512*(INPUT_M * sizeof(float*)/512+1));
    for (int i = 0; i < INPUT_M; ++i) {
        matrix_c[i] = (float *)memkind_malloc(pmem_kind_unlimited, 512*(INPUT_N * sizeof(float)/512+1));
    }

    for (int k = 0; k < INPUT_M; ++k) {
        for (int i = 0; i < INPUT_K; ++i) {
            matrix_a[k][i] = 1.0;
        }
    }

    for (int i = 0; i < INPUT_M; ++i) {
        for (int j = 0; j < INPUT_N; ++j) {
            matrix_c[i][j] = 0.0;
        }
    }

    for (int i = 0; i <INPUT_K ; ++i) {
        for (int j = 0; j < INPUT_N; ++j) {
            matrix_b[i][j] = 2.0;
        }
    }
//record time
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    int bran = atoi(argv[2]);
    cout << bran << endl;
    //分支,普通矩阵乘法和优化后的对比
    if (bran == 2){
        vex_wo(matrix_a,matrix_b,matrix_c,0,0,0,0,INPUT_M,INPUT_K,INPUT_N);
    }else if(bran == 4){
        omp_vec_mat_mult(matrix_a,matrix_b,matrix_c,0,0,0,0,INPUT_M,INPUT_K,INPUT_N);
    }
    else{
        cout<< "argv2 fault"<< endl;
    }

//    matrix_multi(matrix_a,matrix_b,matrix_c,0,0,0,0,INPUT_M,INPUT_K,INPUT_N);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

//输出
    for (int k = 0; k < 16; ++k) {
        for (int i = 0; i < 16; ++i) {
            std::cout<<matrix_c[k][i]<<std::endl;
            std::cout<<matrix_c[INPUT_M-k][INPUT_M-i]<<std::endl;
            std::cout<<matrix_c[k][INPUT_M-i]<<std::endl;
            std::cout<<matrix_c[INPUT_M - k][i]<<std::endl;
        }
    }

    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout<<"TIme"<<time_span.count()<<"ms"<<std::endl;
    return 0;
}
