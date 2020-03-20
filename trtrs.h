#ifndef WONAL_TRTRS_H
#define WONAL_TRTRS_H

#include "trsm.h"

//Purpose:
//solves a triangular system of the form
//A * X = B or A**T * X = B,
//A is a triangular matrix of order M,and B is and M by N matrix.
void trtrs(char UPLO, char TRANS,
           float **mata, float **matb, float **matx,
           int ar, int ac, int xr, int xc,
           int m, int n){
    trsm('L',UPLO,TRANS,mata,matb,matx,ar,ac,xr,xc,m,n);

}

#endif //WONAL_TRTRS_H

//reference : http://www.netlib.org/lapack/explore-html/da/dba/group__double_o_t_h_e_rcomputational_ga7068947990361e55177155d044435a5c.html#ga7068947990361e55177155d044435a5c