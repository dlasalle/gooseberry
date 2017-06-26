/**
 * @file cgd.c
 * @brief Conjugate gradient functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-05-06
 */




#ifndef GOOSEBERRY_CGD_C
#define GOOSEBERRY_CGD_C




#include "cgd.h"
#include "blas.h"




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


static void __print_vector(const char * const name, matrix_t * v)
{
  dim_t i;
  printf("%s: [",name);
  for (i=0;i<v->nrows;++i) {
    printf(PF_REAL_T",",v->rowval[i]);
  }
  printf("]\n");
}

int cgd(
    matrix_t * const a, 
    matrix_t * const b,
    matrix_t * const x,
    const real_t error,
    const size_t niter)
{
  size_t iter;
  real_t alpha, beta, rmse, dr, dnr;
  matrix_t * r, * p, * nr, * nx, *ap;

  if (a->nrows != a->ncols) {
    eprintf("Input matrix for CGD is not square: "PF_DIM_T"x"PF_DIM_T"\n",
        a->nrows,a->ncols);
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (a->nrows != b->nrows) {
    eprintf("Input vector for CGD has wrong number of rows: "PF_DIM_T"\n",
        b->nrows);
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (a->ncols != x->nrows) {
    eprintf("Output vector for CGD has wrong number of rows: "PF_DIM_T"\n",
        x->nrows);
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  
  r = matrix_alloc(1);
  p = matrix_alloc(1);
  nr = matrix_alloc(1);
  nx = matrix_alloc(1);
  ap = matrix_alloc(1);

  matrix_init(MATRIX_TYPE_DENSE_VECTOR,a->nrows,1,0,r);
  matrix_init(MATRIX_TYPE_DENSE_VECTOR,a->nrows,1,0,p);
  matrix_init(MATRIX_TYPE_DENSE_VECTOR,a->nrows,1,0,nr);
  matrix_init(MATRIX_TYPE_DENSE_VECTOR,a->ncols,1,0,nx);
  matrix_init(MATRIX_TYPE_DENSE_VECTOR,a->ncols,1,0,ap);


  /* Assume x_0 is all zeroes : r_0 = b - Ax_0 */
  real_set(x->rowval,0,x->nrows);
  real_copy(r->rowval,b->rowval,b->nrows);
  /* set p_0 = r_0 */
  real_copy(p->rowval,r->rowval,r->nrows);
  dr = blas_dot(r->nrows,r->rowval,r->rowval);

  for (iter=0;iter<niter||niter==0;++iter) {
    blas_spmult(a->nrows,a->ncols,p->ncols,a->rowptr,a->rowind,a->rowval,
        p->rowval,ap->rowval,a->rdist,a->nrdist);
    alpha = dr / blas_dot(p->nrows,p->rowval,ap->rowval);
    /* nx = x + alpha p */
    blas_scale(p->nrows,p->rowval,alpha);
    blas_add(x->nrows,x->ncols,x->rowval,p->rowval,nx->rowval);

    /* nr = r - alpha a p */
    blas_scale(ap->nrows,ap->rowval,alpha);
    blas_sub(r->nrows,r->ncols,r->rowval,ap->rowval,nr->rowval);

    /* calculate beta */
    dnr = blas_dot(nr->nrows,nr->rowval,nr->rowval);
    beta = dnr / dr;

    /* set p_{i+1} = r + beta p_{i} */
    blas_scale(p->nrows,p->rowval,beta/alpha);     
    blas_add(nr->nrows,nr->ncols,nr->rowval,p->rowval,p->rowval);

    /* update vectors */
    real_copy(x->rowval,nx->rowval,x->nrows);
    real_copy(r->rowval,nr->rowval,r->nrows);

    /* update r dot product */
    dr = dnr;

    /* decide if r is small enough to exit */
    rmse = sqrt(dr);
    printf("ITER: "PF_SIZE_T" RMSE: "PF_REAL_T"\n",iter,rmse);
    if (rmse <= error) {
      break;
    }
  }

  matrix_free(r);
  matrix_free(p);
  matrix_free(nr);
  matrix_free(nx);
  matrix_free(ap);

  return GOOSEBERRY_SUCCESS;
}




#endif
