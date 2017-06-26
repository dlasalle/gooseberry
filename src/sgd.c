/**
 * @file sgd.c
 * @brief Functions for Stocastic Gradient Descent
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-04-26
 */




#ifndef GOOSEBERRY_SGD_C
#define GOOSEBERRY_SGD_C




#include "sgd.h"




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int sgd(
    matrix_t * const mat, 
    matrix_t * const u, 
    matrix_t * const v,
    const real_t rate,
    const real_t lambda,
    const size_t niter)
{
  int rv = GOOSEBERRY_SUCCESS;
  size_t iter;
  dim_t i, row, col, k;
  ind_t j, p;
  real_t r;
  ind_t * perm = NULL;
  dim_t * rowid = NULL;
  real_t * wu, * wv;

  const dim_t nrows = mat->nrows;
  const dim_t ncols = mat->ncols;
  const ind_t nnz = mat->rowptr[nrows];
  const ind_t * const ptr = mat->rowptr;
  const dim_t * const ind = mat->rowind;
  const real_t * const val = mat->rowval;
  const dim_t nfactors = u->ncols;

  if (u->nrows != nrows) {
    eprintf("The input matrix U has an invalid number of rows "PF_DIM_T"\n",
        u->nrows);
    rv = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  }
  if (v->nrows != ncols) {
    eprintf("The input matrix V has an invalid number of rows "PF_DIM_T"\n",
        v->nrows);
    rv = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  }
  if (v->ncols != u->ncols) {
    eprintf("U an V have different numbers of columns: "PF_DIM_T" vs "PF_DIM_T
        "\n",u->ncols,v->ncols);
    rv = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  }

  perm = ind_alloc(nnz);
  rowid = dim_alloc(nnz);

  /* build a row lookup */
  for (i=0;i<mat->nrows;++i) {
    for (j=ptr[i];j<ptr[i+1];++j) {
      rowid[j] = i;
    }
  }

  for (iter=0;iter<niter;++iter) {
    ind_incset(perm,0,1,nnz);
    ind_pseudo_shuffle(perm,nnz/8,nnz);

    for (p=0;p<nnz;++p) {
      j = perm[p];

      row = rowid[j];
      col = ind[j];

      wu = u->rowval+(row*ncols);
      wv = v->rowval+(col*nrows);

      /* calculate error */
      r = val[j];
      for (k=0;k<nfactors;++k) {
        r -= wu[k]*wv[k];
      }

      /* update factors */
      for (k=0;k<nfactors;++k) {
        wv[k] += rate*((r*wu[k]) - lambda*wv[k]);
      }
      for (k=0;k<nfactors;++k) {
        wu[k] += rate*((r*wv[k]) - lambda*wu[k]);
      }
    }
    
  }

  END:

  if (perm) {
    dl_free(perm);
  }
  if (rowid) {
    dl_free(rowid);
  }

  return rv;
}




#endif
