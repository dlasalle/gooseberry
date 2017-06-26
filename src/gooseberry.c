/**
 * @file gooseberry.c
 * @brief Library API functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-05-01
 */




#ifndef GOOSEBERRY_C
#define GOOSEBERRY_C



#include "base.h"
#include "io.h"
#include "matrix.h"
#include "sgd.h"
#include "blas.h"
#include "transform.h"



/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


/* BLAS **********************************************************************/


int gooseberry_scale(
    const ind_t n,
    real_t * const val,
    const real_t s)
{
  return blas_scale(n,val,s);
}


int gooseberry_add_scalar(
    const ind_t n,
    real_t * const val,
    const real_t s)
{
  return blas_add_scalar(n,val,s);
}


int gooseberry_spmult(
    const dim_t anrows,
    const dim_t ancols,
    const dim_t bncols,
    const ind_t * const arowptr,
    const dim_t * const arowind,
    const real_t * const arowval,
    const real_t * const bcolval,
    real_t * const crowval,
    const dim_t * const blocks,
    const dim_t nblocks)
{
  return blas_spmult(anrows,ancols,bncols,arowptr,arowind,arowval,bcolval,
      crowval,blocks,nblocks);
}


int gooseberry_spmultsp(
    const dim_t anrows,
    const dim_t ancols,
    const dim_t bncols,
    const ind_t * const arowptr,
    const dim_t * const arowind,
    const real_t * const arowval,
    const ind_t * const browptr,
    const dim_t * const browind,
    const real_t * const browval,
    real_t * const crowval,
    const dim_t * const blocks,
    const dim_t nblocks)
{
  return blas_spmultsp(anrows,ancols,bncols,arowptr,arowind,arowval,browptr,
      browind,browval,crowval,blocks,nblocks);
}


int gooseberry_spmultsp_sp(
    const dim_t anrows,
    const dim_t ancols,
    const dim_t bncols,
    const ind_t * const arowptr,
    const dim_t * const arowind,
    const real_t * const arowval,
    const ind_t * const browptr,
    const dim_t * const browind,
    const real_t * const browval,
    ind_t ** const crowptr,
    dim_t ** const crowind,
    real_t ** const crowval,
    const dim_t * const blocks,
    const dim_t nblocks)
{
  return blas_spmultsp_sp(anrows,ancols,bncols,arowptr,arowind,arowval,browptr,
      browind,browval,crowptr,crowind,crowval,blocks,nblocks);
}



/* TRANSFORM *****************************************************************/

int gooseberry_symmetrify_sparse(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  return transform_symmetrify_sparse(nrows,ncols,rowptr,rowind,rowval,
      r_rowptr,r_rowind,r_rowval);
}


int gooseberry_debipartify_sparse(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  return transform_debipartify_sparse(nrows,ncols,rowptr,rowind,rowval,
      r_rowptr,r_rowind,r_rowval);
}


int gooseberry_rowsplit_sparse(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    const dim_t nparts,
    const dim_t * const dist,
    const dim_t * const map,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  return transform_rowsplit_sparse(nrows,ncols,rowptr,rowind,rowval,
      nparts,dist,map,r_rowptr,r_rowind,r_rowval);
}


int gooseberry_colsplit_sparse(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    const dim_t nparts,
    const dim_t * const dist,
    const dim_t * const map,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  return transform_colsplit_sparse(nrows,ncols,rowptr,rowind,rowval,
      nparts,dist,map,r_rowptr,r_rowind,r_rowval);
}


int gooseberry_transpose_sparse(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    ind_t ** const r_trowptr,
    dim_t ** const r_trowind,
    real_t ** const r_trowval)
{
  return transform_transpose_sparse(nrows,ncols,rowptr,rowind,rowval,
      r_trowptr,r_trowind,r_trowval);
}


/* IO ************************************************************************/

int gooseberry_read_sparse_matrix(
    const int type, 
    const char * const file, 
    dim_t * const nrows,
    dim_t * const ncols,
    ind_t ** const rowptr,
    dim_t ** const rowind,
    real_t ** const rowval)
{
  return read_sparse_matrix(type,file,nrows,ncols,rowptr,rowind,rowval);
}


int gooseberry_read_dense_matrix(
    int type, 
    const char * const file, 
    dim_t * const nrows,
    dim_t * const ncols,
    real_t ** const rowval)
{
  return read_dense_matrix(type,file,nrows,ncols,rowval);
}


int gooseberry_write_dense_matrix(
    const int type, 
    const char * const file, 
    const dim_t nrows,
    const dim_t ncols,
    const real_t * const rowval)
{
  return write_dense_matrix(type,file,nrows,ncols,rowval);
}


int gooseberry_write_sparse_matrix(
    const int type, 
    const char * const file, 
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval)
{
  return write_sparse_matrix(type,file,nrows,ncols,rowptr,rowind,rowval);
}


int gooseberry_read_labels(
    const char * const file,
    dim_t * const nrows,
    dim_t ** const labels)
{
  return read_labels(file,nrows,labels);
}


int gooseberry_read_partition(
    const char * const file,
    dim_t * const r_nrows,
    dim_t * const r_nparts,
    dim_t ** const r_map,
    dim_t ** const r_perm,
    dim_t ** const r_dist)
{
  dim_t nrows, nparts;
  dim_t * perm, * order, * dist, * map;

  if (r_nrows) {
    nrows = *r_nrows;
  }

  read_labels(file,&nrows,&map);

  nparts = dim_max_value(map,nrows)+1;

  perm = dim_alloc(nrows);
  order = dim_alloc(nrows);
  dim_incset(order,0,1,nrows);
  dd_countingsort_kv(map,order,0,nparts-1,nrows,perm,&dist);
  dl_free(order);

  if (r_nrows) {
    /* hanldes the case where *r_nrows was passed in pointing to null */
    *r_nrows = nrows;
  }

  if (r_nparts) {
    *r_nparts = nparts;
  }

  if (r_dist) {
    *r_dist = dist;
  } else {
    dl_free(dist);
  }

  if (r_perm) {
    *r_perm = perm;
  } else {
    dl_free(perm);
  }

  if (r_map) {
    *r_map = map;
  } else {
    dl_free(map);
  }

  return GOOSEBERRY_SUCCESS;
}



#endif

