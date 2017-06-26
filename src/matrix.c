/**
 * @file matrix.c
 * @brief Functions for creating and manipulating matrices 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-04-30
 */




#ifndef GOOSEBERRY_MATRIX_C
#define GOOSEBERRY_MATRIX_C




#include "matrix.h"
#include "permute.h"




/******************************************************************************
* MEMORY FUNCTIONS ************************************************************
******************************************************************************/


#define DLMEM_PREFIX matrix
#define DLMEM_TYPE_T matrix_t
#include "dlmem_funcs.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int matrix_init(
    const int type, 
    const dim_t nrows, 
    const dim_t ncols,
    const ind_t nnz,
    matrix_t * const mat)
{
  mat->type = type;
  mat->nrows = nrows;
  mat->ncols = ncols;
  mat->ncdist = 0;
  mat->nrdist = 0;
  mat->rdist = NULL;
  mat->cdist = NULL;

  if (nrows == 0 || ncols == 0) {
    eprintf("Invalid dimensions for all types: "PF_DIM_T"x"PF_DIM_T"\n",
        nrows,ncols);
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }

  switch (type) {
    case MATRIX_TYPE_DENSE_VECTOR:
      if (nrows > 1 && ncols != 1) {
          eprintf("Invalid dense vector dimensions: "PF_DIM_T"x"PF_DIM_T"\n",
              nrows, ncols);
        return GOOSEBERRY_ERROR_INVALIDINPUT;
      }
      /* nrows or ncols has to equal 1 */
      mat->rowval = mat->colval = real_alloc(nrows*ncols);
      mat->rowptr = mat->colptr = NULL;
      mat->rowind = mat->colind = NULL;
      break;
    case MATRIX_TYPE_SPARSE_VECTOR:
      if (nrows > 1) {
        if (ncols != 1) {
          eprintf("Invalid sparse vector dimensions: "PF_DIM_T"x"PF_DIM_T"\n",
              nrows, ncols);
          return GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        mat->rowptr = ind_alloc(2);
        mat->rowptr[0] = 0;
        mat->rowptr[1] = nnz;
        if (nnz != NULL_IND) {
          mat->rowind = dim_alloc(nnz);
          mat->rowval = real_alloc(nnz);
        } else {
          mat->rowind = NULL;
          mat->rowval = NULL;
        }
        mat->colptr = NULL;
        mat->colind = NULL;
        mat->colval = NULL;
      } else {
        if (nrows != 1) {
          eprintf("Invalid sparse vector dimensions: "PF_DIM_T"x"PF_DIM_T"\n",
              nrows, ncols);
          return GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        mat->colptr = ind_alloc(2);
        mat->colptr[0] = 0;
        mat->colptr[1] = nnz;
        if (nnz != NULL_IND) {
          mat->colind = dim_alloc(nnz);
          mat->colval = real_alloc(nnz);
        }  else {
          mat->colind = NULL;
          mat->colval = NULL;
        }
        mat->rowptr = NULL;
        mat->rowind = NULL;
        mat->rowval = NULL;
      } 
      break;
    case MATRIX_TYPE_DENSE_MATRIX:
      mat->rowval = real_alloc(nrows*ncols);
      mat->colval = NULL;
      mat->rowptr = mat->colptr = NULL;
      mat->rowind = mat->colind = NULL;
      break;
    case MATRIX_TYPE_SPARSE_MATRIX:
      if (nnz != NULL_IND) {
        mat->rowptr = ind_alloc(nrows+1);
        mat->rowind = dim_alloc(nnz);
        mat->rowval = real_alloc(nnz);
      } else {
        mat->rowptr = NULL;
        mat->rowind = NULL;
        mat->rowval = NULL;
      }
      mat->colptr = NULL;
      mat->colind = NULL;
      mat->colval = NULL;
      break;
    default:
      eprintf("Unknown matrix type: %d\n",type);
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  return GOOSEBERRY_SUCCESS;
}


int matrix_zero(
    matrix_t * const mat)
{
  switch (mat->type) {
    case MATRIX_TYPE_DENSE_VECTOR:
      real_set(mat->rowval,0,mat->nrows*mat->ncols);
      break;
    case MATRIX_TYPE_DENSE_MATRIX:
      if (mat->rowval) {
        real_set(mat->rowval,0,mat->nrows*mat->ncols);
      } 
      if (mat->colval) {
        real_set(mat->colval,0,mat->nrows*mat->ncols);
      }
      break;
    default:
      eprintf("Cannot zero sparse structures\n");
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  return GOOSEBERRY_SUCCESS;
}


int matrix_identify(
    matrix_t * const mat)
{
  switch (mat->type) {
    case MATRIX_TYPE_DENSE_VECTOR:
    case MATRIX_TYPE_SPARSE_VECTOR:
    case MATRIX_TYPE_DENSE_MATRIX:
    case MATRIX_TYPE_SPARSE_MATRIX:
    default:
      eprintf("Unknown matrix type: %d\n",mat->type);
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  return GOOSEBERRY_SUCCESS;
}


int matrix_sparsify(
    matrix_t * const mat)
{  
  switch (mat->type) {
    case MATRIX_TYPE_SPARSE_VECTOR:
    case MATRIX_TYPE_SPARSE_MATRIX:
      wprintf("Cannot sparsify sparse structures\n");
      break;
    case MATRIX_TYPE_DENSE_VECTOR:
    case MATRIX_TYPE_DENSE_MATRIX:
    default:
      eprintf("Unknown matrix type: %d\n",mat->type);
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  return GOOSEBERRY_SUCCESS;
}


int matrix_densify(
    matrix_t * const mat)
{
  switch (mat->type) {
    case MATRIX_TYPE_SPARSE_VECTOR:
    case MATRIX_TYPE_SPARSE_MATRIX:
      wprintf("Cannot sparsify sparse structures\n");
      break;
    case MATRIX_TYPE_DENSE_VECTOR:
    case MATRIX_TYPE_DENSE_MATRIX:
    default:
      eprintf("Unknown matrix type: %d\n",mat->type);
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  return GOOSEBERRY_SUCCESS;
}


int matrix_transpose(
    matrix_t * const mat)
{
  dl_swap(mat->nrows,mat->ncols);
  dl_swap(mat->rowval,mat->colval);
  /* these are no-ops for dense structures, but easier to do it than check */
  dl_swap(mat->rowptr,mat->colptr);
  dl_swap(mat->rowind,mat->colind);

  return GOOSEBERRY_SUCCESS;
}



int matrix_permute(
    matrix_t * const mat,
    const dim_t * const rperm,
    const dim_t * const cperm)
{
  int rv;

  const dim_t nrows = mat->nrows;
  const dim_t ncols = mat->ncols;

  ind_t * const rowptr = mat->rowptr;
  ind_t * const colptr = mat->colptr;
  dim_t * const rowind = mat->rowind;
  dim_t * const colind = mat->colind;
  real_t * const rowval = mat->rowval;
  real_t * const colval = mat->colval;

  switch (mat->type) {
    case MATRIX_TYPE_DENSE_VECTOR:
      /* re-order row representation only */
      if (rowval) {
        if ((rv = permute_dense(nrows,ncols,rowval,rperm,cperm)) != 
            GOOSEBERRY_SUCCESS) {
          return rv;
        }
      }
      break;
    case MATRIX_TYPE_DENSE_MATRIX:
      /* re-order row representation */
      if (rowval) {
        if ((rv = permute_dense(nrows,ncols,rowval,rperm,cperm)) != 
            GOOSEBERRY_SUCCESS) {
          return rv;
        }
      }
      /* re-order column representation */
      if (colval) {
        if ((rv = permute_dense(ncols,nrows,colval,cperm,rperm)) != 
            GOOSEBERRY_SUCCESS) {
          return rv;
        }
      }
      break;
    case MATRIX_TYPE_SPARSE_VECTOR:
    case MATRIX_TYPE_SPARSE_MATRIX:
      /* re-order row representation */
      if (rowptr) {
        if ((rv = permute_sparse(nrows,ncols,rowptr,rowind,rowval,rperm,
            cperm)) != GOOSEBERRY_SUCCESS) {
          return rv;
        }
      }
      /* re-order column representation */
      if (colptr) {
        if ((rv = permute_sparse(ncols,nrows,colptr,colind,colval,cperm,
            rperm)) != GOOSEBERRY_SUCCESS) {
          return rv;
        }
      }
      break;
    default:
      eprintf("Unknown matrix type: %d\n",mat->type);
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  return GOOSEBERRY_SUCCESS;
}


int matrix_unpermute(
    matrix_t * mat,
    const dim_t * rperm,
    const dim_t * cperm)
{
  int rv = GOOSEBERRY_SUCCESS;
  dim_t i, v;
  dim_t * rren = NULL, * cren = NULL;

  const dim_t nrows = mat->nrows;
  const dim_t ncols = mat->ncols;

  ind_t * const rowptr = mat->rowptr;
  ind_t * const colptr = mat->colptr;
  dim_t * const rowind = mat->rowind;
  dim_t * const colind = mat->colind;
  real_t * const rowval = mat->rowval;
  real_t * const colval = mat->colval;

  switch (mat->type) {
    case MATRIX_TYPE_DENSE_VECTOR:
      /* re-order row representation only */
      if (rowval) {
        if ((rv = permute_dense_rev(nrows,ncols,rowval,rperm,NULL)) 
            != GOOSEBERRY_SUCCESS) {
          goto END;
        }
      }
      break;
    case MATRIX_TYPE_DENSE_MATRIX:
      /* re-order row representation */
      if (rowval) {
        if ((rv = permute_dense_rev(nrows,ncols,rowval,rperm,cperm)) 
            != GOOSEBERRY_SUCCESS) {
          goto END;
        }
      }
      /* re-order column representation */
      if (colval) {
        if ((rv = permute_dense_rev(ncols,nrows,colval,cperm,rperm)) 
            != GOOSEBERRY_SUCCESS) {
          goto END;
        }
      }
      break;
    case MATRIX_TYPE_SPARSE_VECTOR:
    case MATRIX_TYPE_SPARSE_MATRIX:
      if (rperm) {
        rren = dim_alloc(nrows);
        for (i=0;i<nrows;++i) {
          v = rperm[i];
          rren[v] = i;
        }
      }
      if (cperm) {
        cren = dim_alloc(ncols);
        for (i=0;i<ncols;++i) {
          v = cperm[i];
          cren[v] = i;
        }
      }
      /* re-order row representation */
      if (rowptr) {
        if ((rv = permute_sparse(nrows,ncols,rowptr,rowind,rowval,
            rren,cren)) != GOOSEBERRY_SUCCESS) {
          goto END;
        }
      }
      /* re-order column representation */
      if (colptr) {
        if ((rv = permute_sparse(ncols,nrows,colptr,colind,colval,
            cren,rren)) != GOOSEBERRY_SUCCESS) {
          goto END;
        }
      }
      break;
    default:
      eprintf("Unknown matrix type: %d\n",mat->type);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto END;
  }

  END:

  if (rren) {
    dl_free(rren);
  }
  if (cren) {
    dl_free(cren);
  }

  return rv;
}


int matrix_buildindex(
    matrix_t * const mat)
{
  dim_t i;
  ind_t j,k,nnz;

  switch (mat->type) {
    case MATRIX_TYPE_SPARSE_VECTOR:
    case MATRIX_TYPE_SPARSE_MATRIX:
      if (mat->rowptr && mat->colptr) {
        /* do nothing */
      } else if (mat->rowptr) {
        nnz = mat->rowptr[mat->nrows];
        mat->colptr = ind_calloc(mat->ncols+1);
        mat->colind = dim_alloc(nnz); 
        mat->colval = real_alloc(nnz); 
        for (i=0;i<mat->nrows;++i) {
          for (j=mat->rowptr[i];j<mat->rowptr[i+1];++j) {
            ++mat->colptr[mat->rowind[j]+1];
          }
        }
        ind_prefixsum_exc(mat->colptr+1,mat->ncols);
        for (i=0;i<mat->nrows;++i) {
          for (j=mat->rowptr[i];j<mat->rowptr[i+1];++j) {
            k = mat->colptr[mat->rowind[j]+1]++;
            mat->colind[k] = i;
            mat->colval[k] = mat->rowval[j];
          }
        }
      } else {
        nnz = mat->colptr[mat->ncols];
        mat->rowptr = ind_calloc(mat->ncols+1);
        mat->rowind = dim_alloc(nnz); 
        mat->rowval = real_alloc(nnz); 
        for (i=0;i<mat->ncols;++i) {
          for (j=mat->colptr[i];j<mat->colptr[i+1];++j) {
            ++mat->rowptr[mat->colind[j]+1];
          }
        }
        ind_prefixsum_exc(mat->rowptr+1,mat->nrows);
        for (i=0;i<mat->ncols;++i) {
          for (j=mat->colptr[i];j<mat->colptr[i+1];++j) {
            k = mat->rowptr[mat->colind[j]+1]++;
            mat->rowind[k] = (dim_t)i;
            mat->rowval[k] = mat->colval[j];
          }
        }
      }
      break;
    case MATRIX_TYPE_DENSE_VECTOR:
      if (mat->rowval && mat->colval) {
        /* do nothing */
      } else if (mat->rowval) {
        mat->colval = mat->rowval;
      } else {
        mat->rowval = mat->colval;
      }
      break;
    case MATRIX_TYPE_DENSE_MATRIX:
      if (mat->rowval && mat->colval) {
        /* do nothing */
      } else if (mat->rowval) {
        mat->colval = real_alloc(mat->nrows*mat->ncols);
        for (i=0;i<mat->nrows;++i) {
          for (j=0;j<mat->ncols;++j) {
            mat->colval[j*mat->nrows+i] = mat->rowval[i*mat->ncols+j];
          }
        }
      } else {
        mat->rowval = real_alloc(mat->nrows*mat->ncols);
        for (i=0;i<mat->ncols;++i) {
          for (j=0;j<mat->nrows;++j) {
            mat->rowval[j*mat->ncols+i] = mat->colval[i*mat->nrows+j];
          }
        }
      }
      break;
    default:
      eprintf("Unknown matrix type: %d\n",mat->type);
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  return GOOSEBERRY_SUCCESS;
}


int matrix_free(
    matrix_t * mat)
{
  if (mat->rowptr) {
    dl_free(mat->rowptr);
    mat->rowptr = NULL;
  }
  if (mat->colptr) {
    dl_free(mat->colptr);
    mat->colptr = NULL;
  }
  if (mat->rowind) {
    dl_free(mat->rowind);
    mat->rowind = NULL;
  }
  if (mat->colind) {
    dl_free(mat->colind);
    mat->colind = NULL;
  }
  if (mat->rowval == mat->colval) {
    if (mat->rowval) {
      dl_free(mat->rowval);
      mat->rowval = NULL;
    }
  } else {
    if (mat->rowval) {
      dl_free(mat->rowval);
      mat->rowval = NULL;
    } 
    if (mat->colval) {
      dl_free(mat->colval);
      mat->colval = NULL;
    }
  }
  if (mat->rdist != mat->cdist) {
    if (mat->rdist) {
      dl_free(mat->rdist);
      mat->rdist = NULL;
    }
    if (mat->cdist) {
      dl_free(mat->cdist);
      mat->cdist = NULL;
    }
  } else if (mat->rdist) {
    dl_free(mat->rdist);
    mat->rdist = NULL;
    mat->cdist = NULL;
  }

  dl_free(mat);

  return GOOSEBERRY_SUCCESS;
}


#endif
