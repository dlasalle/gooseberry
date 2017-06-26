/**
 * @file permute.c
 * @brief Functions for re-ordering a matrix
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-06-05
 */




#ifndef GOOSEBERRY_PERMUTE_C
#define GOOSEBERRY_PERMUTE_C




#include "permute.h"



/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLPQ_PREFIX dd
#define DLPQ_KEY_T dim_t
#define DLPQ_VAL_T dim_t
#define DLPQ_MIN 1
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_MIN
#undef DLPQ_KEY_T
#undef DLPQ_VAL_T
#undef DLPQ_PREFIX





/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __cuthillmckee_block(
    const dim_t start,
    const dim_t end,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    dim_t * const perm)
{
  dim_t i,k,d,nordered,sr;
  ind_t j;
  dim_t * deg;
  dd_pq_t * q, * rem;

  q = dd_pq_create(start,end);
  rem = dd_pq_create(start,end);
  /* offset pointer */
  deg = dim_alloc(end-start)-start;

  /* find my lowest degree vertex */
  for (i=start;i<end;++i) {
    d = 0;
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = rowind[j]; 
      if (k < end && k >= start) {
        ++d;
      }
    }
    deg[i] = d;
    dd_pq_push(d,i,rem);
  }

  sr = nordered = start;

  /* loop through connected components */
  while (rem->size > 0) {
    i = dd_pq_pop(rem);
    perm[nordered++] = i;

    /* perform bfs */
    while (sr < nordered) {
      i = perm[sr++];
      for (j=rowptr[i];j<rowptr[i+1];++j) {
        k = rowind[j]; 
        if (k < end && k >= start && dd_pq_contains(k,rem)) {
          /* local non-zero */
          dd_pq_remove(k,rem);
          dd_pq_push(deg[k],k,q);
        }
      }
      /* add rows/vertices in ascending order of local degree */
      while (q->size > 0) {
        k = dd_pq_pop(q);
        perm[nordered++] = k;
      }
    }
  }

  dd_pq_free(q);
  dd_pq_free(rem);
  /* un-offset */
  dl_free((deg+start));
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int permute_sparse(
    const dim_t nrows,
    const dim_t ncols,
    ind_t * const rowptr,
    dim_t * const rowind,
    real_t * const rowval,
    const dim_t * const rowperm,
    const dim_t * const colperm)
{
  dim_t i, k, v;
  ind_t j;
  dim_t * rename;
  ind_t * prowptr;
  dim_t * prowind;
  real_t * prowval;

  if (rowperm) {
    /* perform a row and possibly a column permutation */
    prowptr = ind_alloc(nrows+1);
    prowind = dim_alloc(rowptr[nrows]);
    prowval = real_alloc(rowptr[nrows]);
    prowptr[0] = 0;
    for (i=0;i<nrows;++i) {
      v = rowperm[i];
      prowptr[i+1] = prowptr[i] + (rowptr[v+1] - rowptr[v]);
    }
    
    k = 0;
    if (colperm) {
      rename = dim_alloc(ncols);
      /* reverse the colperm */
      for (i=0;i<ncols;++i) {
        v = colperm[i];
        rename[v] = i;
      }
      /* apply the colperm while filling rowindex and rowvalue arrays */
      for (i=0;i<nrows;++i) {
        v = rowperm[i];
        for (j=rowptr[v];j<rowptr[v+1];++j) {
          prowind[k] = rename[rowind[j]];
          prowval[k] = rowval[j];
          ++k;
        }
      }
      dl_free(rename);
    } else {
      /* fill rowindex and rowvalue arrays */
      for (i=0;i<nrows;++i) {
        v = rowperm[i];
        for (j=rowptr[v];j<rowptr[v+1];++j) {
          prowind[k] = rowind[j];
          prowval[k] = rowval[j];
          ++k;
        }
      }
    }
    ind_copy(rowptr,prowptr,nrows);
    dim_copy(rowind,prowind,rowptr[nrows]);
    real_copy(rowval,prowval,rowptr[nrows]);

    dl_free(prowind);
    dl_free(prowptr);
    dl_free(prowval);
  } else if (colperm) {
    /* perform only a column permutation */
    prowind = dim_alloc(rowptr[nrows]);
    rename = dim_alloc(ncols);
    /* reverse the colperm */
    for (i=0;i<ncols;++i) {
      v = colperm[i];
      rename[v] = i;
    }
    /* apply the colperm while filling rowindex and rowvalue arrays */
    for (i=0;i<nrows;++i) {
      for (j=rowptr[i];j<rowptr[i+1];++j) {
        prowind[j] = rename[rowind[j]];
      }
    }

    dim_copy(rowind,prowind,rowptr[nrows]);
    dl_free(rename);
    dl_free(prowind);
  }

  return GOOSEBERRY_SUCCESS;
}


int permute_dense(
    const dim_t nrows,
    const dim_t ncols,
    real_t * const rowval,
    const dim_t * const rowperm,
    const dim_t * const colperm)
{
  dim_t i,v,j,u;
  real_t * prowval;

  if (rowperm || colperm) {
    prowval = real_alloc(nrows*ncols);

    for (i=0;i<nrows;++i) {
      if (rowperm) {
        v = rowperm[i];
      } else {
        v = i;
      }
      if (colperm) {
        for (j=0;j<ncols;++j) {
          u = colperm[j];
          prowval[(i*ncols)+j] = rowval[(v*ncols)+u];
        }
      } else {
        for (j=0;j<ncols;++j) {
          prowval[(i*ncols)+j] = rowval[(v*ncols)+j];
        }
      }
    }

    real_copy(rowval,prowval,nrows*ncols);

    dl_free(prowval);
  }
  
  return GOOSEBERRY_SUCCESS;
}


int permute_dense_rev(
    const dim_t nrows,
    const dim_t ncols,
    real_t * const rowval,
    const dim_t * const rowperm,
    const dim_t * const colperm)
{
  dim_t i,v,j,u;
  real_t * prowval;

  if (rowperm || colperm) {
    prowval = real_alloc(nrows*ncols);

    for (i=0;i<nrows;++i) {
      if (rowperm) {
        v = rowperm[i];
      } else {
        v = i;
      }
      if (colperm) {
        for (j=0;j<ncols;++j) {
          u = colperm[j];
          prowval[(v*ncols)+u] = rowval[(i*ncols)+j];
        }
      } else {
        for (j=0;j<ncols;++j) {
          prowval[(v*ncols)+j] = rowval[(i*ncols)+j];
        }
      }
    }

    real_copy(rowval,prowval,nrows*ncols);

    dl_free(prowval);
  }
  
  return GOOSEBERRY_SUCCESS;
}


int permute_cuthillmckee(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    const dim_t * const blocks,
    const dim_t nblocks,
    dim_t * const perm)
{
  int err;

  if (nrows != ncols) {
    eprintf("Cuthill-McKee requires a structurally symmetric matrix: Given "
        PF_DIM_T"x"PF_DIM_T"\n",nrows,ncols);
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  }

  if (blocks) {
    /* if block is supplied, bandwidth reduction happens per block */
    #pragma omp parallel default(none)
    {
      dim_t c, mystart, myend;

      const dim_t myid = omp_get_thread_num();
      const dim_t nthreads = omp_get_num_threads();

      for (c=myid;c<nblocks;c+=nthreads) {
        mystart = blocks[c];
        myend = blocks[c+1];

        __cuthillmckee_block(mystart,myend,rowptr,rowind,rowval,perm);
      }
    }
  } else {
    __cuthillmckee_block(0,nrows,rowptr,rowind,rowval,perm);
  }

  err = GOOSEBERRY_SUCCESS;

  END:

  return err;
}




#endif
