/**
 * @file analyze.c
 * @brief Functions for determining matrix statistics
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-06-28
 */




#ifndef GOOSEBERRY_ANALYZE_C
#define GOOSEBERRY_ANALYZE_C




#include "analyze.h"
#include "cholesky.h"




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static dim_t __min_active_set_block(
    const dim_t start,
    const dim_t end,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval)
{
  dim_t i, k, nactive, maxactive;
  ind_t j;
  dim_t * first;
  dim_t * last;

  first = dim_init_alloc(NULL_DIM,ncols);
  last = dim_init_alloc(NULL_DIM,ncols);
  
  /* do a first pass to find the first and last use of each element */
  for (i=start;i<end;++i) {
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = rowind[j];
      if (first[k] == NULL_DIM) {
        first[k] = i;
      }
      last[k] = i;
    }
  }

  /* track the size of the active set */
  nactive = 0;
  maxactive = 0;
  for (i=start;i<end;++i) {
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = rowind[j];
      if (first[k] == i) {
        ++nactive;
      }
      if (last[k] == i) {
        --nactive;
      }
      if (nactive > maxactive) {
        maxactive = nactive;
      }
    }
  }

  dl_free(first);
  dl_free(last);

  return maxactive;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int analyze_min_active_set(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    const dim_t * const blocks,
    const dim_t nblocks,
    dim_t * const mas)
{
  #pragma omp parallel default(none)
  {
    dim_t c,mystart,myend;

    const size_t myid = omp_get_thread_num();
    const size_t nthreads = omp_get_num_threads();

    if (blocks) {
      for (c=myid;c<nblocks;c+=nthreads) {
        mystart = blocks[c];
        myend = blocks[c+1];

        mas[c] = __min_active_set_block(mystart,myend,ncols,rowptr,rowind,
            rowval);
      }
    } else {
      mas[0] = __min_active_set_block(0,nrows,ncols,rowptr,rowind,rowval);
    }
  }

  return GOOSEBERRY_SUCCESS;
}


void analyze_cholesky(
    dim_t const nrows,
    ind_t const * const rowptr,
    dim_t const * const rowind,
    ind_t * const r_nnz,
    double * const r_nops)
{
  dim_t i;
  double nops;
  ind_t nnz;
  dim_t * counts;

  counts = dim_alloc(nrows);

  nops = 0;
  nnz = 0;

  cholesky_rowcounts(nrows,rowptr,rowind,counts);

  for (i=0;i<nrows;++i) {
    nops += counts[i]*counts[i];
    nnz += counts[i];
  }

  dl_free(counts);

  *r_nnz = nnz;
  *r_nops = nops;
}




#endif
