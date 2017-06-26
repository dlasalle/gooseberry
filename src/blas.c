/**
 * @file dlblas_headers.h
 * @brief Basic linear algebra functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-09-29
 */




#ifndef GOOSEBERRY_BLAS_C
#define GOOSEBERRY_BLAS_C




#include "blas.h"

#ifndef NO_OMP
#include "omp.h"
#endif




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define IDX(i,j,m) ((((ind_t)i)*((ind_t)m))+((ind_t)j))




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __dense_mult_dense_block(
    const dim_t start,
    const dim_t end,
    const dim_t ancols,
    const dim_t bncols,
    const real_t * const arowval,
    const real_t * const bcolval,
    real_t * const crowval)
{
  dim_t i,j,k;
  real_t x;

  for (i=start;i<end;++i) {
    for (j=0;j<bncols;++j) {
      x = 0;
      for (k=0;k<ancols;++k) {
        x += arowval[IDX(i,k,ancols)]*bcolval[IDX(j,k,ancols)];
      }
      crowval[IDX(i,j,bncols)] = x; 
    }
  }
}


static void __sparse_mult_dense_block(
    const dim_t start,
    const dim_t end,
    const dim_t ancols,
    const dim_t bncols,
    const ind_t * const arowptr,
    const dim_t * const arowind,
    const real_t * const arowval,
    const real_t * const bcolval,
    real_t * const crowval)
{
  dim_t i,j,l;
  ind_t k;
  real_t x;

  for (i=start;i<end;++i) {
    for (j=0;j<bncols;++j) {
      x = 0;
      for (k=arowptr[i];k<arowptr[i+1];++k) {
        l = arowind[k]; 
        x += arowval[k]*bcolval[IDX(j,l,ancols)];
      }
      crowval[IDX(i,j,bncols)] = x; 
    }
  }
}


static void __sparse_mult_sparse_block(
    const dim_t start,
    const dim_t end,
    const dim_t bncols,
    const ind_t * const arowptr,
    const dim_t * const arowind,
    const real_t * const arowval,
    const ind_t * const browptr,
    const dim_t * const browind,
    const real_t * const browval,
    real_t * const crowval)
{
  dim_t i,l;
  ind_t j,k;
  real_t x;

  for (i=start;i<end;++i) {
    real_set(crowval+(i*bncols),0,bncols);
    for (j=arowptr[i];j<arowptr[i+1];++j) {
      l = arowind[j];
      x = arowval[j];
      for (k=browptr[l];k<browptr[l+1];++k) {
        crowval[IDX(i,browind[k],bncols)] += x*browval[k];
      }
    }
  }
}


static void __sparse_mult_sparse_sparse_count(
    const dim_t start,
    const dim_t end,
    const dim_t bncols,
    const ind_t * const arowptr,
    const dim_t * const arowind,
    const ind_t * const browptr,
    const dim_t * const browind,
    ind_t * const crowptr)
{
  dim_t i,l, rowsize, m;
  ind_t j,k;

  int * row;
  dim_t * mark;

  if (start < end) {
    row = int_calloc(bncols);
    mark = dim_alloc(bncols);

    /* do the first iteration and avoid the offet */
    rowsize = 0;
    for (j=arowptr[start];j<arowptr[start+1];++j) {
      l = arowind[j];
      for (k=browptr[l];k<browptr[l+1];++k) {
        m = browind[k];
        if (!row[m]) {
          row[m] = 1;
          mark[rowsize++] = m;
        }
      }
    }
    crowptr[start+1] = rowsize;
    /* clear row */
    for (m=0;m<rowsize;++m) {
      row[mark[m]] = 0;
    }
    for (i=start+1;i<end;++i) {
      rowsize = 0;
      for (j=arowptr[i];j<arowptr[i+1];++j) {
        l = arowind[j];
        for (k=browptr[l];k<browptr[l+1];++k) {
          m = browind[k];
          if (!row[m]) {
            row[m] = 1;
            mark[rowsize++] = m;
          }
        }
      }
      crowptr[i+1] = rowsize + crowptr[i];
      /* clear row */
      for (m=0;m<rowsize;++m) {
        row[mark[m]] = 0;
      }
    }

    dl_free(row);
    dl_free(mark);
  }
}


static void __sparse_mult_sparse_sparse_fill(
    const dim_t start,
    const dim_t end,
    const dim_t bncols,
    const ind_t * const arowptr,
    const dim_t * const arowind,
    const real_t * const arowval,
    const ind_t * const browptr,
    const dim_t * const browind,
    const real_t * const browval,
    ind_t * const crowptr,
    dim_t * const crowind,
    real_t * const crowval)
{
  dim_t i,l,m;
  ind_t j,k,n,ptr;
  real_t x;

  ind_t * row;

  row = ind_init_alloc(NULL_IND,bncols);

  /* fill matrix */
  for (i=start;i<end;++i) {
    ptr = crowptr[i];
    for (j=arowptr[i];j<arowptr[i+1];++j) {
      l = arowind[j];
      x = arowval[j];
      for (k=browptr[l];k<browptr[l+1];++k) {
        m = browind[k];
        if ((n = row[m]) == NULL_IND) {
          n = ptr++;
          crowind[n] = m;
          crowval[n] = x * browval[k];
          row[m] = n;
        } else {
          crowval[n] += x * browval[k];
        }
      }
    }
    /* clear row and insert intro matrix */
    for (j=crowptr[i];j<crowptr[i+1];++j) {
      l = crowind[j];
      row[l] = NULL_IND;
    }
  }

  dl_free(row);
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int blas_scale(
    const ind_t n,
    real_t * const rowval,
    const real_t s)
{
  ind_t i;

  #pragma omp parallel for default(none) schedule(static,OMP_BIG_BLOCK)
  for (i=0;i<n;++i) {
    rowval[i] *= s;
  }

  return GOOSEBERRY_SUCCESS;
}


int blas_negate(
    const ind_t n,
    real_t * const rowval)
{
  ind_t i;

  #pragma omp parallel for default(none) schedule(static,OMP_BIG_BLOCK)
  for (i=0;i<n;++i) {
    rowval[i] = -rowval[i];
  }

  return GOOSEBERRY_SUCCESS;
}


int blas_add(
    const dim_t anrows,
    const dim_t ancols,
    const real_t * const arowval,
    const real_t * const browval,
    real_t * const crowval)
{
  ind_t i;

  #pragma omp parallel for default(none) schedule(static,OMP_BIG_BLOCK)
  for (i=0;i<(ind_t)anrows*(ind_t)ancols;++i) {
    crowval[i] = arowval[i]+browval[i];
  }

  return GOOSEBERRY_SUCCESS;
}


int blas_addsp(
    const dim_t anrows,
    const dim_t ancols,
    const real_t * const arowval,
    const ind_t * const browptr,
    const dim_t * const browind,
    const dim_t * const browval,
    real_t * const crowval)
{
  dim_t i;
  ind_t j;

  if (crowval != arowval) {
    real_copy(crowval,arowval,(ind_t)anrows*(ind_t)ancols);
  }

  #pragma omp parallel for default(none) private(j) \
    schedule(static,OMP_SMALL_BLOCK)
  for (i=0;i<anrows;++i) {
    for (j=browptr[i];j<browptr[i+1];++j) {
      crowval[(i*ancols)+browind[j]] += browval[j];
    }
  }

  return GOOSEBERRY_SUCCESS;
}


int blas_sub(
    const dim_t anrows,
    const dim_t ancols,
    const real_t * const arowval,
    const real_t * const browval,
    real_t * const crowval)
{
  ind_t i;

  #pragma omp parallel for default(none) schedule(static,OMP_BIG_BLOCK)
  for (i=0;i<(ind_t)anrows*(ind_t)ancols;++i) {
    crowval[i] = arowval[i]-browval[i];
  }

  return GOOSEBERRY_SUCCESS;
}


int blas_subsp(
    const dim_t anrows,
    const dim_t ancols,
    const real_t * const arowval,
    const ind_t * const browptr,
    const dim_t * const browind,
    const real_t * const browval,
    real_t * const crowval)
{
  dim_t i;
  ind_t j;

  if (crowval != arowval) {
    real_copy(crowval,arowval,anrows*ancols);
  }

  #pragma omp parallel for default(none) private(j) \
    schedule(static,OMP_SMALL_BLOCK)
  for (i=0;i<anrows;++i) {
    for (j=browptr[i];j<browptr[i+1];++j) {
      crowval[(i*ancols)+browind[j]] -= browval[j];
    }
  }

  return GOOSEBERRY_SUCCESS;
}


int blas_spsub(
    const dim_t anrows,
    const dim_t ancols,
    const ind_t * const arowptr,
    const dim_t * const arowind,
    const real_t * const arowval,
    const real_t * const browval,
    real_t * const crowval)

{
  dim_t i;
  ind_t j;

  #pragma omp parallel for default(none) schedule(static,OMP_BIG_BLOCK)
  for (i=0;i<(ind_t)anrows*(ind_t)ancols;++i) {
    crowval[i] = -browval[i];
  }

  #pragma omp parallel for default(none) private(j) \
    schedule(static,OMP_SMALL_BLOCK)
  for (i=0;i<anrows;++i) {
    for (j=arowptr[i];j<arowptr[i+1];++j) {
      crowval[(i*ancols)+arowind[j]] += arowval[j];
    }
  }

  return GOOSEBERRY_SUCCESS;
}


int blas_add_scalar(
    const ind_t n,
    real_t * const rowval,
    const real_t s)
{
  ind_t i;

  #pragma omp parallel for default(none) schedule(static,OMP_BIG_BLOCK)
  for (i=0;i<n;++i) {
    rowval[i] += s;
  }

  return GOOSEBERRY_SUCCESS;
}


real_t blas_dot(
    const dim_t n,
    const real_t * const arowval,
    const real_t * const bcolval)
{
  dim_t i;
  real_t x;

  x = 0;

  #pragma omp parallel for default(none) schedule(static,OMP_BIG_BLOCK) \
    reduction(+:x)
  for (i=0;i<n;++i) {
    x += arowval[i]*bcolval[i];
  }

  return x;
}


real_t blas_spdot(
    const dim_t annz,
    const real_t * const arowind, 
    const real_t * const arowval,
    const real_t * const bcolval)
{
  dim_t i,j;
  real_t x;

  x = 0;
  #pragma omp parallel for default(none) schedule(static,OMP_BIG_BLOCK) \
    private(j) reduction(+:x)
  for (i=0;i<annz;++i) {
    j = arowind[i];
    x += arowval[i]*bcolval[j];
  }

  return x;
}


real_t blas_spdotsp(
    const dim_t annz,
    const dim_t bnnz,
    const real_t * const arowind, 
    const real_t * const arowval,
    const real_t * const bcolind, 
    const real_t * const bcolval)
{
  dim_t i,j,k,l;
  real_t x = 0;

  if (dl_min(annz,bnnz) > OMP_BIG_BLOCK) {
    #pragma omp parallel default(none) reduction(+:x) private(i,j,k,l)
    {
      const size_t myid = omp_get_thread_num();
      const size_t nthreads = omp_get_num_threads();

      const dim_t mystart = size_chunkstart(myid,nthreads,annz);
      const dim_t myend = mystart + size_chunksize(myid,nthreads,annz);

      k = 0;
      for (i=mystart;i<myend;++i) {
        j = arowind[i];
        while (bcolind[k] < j) {
          ++k;
          if (k >= bnnz) {
            goto END_PAR;
          }
        }
        l = bcolind[k];
        if (j == l) {
          x += arowval[j]*bcolval[l];
        }
      }

      END_PAR:;
    }
  } else {
    k = 0;
    for (i=0;i<annz;++i) {
      j = arowind[i];
      while (bcolind[k] < j) {
        ++k;
        if (k >= bnnz) {
          goto END_SER;
        }
      }
      l = bcolind[k];
      if (j == l) {
        x += arowval[j]*bcolval[l];
      }
    }

    END_SER:;
  }

  return x;
}


int blas_mult(
    const dim_t anrows,
    const dim_t ancols,
    const dim_t bncols, 
    const real_t * const arowval,
    const real_t * const bcolval,
    real_t * const crowval,
    const dim_t * const blocks,
    const dim_t nblocks)
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

        __dense_mult_dense_block(mystart,myend,ancols,bncols,arowval,bcolval,
            crowval);
      }
    } else {
      mystart = size_chunkstart(myid,nthreads,anrows);
      myend = mystart + size_chunksize(myid,nthreads,anrows);

      __dense_mult_dense_block(mystart,myend,ancols,bncols,arowval,bcolval,
          crowval);
    }
  }

  return GOOSEBERRY_SUCCESS;
}


int blas_spmult(
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
  /* someday use block decomposition if b is a thick enough matrix */
  #pragma omp parallel default(none)
  {
    dim_t c,mystart,myend;

    #ifdef CONTINUOUS_CHUNKS
    dim_t bstart,bend;
    #endif

    const size_t myid = omp_get_thread_num();
    const size_t nthreads = omp_get_num_threads();

    if (blocks) {
      #ifdef CONTINUOUS_CHUNKS
      bstart = dim_chunkstart(myid,nblocks,nthreads);
      bend = bstart + dim_chunksize(myid,nblocks,nthreads);
      for (c=bstart;c<bend;++c) {
        mystart = blocks[c];
        myend = blocks[c+1];

        __sparse_mult_dense_block(mystart,myend,ancols,bncols,arowptr,arowind,
            arowval,bcolval,crowval);
      }
      #else
      for (c=myid;c<nblocks;c+=nthreads) {
        mystart = blocks[c];
        myend = blocks[c+1];

        __sparse_mult_dense_block(mystart,myend,ancols,bncols,arowptr,arowind,
            arowval,bcolval,crowval);
      }
      #endif
    } else {
      mystart = size_chunkstart(myid,nthreads,anrows);
      myend = mystart + size_chunksize(myid,nthreads,anrows);

      __sparse_mult_dense_block(mystart,myend,ancols,bncols,arowptr,arowind,
          arowval,bcolval,crowval);
    }
  }

  return GOOSEBERRY_SUCCESS;
}


int blas_spmultsp(
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
  #pragma omp parallel default(none)
  {
    dim_t c,mystart,myend;

    const size_t myid = omp_get_thread_num();
    const size_t nthreads = omp_get_num_threads();

    if (blocks) {
      for (c=myid;c<nblocks;c+=nthreads) {
        mystart = blocks[c];
        myend = blocks[c+1];

        __sparse_mult_sparse_block(mystart,myend,bncols,arowptr,arowind,
            arowval,browptr,browind,browval,crowval);
      }
    } else {
      mystart = size_chunkstart(myid,nthreads,anrows);
      myend = mystart + size_chunksize(myid,nthreads,anrows);

      __sparse_mult_sparse_block(mystart,myend,bncols,arowptr,arowind,arowval,
          browptr,browind,browval,crowval);
    }
  }

  return GOOSEBERRY_SUCCESS;
}


int blas_spmultsp_sp(
    const dim_t anrows,
    const dim_t ancols,
    const dim_t bncols,
    const ind_t * const arowptr,
    const dim_t * const arowind,
    const real_t * const arowval,
    const ind_t * const browptr,
    const dim_t * const browind,
    const real_t * const browval,
    ind_t ** const r_crowptr,
    dim_t ** const r_crowind,
    real_t ** const r_crowval,
    const dim_t * const blocks,
    const dim_t nblocks)
{
  dim_t * crowind;
  ind_t * crowptr;
  real_t * crowval;
  ind_t * prefix;

  /* the new matrix */
  crowptr = ind_init_alloc(0,anrows+1);

  #pragma omp parallel default(none) shared(crowptr,crowind,crowval,prefix)
  {
    dim_t c,mystart,myend,i;
    ind_t offset;

    const size_t nthreads = omp_get_num_threads();
    const size_t myid = omp_get_thread_num();

    if (blocks) {
      /* verifiy this still */
      #pragma omp barrier
      #pragma omp master
      {
        prefix = ind_alloc(nblocks);
      }
      #pragma omp barrier

      /* generate initial offset */
      for (c=myid;c<nblocks;c+=nthreads) {
        mystart = blocks[c];
        myend = blocks[c+1];

        __sparse_mult_sparse_sparse_count(mystart,myend,bncols,arowptr,arowind,
            browptr,browind,crowptr);
        prefix[c] = crowptr[myend];
      }
      #pragma omp barrier

      /* finish prefix sum on crowptr */
      #pragma omp master
      {
        ind_prefixsum_exc(prefix,nblocks); 
      }
      #pragma omp barrier
      for (c=myid;c<nblocks;c+=nthreads) {
        mystart = blocks[c];
        myend = blocks[c+1];

        offset = prefix[c];

        for (i=mystart;i<myend;++i) {
          crowptr[i+1] += offset;
        }
      }
      #pragma omp barrier

      /* allocate indices and values */
      #pragma omp master
      {
        /* allocate the rest of the output matrix */
        crowind = dim_alloc(crowptr[anrows]);
        crowval = real_alloc(crowptr[anrows]);
      }
      #pragma omp barrier

      /* fill in indicies and values */
      for (c=myid;c<nblocks;c+=nthreads) {
        mystart = blocks[c];
        myend = blocks[c+1];

        __sparse_mult_sparse_sparse_fill(0,anrows,bncols,arowptr,arowind,
            arowval,browptr,browind,browval,crowptr,crowind,crowval);
      }
    } else {
      #pragma omp barrier
      #pragma omp master
      {
        prefix = ind_alloc(nthreads);
      }
      #pragma omp barrier

      mystart = size_chunkstart(myid,nthreads,anrows);
      myend = mystart + size_chunksize(myid,nthreads,anrows);

      __sparse_mult_sparse_sparse_count(mystart,myend,bncols,arowptr,arowind,
          browptr,browind,crowptr);
      if (myend > mystart) {
        prefix[myid] = crowptr[myend];
      } else {
        prefix[myid] = 0;
      }
      #pragma omp barrier

      /* finish prefix sum on crowptr */
      #pragma omp master
      {
        ind_prefixsum_exc(prefix,nthreads); 
      }
      #pragma omp barrier
      offset = prefix[myid];

      for (i=mystart;i<myend;++i) {
        crowptr[i+1] += offset;
      }

      #pragma omp barrier
      #pragma omp master
      {
        /* allocate the rest of the output matrix */
        crowind = dim_alloc(crowptr[anrows]);
        crowval = real_alloc(crowptr[anrows]);
      }
      #pragma omp barrier

      __sparse_mult_sparse_sparse_fill(mystart,myend,bncols,arowptr,arowind,
          arowval,browptr,browind,browval,crowptr,crowind,crowval);
    }
  }

  dl_free(prefix);

  *r_crowptr = crowptr;
  *r_crowind = crowind;
  *r_crowval = crowval;

  return GOOSEBERRY_SUCCESS;
}





#endif
