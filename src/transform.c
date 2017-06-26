/**
 * @file transform.c
 * @brief Matrix transformation functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-06-10
 */




#ifndef GOOSEBERRY_TRANSFORM_C
#define GOOSEBERRY_TRANSFORM_C




#include "transform.h"




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int transform_symmetrify_sparse(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  dim_t i, k, snrows;
  ind_t j, l, nnz, m;
  ind_t * srowptr, * shash;
  dim_t * srowind;
  real_t * srowval;

  snrows = dl_max(nrows,ncols);
  shash = ind_init_alloc(NULL_IND,snrows);

  srowptr = ind_init_alloc(0,snrows+1);
  srowind = dim_alloc(rowptr[nrows]*2);
  srowval = real_alloc(rowptr[nrows]*2);

  /* count nnz per row */
  for (i=0;i<nrows;++i) {
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = rowind[j];
      ++srowptr[i+1];
      if (k != i) {
        ++srowptr[k+1];
      }
    }
  }

  ind_prefixsum_exc(srowptr+1,snrows);

  /* fill rowind */
  for (i=0;i<nrows;++i) {
    l = srowptr[i+1];
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = rowind[j];
      srowind[l] = k;
      srowval[l] = rowval[j];
      ++l;
      if (k != i) {
        srowind[srowptr[k+1]] = i;
        srowval[srowptr[k+1]] = rowval[j];
        ++srowptr[k+1];
      }
    }
    srowptr[i+1] = l;
  }

  /* remove/combine duplicate entries and fill rowval */
  nnz = 0;
  j = 0;
  for (i=0;i<snrows;++i) {
    for (;j<srowptr[i+1];++j) {
      k = srowind[j];
      m = shash[k];
      /* someday optimize for null values passed in r_* */
      if (m != NULL_IND) {
        /* combine entries */
        srowval[m] += srowval[j];
      } else {
        /* fill entry */
        shash[k] = nnz;
        srowind[nnz] = k;
        srowval[nnz] = srowval[j];
        ++nnz;
      }
    }
    srowptr[i+1] = nnz;
    /* clear shash */
    for (l=srowptr[i];l<srowptr[i+1];++l) {
      k = srowind[l];
      shash[k] = NULL_IND;
    }
  }

  if (nnz < 1.5 * rowptr[nrows]) {
    /* if there are real gains to be made by shrinking */
    srowind = dim_realloc(srowind,nnz);
    srowval = real_realloc(srowval,nnz);
  }

  dl_free(shash);

  if (r_rowptr) {
    *r_rowptr = srowptr;
  } else {
    dl_free(srowptr);
  }
  if (r_rowind) {
    *r_rowind = srowind;
  } else {
    dl_free(srowind);
  }
  if (r_rowval) {
    *r_rowval = srowval;
  } else {
    dl_free(srowval);
  }

  return GOOSEBERRY_SUCCESS;
}


int transform_debipartify_sparse(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  dim_t i, k, snrows;
  ind_t j, nnz;
  ind_t * srowptr, * colptr;
  dim_t * srowind, * colind;
  real_t * srowval, * colval;

  snrows = nrows + ncols;

  srowptr = ind_alloc(snrows+1);
  srowind = dim_calloc(rowptr[nrows]*2);
  srowval = real_calloc(rowptr[nrows]*2);

  colptr = srowptr + nrows+1;
  colind = srowind + rowptr[nrows];
  colval = srowval + rowptr[nrows];

  ind_set(colptr,0,ncols);

  srowptr[0] = 0;
  for (i=0;i<nrows;++i) {
    srowptr[i+1] = rowptr[i+1];
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = rowind[j];
      /* insert A */
      srowind[j] = k+nrows;
      srowval[j] = rowval[j];
      /* count A^T */
      ++colptr[k];
    }
  }

  ind_prefixsum_exc(colptr,ncols);

  /* complete A^T */
  for (i=0;i<nrows;++i) {
    srowptr[i] = rowptr[i];
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = rowind[j];
      nnz = colptr[k]++;
      colind[nnz] = i;
      colval[nnz] = rowval[j];
    }
  }

  for (i=0;i<ncols;++i) {
    colptr[i] += srowptr[nrows];
  }

  if (r_rowptr) {
    *r_rowptr = srowptr;
  } else {
    dl_free(srowptr);
  }
  if (r_rowind) {
    *r_rowind = srowind;
  } else {
    dl_free(srowind);
  }
  if (r_rowval) {
    *r_rowval = srowval;
  } else {
    dl_free(srowval);
  }

  return GOOSEBERRY_SUCCESS;
}


int transform_rowsplit_sparse(
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
  dim_t i, p, o;
  ind_t j, l;
  ind_t * srowptr;
  dim_t * srowind = NULL, *ldist, * rename;
  real_t * srowval = NULL;

  const ind_t nnz = rowptr[nrows];
  const int doz = r_rowind || r_rowval;

  /* make a modifiable copy of dist */
  ldist = dim_duplicate(dist,nparts);
  rename = dim_alloc(nrows);

  /* allocate row pointers */
  srowptr = ind_alloc(nrows+1);

  /* fill srowptr */
  srowptr[0] = 0;
  for (i=0;i<nrows;++i) {
    p = map[i];
    o = ldist[p]++;
    if (doz) {
      rename[i] = o;
    }
    srowptr[o+1] = rowptr[i+1] - rowptr[i];
  }
  ind_prefixsum_inc(srowptr+1,nrows);

  dl_free(ldist);

  if (doz) {
    /* allocate non-zero elements */
    if (r_rowind) {
      srowind = dim_alloc(nnz);
    }
    if (r_rowval) {
      srowval = real_alloc(nnz);
    }

    for (i=0;i<nrows;++i) {
      o = rename[i];
      l = srowptr[o];
      for (j=rowptr[i];j<rowptr[i+1];++j) {
        if (r_rowind) {
          srowind[l] = rowind[j];
        }
        if (r_rowval) {
          srowval[l] = rowval[j];
        }
        ++l;
      }
    }

    dl_free(rename);

    /* only get allocated if supplied */
    if (r_rowind) {
      *r_rowind = srowind;
    }
    if (r_rowval) {
      *r_rowval = srowval;
    }
  }

  if (r_rowptr) {
    *r_rowptr = srowptr;
  } else {
    dl_free(srowptr);
  }

  return GOOSEBERRY_SUCCESS;
}


int transform_colsplit_sparse(
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
  dim_t i, k, p;
  ind_t j, l;
  ind_t * srowptr;
  dim_t * srowind, * perm, * nperm;
  real_t * srowval;

  nperm = dim_calloc(nparts);
  perm = dim_alloc(ncols);

  srowptr = ind_calloc((nrows*nparts)+1);
  srowind = dim_alloc(rowptr[nrows]);
  srowval = real_alloc(rowptr[nrows]);

  /* create local index numbering */
  for (i=0;i<ncols;++i) {
    perm[i] = nperm[map[i]]++;
  }

  /* determine the length of each sub-row */
  srowptr[0] = 0;
  for (i=0;i<nrows;++i) {
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = rowind[j]; 
      p = map[k];
      ++srowptr[(p*nrows)+i+1];
    }
  }
  ind_prefixsum_exc(srowptr+1,nrows*nparts);

  /* fill in matrices */
  for (i=0;i<nrows;++i) {
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = rowind[j]; 
      p = map[k];
      l = srowptr[(p*nrows)+i+1]++;
      /* fix index */
      srowind[l] = perm[rowind[j]];
      srowval[l] = rowval[j];
    }
  }

  dl_free(perm);
  dl_free(nperm); 

  /* only get allocated if supplied */
  if (r_rowind) {
    *r_rowind = srowind;
  }
  if (r_rowval) {
    *r_rowval = srowval;
  }
  if (r_rowptr) {
    *r_rowptr = srowptr;
  } else {
    dl_free(srowptr);
  }

  return GOOSEBERRY_SUCCESS;
}


int transform_transpose_sparse(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval,
    ind_t ** const r_trowptr,
    dim_t ** const r_trowind,
    real_t ** const r_trowval)
{
  dim_t i;
  ind_t j, k;

  ind_t * trowptr;
  dim_t * trowind;
  real_t * trowval = NULL;

  const ind_t nnz = rowptr[nrows];

  trowptr = ind_calloc(ncols+1);
  trowind = dim_alloc(nnz); 
  if (r_trowval) {
    trowval = real_alloc(nnz); 
  }
  for (i=0;i<nrows;++i) {
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      ++trowptr[rowind[j]+1];
    }
  }
  ind_prefixsum_exc(trowptr+1,ncols);
  for (i=0;i<nrows;++i) {
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = trowptr[rowind[j]+1]++;
      trowind[k] = i;
      if (r_trowval) {
        trowval[k] = rowval[j];
      }
    }
  }

  DL_ASSERT_EQUALS(trowptr[ncols],rowptr[nrows],PF_IND_T);

  if (r_trowptr) {
    *r_trowptr = trowptr;
  } else {
    dl_free(trowptr);
  }
  if (r_trowind) {
    *r_trowind = trowind;
  } else {
    dl_free(trowind);
  }
  if (r_trowval) {
    *r_trowval = trowval;
  }

  return GOOSEBERRY_SUCCESS;
}


#endif
