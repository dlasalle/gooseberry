/**
 * @file cholesky.c
 * @brief Cholesky factorization functions.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2015-02-02
 */





#ifndef GOOSEBERRY_CHOLESKY_C
#define GOOSEBERRY_CHOLESKY_C




#include "cholesky.h"




/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


static inline int __cs_flip(
    int const i)
{
  return - 2 - i;
}


static inline int __cs_unflip(
    int const i)
{
  return (i < 0) ? __cs_flip(i) : i;
}


static inline int __cs_marked(
    int const * const w,
    int const i)
{
  return w[i] < 0;
}


static inline void __cs_mark(
    int * const w,
    int const i)
{
  w[i] = __cs_flip(w[i]);
}




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __firstdesc(
    dim_t const n,
    dim_t const * const parent,
    dim_t const * const post,
    dim_t * const first,
    dim_t * const level)
{
  dim_t len, i, k, r, s;

  dim_set(first,NULL_DIM,n);
  
  for (k=0;k<n;++k) {
    i = post[k];
    len = 0;
    for (r=i;r!=NULL_DIM&&first[r]==NULL_DIM;r=parent[r],++len) {
      first[r] = k;
    }
    if (r == NULL_DIM) {
      --len;
    } else {
      len += level[r];
    }

    for (s=i;s!=r;s=parent[s]) {
      level[s] = len--;
    }
  }
}



/**
 * @brief 
 *
 * @param i
 * @param k
 * @param head
 * @param next
 * @param post
 * @param stack
 *
 * @return 
 */
static dim_t __traverse_dfs(
    dim_t const i,
    dim_t k,
    dim_t * const head,
    dim_t const * const next,
    dim_t * const post,
    dim_t * const stack)
{
  dim_t j, m;
  ssize_t top;

  /* traverse tree in dfs */
  top = 0;
  stack[0] = i; 
  while (top >= 0) {
    j = stack[top];
    m = head[j];
    if (m == NULL_DIM) {
      --top;
      post[k++] = j;
    } else {
      head[j] = next[m];
      stack[++top] = m;
    }
  }

  return k;
}


/**
 * @brief Compute the elimination tree for symbolic cholesky factorization.
 *
 * A = LL^T
 *
 * @param nrows The number of rows of A.
 * @param rowptr The pointer to the start of each row of A.
 * @param rowind The column indexes stored in each row of A.
 * @param parent The parent of column element (output, nrows in length).
 * @param post The post ordering of the trees (output, nrows in length).
 */
static void __compute_elim_tree(
    dim_t const nrows,
    ind_t const * const rowptr,
    dim_t const * const rowind,
    dim_t * const parent,
    dim_t * const post)
{
  dim_t i, k, next;
  ind_t j;
  dim_t * wnext, * whead, * wstack, * scratch;

  scratch = dim_alloc(nrows*3);

  for (i=0;i<nrows;++i) {
    parent[i] = NULL_DIM;
    scratch[i] = NULL_DIM;
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      k = rowind[j];
      while (k != NULL_DIM && k < i) { /* stay in the lower triangle */
        next = scratch[k];
        scratch[k] = i;
        if (next == NULL_DIM) {
          parent[k] = i;
        }
        k = next;
      }
    }
  }

  /* post order the tree */
  whead = scratch;
  wnext = whead + nrows;
  wstack = wnext + nrows; 

  /* mark all whead empty */
  dim_set(whead,NULL_DIM,nrows);

  for (i=nrows;i>0;) {
    --i;
    if (parent[i] == NULL_DIM) {
      /* root vertex */
      continue; 
    }
    /* store parents pervious chain */
    wnext[i] = whead[parent[i]];
    /* add myself to parents chain */
    whead[parent[i]] = i;
  }

  k = 0;
  for (i=0;i<nrows;++i) {
    /* find tree roots */
    if (parent[i] == NULL_DIM) {
      k = __traverse_dfs(i,k,whead,wnext,post,wstack);
    }
  }
  
  /* free whead, wnext, and wstack */
  dl_free(scratch);
}


/**
 * @brief Compute the non-zero pattern of column i of the matrix L in cholesky
 * decomposition.
 *
 * A = LL^T
 *
 * @param i The column to compute.
 * @param nrows The number of rows in A.
 * @param rowptr The pointer to the start of each row of A.
 * @param rowind The column indexes stored in each row of A.
 * @param mrk The marking array from previous rows (lenght of nrows).
 * @param parent The parent array of the elimination tree.
 * @param symrow The row to populate with non-zero indices.
 *
 * @return The number of non-zeros in this row
 */
static dim_t __compute_colpattern(
    dim_t const i,
    dim_t const nrows,
    ind_t const * const rowptr,
    dim_t const * const rowind,
    int * const mrk,
    dim_t const * const parent,
    dim_t * const symrow)
{
  dim_t k, top, len;
  ind_t j;

  /* this comes from Tim Davis's book "Direct Methods for Sparse Linear 
   * Systems", page 53. */

  top = nrows;

  __cs_mark(mrk,i);

  /* go through the ith row of A */
  for (j=rowptr[i];j<rowptr[i+1];++j) {
    k = rowind[j];
    if (k > i) {
      /* only use upper triangular part of A */
      continue;
    }
    for (len=0;!__cs_marked(mrk,k);k=parent[k]) {
      /* traverse up the tree storing the path and marking vertices */
      symrow[len++] = k;
      __cs_mark(mrk,k);
    }
    /* push the path onto the stack */
    while (len > 0) {
      symrow[--top] = symrow[--len];
    }
  }
  /* unmark the excess rows */
  for (k=top;k<nrows;++k) {
    __cs_mark(mrk,symrow[k]);
  }
  __cs_mark(mrk,i);

  return top;
}


static dim_t __leaf(
    dim_t i,
    dim_t const k,
    dim_t const * const first,
    dim_t * const maxfirst,
    dim_t * const prevleaf,
    dim_t * const ancestor,
    dim_t * const jleaf)
{
  dim_t s, sparent, jprev;

  *jleaf = 0;

	if (i <= k || (maxfirst[i] != NULL_DIM && first[k] <= maxfirst[i])) {
    return -1;
  }

  maxfirst[i] = first[k];
  jprev = prevleaf[i];
  prevleaf[i] = k;

  if (jprev == NULL_DIM) {
    *jleaf = 1;
  } else {
    *jleaf = 2;
    /* this seems stupid and should just be an if statement */
    for (i=jprev;i!=ancestor[i];i=ancestor[i]);

    for (s=jprev;s!=i;s=sparent) {
      sparent = ancestor[s];
      ancestor[s] = i;
    }
  }

  return i;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void cholesky_rowcounts(
    dim_t const nrows,
    ind_t const * const rowptr,
    dim_t const * const rowind,
    dim_t * const counts)
{
  dim_t i, k, m, q, kleaf;
  ssize_t * delta;
  ind_t j;
  dim_t * ancestor, * maxfirst, * prevleaf, * first, * parent, * post;

  parent = dim_alloc(nrows);
  post = dim_alloc(nrows);

  __compute_elim_tree(nrows,rowptr,rowind,parent,post);

  ancestor = dim_alloc(nrows);
  maxfirst = dim_init_alloc(NULL_DIM,nrows);
  prevleaf = dim_init_alloc(NULL_DIM,nrows);
  first = dim_init_alloc(NULL_DIM,nrows);

  delta = ssize_alloc(nrows);

  /* compute first */
  for (i=0;i<nrows;++i) {
    k = post[i];
    delta[k] = (first[k] == NULL_DIM);
    while (k != NULL_DIM && first[k] == NULL_DIM) {
      first[k] = i;
      /* go up the tree */
      k = parent[k];
    }
  }

  dim_incset(ancestor,0,1,nrows);
  
  for (i=0;i<nrows;++i) {
    k = post[i];
    if (parent[k] != NULL_DIM) {
      /* decrement non-root vertices */
      --delta[parent[k]];
    }
    for (j=rowptr[k];j<rowptr[k+1];++j) {
      m = rowind[j];
      q = __leaf(m,k,first,maxfirst,prevleaf,ancestor,&kleaf);
      if (kleaf == 1) {
        /* [p,k] exits in L */
        ++delta[k];
      } else if (kleaf == 2) {
        /* [p,k] exits in L */
        ++delta[k];
        /* remove overlap */
        --delta[q];
      }
    }
    if (parent[k] != NULL_DIM) {
      ancestor[k] = parent[k];
    }
  }

  /* copy deltas over to counts */
  for (i=0;i<nrows;++i) {
    counts[i] = delta[i];
  }

  /* sum counts for each child */
  for (i=0;i<nrows;++i) {
    if (parent[i] != NULL_DIM) {
      counts[parent[i]] += counts[i];
    }
  }


  dl_free(delta);
  dl_free(ancestor);
  dl_free(maxfirst);
  dl_free(prevleaf);
  dl_free(first);
  dl_free(post);
  dl_free(parent);
}


void cholesky_colcounts(
    dim_t const nrows,
    ind_t const * const rowptr,
    dim_t const * const rowind,
    dim_t * const counts)
{
  dim_t i, k, m, q, kleaf;
  ind_t j;
  dim_t * ancestor, * maxfirst, * prevleaf, * first, * level, * parent, * post;

  parent = dim_alloc(nrows);
  post = dim_alloc(nrows);

  __compute_elim_tree(nrows,rowptr,rowind,parent,post);

  ancestor = dim_alloc(nrows);
  maxfirst = dim_init_alloc(NULL_DIM,nrows);
  prevleaf = dim_init_alloc(NULL_DIM,nrows);
  first = dim_alloc(nrows);
  level = dim_alloc(nrows);

  dim_incset(ancestor,0,1,nrows);

  __firstdesc(nrows,parent,post,first,level);

  /* assume its at least a diagonal matrix */
  dim_set(counts,1,nrows);

  for (i=0;i<nrows;++i) {
    k = post[i];
    for (j=rowptr[k];j<rowptr[k+1];++j) {
      m = rowind[j];
      q = __leaf(m,k,first,maxfirst,prevleaf,ancestor,&kleaf);
      if (kleaf) {
        counts[m] += level[k] - level[q];
      }
    }
    if (parent[k] != NULL_DIM) {
      ancestor[k] = parent[k];
    }
  }

  dl_free(ancestor);
  dl_free(maxfirst);
  dl_free(prevleaf);
  dl_free(first);
  dl_free(level);
  dl_free(parent);
  dl_free(post);
}


void cholesky_symbolic(
    dim_t const nrows,
    ind_t const * const rowptr,
    dim_t const * const rowind,
    ind_t ** const r_symptr,
    dim_t ** const r_symind)
{
  #ifdef XXX
  dim_t i, k, m, nr;
  ind_t j, nnz, maxnnz;
  ind_t * symptr;
  dim_t * row, * mrk, * min, * symind;

  symptr = ind_alloc(nrows+1);

  /* guess how many nnz's in L */
  maxnnz = NNZ_FACTOR*(nrows*rowptr[nrows]); 
  symind = dim_alloc(maxnnz);

  /* Here's the following algorithm:
   *
   * mrk[0] = 0 for all i;
   * for i = 1 to n:
   *   row = A[i];
   *   for all k such that mrk[j] == i:
   *     row = (row \union L[j]) \setminus {j};
   *   mrk[i] = min(row \setminus {i});
   */

  mrk = dim_init_alloc(NULL_DIM,nrows);

  dl_free(mrk);

  *r_symptr = symptr;
  *r_symind = symind;
  #endif
}





#endif
