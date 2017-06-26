/**
 * @file matrix.h
 * @brief Functions for creating and manipulating matrices
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-04-30
 */




#ifndef GOOSEBERRY_MATRIX_H
#define GOOSEBERRY_MATRIX_H




#include "base.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef enum matrix_type_t {
  MATRIX_TYPE_DENSE_VECTOR,
  MATRIX_TYPE_DENSE_MATRIX,
  MATRIX_TYPE_SPARSE_VECTOR,
  MATRIX_TYPE_SPARSE_MATRIX
} matrix_type_t;


typedef struct matrix_t {
  matrix_type_t type;
  dim_t nrows, ncols;
  ind_t * rowptr, * colptr;
  dim_t * rowind, * colind;
  real_t * rowval, * colval;
  dim_t * rdist, * cdist;
  dim_t nrdist, ncdist;
} matrix_t;




/******************************************************************************
* MEMORY FUNCTION PROTOTYPES **************************************************
******************************************************************************/


#define DLMEM_PREFIX matrix
#define DLMEM_TYPE_T matrix_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


int matrix_init(
    int type, 
    dim_t nrows, 
    dim_t ncols,
    ind_t nnz,
    matrix_t * mat);


int matrix_zero(
    matrix_t * mat);


int matrix_identify(
    matrix_t * mat);


int matrix_sparsify(
    matrix_t * mat);


int matrix_densify(
    matrix_t * mat);


int matrix_transpose(
    matrix_t * mat);


int matrix_permute(
    matrix_t * mat,
    const dim_t * rperm,
    const dim_t * cperm);


int matrix_unpermute(
    matrix_t * mat,
    const dim_t * rperm,
    const dim_t * cperm);


int matrix_buildindex(
    matrix_t * mat);


int matrix_free(
    matrix_t * mat);




#endif
