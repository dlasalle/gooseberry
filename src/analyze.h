/**
 * @file analyze.h
 * @brief Function prototypes for determining matrix statistics
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-06-28
 */




#ifndef GOOSEBERRY_ANALYZE_H
#define GOOSEBERRY_ANALYZE_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define analyze_min_active_set __gooseberry_min_active_set
/**
 * @brief This function deterimnes the size of the Minimum Active Set of
 * elements in the matrix. That is, the minimum number of elements that must be
 * able to be stored in cache of a vector if we traverse the rows, so that no
 * element is loaded twice.
 *
 * @param nrows The number of rows of the matrix
 * @param ncols The number of columns of the matrix
 * @param rowptr The index in the adjacency lists of the start of each row
 * @param rowind The column index of each element
 * @param rowval The value of each element
 * @param blocks The row-wise blocking of the matrix (can be NULL)
 * @param nblocks The number of blocks
 * @param mer The output array of size nblocks (1 if block is NULL), speficying
 * the MAS for each block.
 *
 * @return GOOSEBERRY_SUCCESS if successful. 
 */
int analyze_min_active_set(
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval,
    const dim_t * blocks,
    dim_t nblocks,
    dim_t * mas);


#define analyze_cholesky __gooseberry_analyze_cholesky
/**
 * @brief Compute the number of non-zeroes in L and the number of operations
 * required for the factorization.
 *
 * A = LL^T
 *
 * @param nrows The number of rows in A.
 * @param rowptr The pointer to the start of each row of A.
 * @param rowind The column indexes of each elment.
 * @param r_nnz A reference to the number of non-zeroes in L.
 * @param r_nops A reference to the number of operationes requreid to compute
 *   L.
 */
void analyze_cholesky(
    dim_t nrows,
    ind_t const * rowptr,
    dim_t const * rowind,
    ind_t * r_nnz,
    double * r_nops);




#endif
