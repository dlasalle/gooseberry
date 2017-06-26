/**
 * @file cholesky.h
 * @brief CHolesky factorization functions prototypes.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2015-02-02
 */




#ifndef GOOSEBERRY_CHOLESKY_H
#define GOOSEBERRY_CHOLESKY_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define cholesky_rowcounts __gooseberry_cholesky_rowcounts
/**
 * @brief Count the number of non-zeroes per row in the cholesky decomposition.
 * A = LL^T
 *
 * @param nrows The number of rows of the matrix A.
 * @param rowptr The pointer to the start of each row in A.
 * @param rowind The column index for each element in the rows A.
 * @param counts The number of non-zeroes in each row (output, of length 
 *   nrows).
 */
void cholesky_rowcounts(
    dim_t nrows,
    ind_t const * rowptr,
    dim_t const * rowind,
    dim_t * counts);


#define cholesky_colcounts __gooseberry_cholesky_colcounts
/**
 * @brief Count the number of non-zeroes per column in the cholesky 
 * decomposition.
 * A = LL^T
 *
 * @param nrows The number of rows of the matrix A.
 * @param rowptr The pointer to the start of each row in A.
 * @param rowind The column index for each element in the rows A.
 * @param counts The number of non-zeroes in each column (output, of length 
 *   nrows).
 */
void cholesky_colcounts(
    dim_t nrows,
    ind_t const * rowptr,
    dim_t const * rowind,
    dim_t * counts);


#define cholesky_symbolic __gooseberry_cholesky_symbolic
/**
 * @brief Perform a symbolic cholesky decomposition on a symetric matrix.
 * A = LL^T
 *
 * @param nrows The number of rows of the matrix A.
 * @param rowptr The pointer to the start of each row in A.
 * @param rowind The column index for each element in the rows A.
 * @param r_symptr A reference to the pointer to the start of each row in L 
 *   (output). 
 * @param r_symind A reference to the column index for each element in the 
 *   rows of L (output).
 */
void cholesky_symbolic(
    dim_t nrows,
    ind_t const * rowptr,
    dim_t const * rowind,
    ind_t ** r_symptr,
    dim_t ** r_symind);




#endif
