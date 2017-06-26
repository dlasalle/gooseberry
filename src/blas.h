/**
 * @file dlblas_headers.h
 * @brief Basic linear algebra prototypes
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-09-29
 */




#ifndef GOOSEBERRY_BLAS_H
#define GOOSEBERRY_BLAS_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define blas_scale __goseeberry_blas_scale
/**
 * @brief Scale a matrix/vector by a scalar: A = sA
 * This function works on sparse and dense matrices.
 *
 * @param n The total length of the rowval array
 * @param rowval Tthe values to scale
 * @param s The scalar
 *
 * @return GOOSEBERRY_SUCCESS unless the input is invalid
 */
int blas_scale(
    ind_t n,
    real_t * rowval,
    real_t s);


#define blas_negate __gooseberry_blas_negate
/**
 * @brief Negate all of the values in a matrix/vector: A = -A
 * This functions works on both sparse and dense matrices.
 *
 * @param n The total length of the rowval array
 * @param rowval The values to negate
 *
 * @return GOOSEBERRY_SUCCESS unless there is invalid input 
 */
int blas_negate(
    ind_t n,
    real_t * rowval);


#define blas_add __gooseberry_blas_add
/**
 * @brief Add two dense matrices/vectors: A + B = C
 *
 * @param anrows The number of rows
 * @param ancols The number columns
 * @param arowval The values of matrix A in row major order
 * @param browval The values of matrix B in row major order
 * @param crowval The values of the output matrix C in row major order
 *
 * @return GOOSEBERRY_SUCCESS unless the input is invalid
 */
int blas_add(
    dim_t anrows,
    dim_t ancols,
    const real_t * arowval,
    const real_t * browval,
    real_t * crowval);


#define blas_addsp __gooseberry_blas_addsp
/**
 * @brief Add a dense matrix/vector and sparse matrix: A + B = C
 *
 * @param anrows The number of rows
 * @param ancols The number of columns
 * @param arowval The values of matrix A in row major order
 * @param browptr The starting index of each row of B 
 * @param browind The column index of each nonzero value of B
 * @param browval The nonzero values of B
 * @param crowval The values of the output matrix C in row major order
 *
 * @return GOOSEBERRY_SUCCESS unless the input is invalid
 */
int blas_addsp(
    dim_t anrows,
    dim_t ancols,
    const real_t * arowval,
    const ind_t * browptr,
    const dim_t * browind,
    const dim_t * browval,
    real_t * crowval);


#define blas_sub __gooseberry_blas_sub
/**
 * @brief Subtract a dense matrix/vector from a dense matrix/vector: A - B = C
 *
 * @param anrows The number of rows
 * @param ancols The number of columns
 * @param arowval The values of matrix A in row major order
 * @param browval The values of matrix B in row major order
 * @param crowval The values of the output matrix C in row major order
 *
 * @return GOOSEBERRY_SUCCESS unless the input is invalid
 */
int blas_sub(
    dim_t anrows,
    dim_t ancols,
    const real_t * arowval,
    const real_t * browval,
    real_t * crowval);


#define blas_subsp __goseeberry_blas_subsp
/**
 * @brief Subtract a sparse matrix from a dense matrix: A - B = C
 *
 * @param anrows The number of rows
 * @param ancols The number of columns
 * @param arowval The values of matrix A in row major order
 * @param browptr The starting index of each row of B 
 * @param browind The column index of each nonzero value of B
 * @param browval The nonzero values of B
 * @param crowval The values of the output matrix C in row major order
 *
 * @return GOOSEBERRY_SUCCESS unless the input is invalid
 */
int blas_subsp(
    dim_t anrows,
    dim_t ancols,
    const real_t * arowval,
    const ind_t * browptr,
    const dim_t * browind,
    const real_t * browval,
    real_t * crowval);


#define blas_spsub __gooseberry_blas_spsub
/**
 * @brief Subtract a dense matrix from a sparse matrix: A - B = C
 *
 * @param anrows The number of rows
 * @param ancols The number of columns
 * @param arowptr The starting index of each row of A
 * @param arowind The column index of each nonzero value of A
 * @param arowval The nonzero values of A
 * @param browval The values of matrix B in row major order
 * @param crowval The values of the output matrix C in row major order
 *
 * @return GOOSEBERRY_SUCCESS unless the input is invalid
 */
int blas_spsub(
    dim_t anrows,
    dim_t ancols,
    const ind_t * arowptr,
    const dim_t * arowind,
    const real_t * arowval,
    const real_t * browval,
    real_t * crowval);


#define blas_add_scalar __gooseberry_blas_add_scalar
/**
 * @brief Add a scalar to all elements of a matrix vector
 *
 * @param n The number of elements in the matrix/vector. For dense structures
 * this is nrows*ncols, for sparse structures, this is nnz.
 * @param val The elements of the matrix/vector
 * @param s The scalar to add
 *
 * @return GOOSEBERRY_SUCCESS
 */
int blas_add_scalar(
    ind_t n,
    real_t * val,
    real_t s);


#define blas_dot __gooseberry_blas_dot
/**
 * @brief Performs the dot product of two dense vectors: a^t b = s
 *
 * @param n The number of elements in both vectors
 * @param arowval The values of vector a in row major order 
 * @param browval The values of vector b in column major order
 *
 * @return The resulting dot product
 */
real_t blas_dot(
    dim_t n,
    const real_t * arowval,
    const real_t * browval);


#define blas_spdot __gooseberry_blas_spdot 
/**
 * @brief Performs the dot product of a sparse vector a and a dense vector b:
 *        a^t b = s
 *
 * @param annz The number of non-zero elements in a
 * @param arowind The column index of each non-zero element in a
 * @param arowval The value of each non-zero element in a
 * @param bcolval The values of vector b in column major order
 *
 * @return The resulting dot product
 */
real_t blas_spdot(
    dim_t annz,
    const real_t * arowind, 
    const real_t * arowval,
    const real_t * bcolval);


#define blas_spdotsp __gooseberry_blas_spdotsp
/**
 * @brief Performs the dot product of two sparse vectors: a^t b = s 
 *
 * @param annz The number of non-zero elements in a
 * @param bnnz The number of non-zero elements in b
 * @param arowind The column index of each non-zero element in a
 * @param arowval The value of each non-zero element in a
 * @param bcolind The column index of each non-zero elmeent in b
 * @param bcolval The value of each non-zero element in b
 *
 * @return The resulting dot prodcut
 */
real_t blas_spdotsp(
    dim_t annz,
    dim_t bnnz,
    const real_t * arowind, 
    const real_t * arowval,
    const real_t * bcolind, 
    const real_t * bcolval);


#define blas_mult __gooseberry_blas_mult
/**
 * @brief Performs dense matrix/vector dense matrix/vector multiplication: 
 *        A x B = C
 *
 * @param anrows The number of rows of the matrix A
 * @param ancols The number of columns of the matrix A (should also be the
 * number of rows of matrix B).
 * @param arowval The values matrix A in row major format
 * @param bcolval The values matrix B in column major format (B^t)
 * @param crowval The values of the output matrix C in row major format
 * @param blocks The starting index of rows for each block of matrix A and C
 * (leave NULL for automatic assignment).
 * @param nblocks The number of blocks of A and C
 *
 * @return GOOSEBERRY_SUCCESS unless there is a problem with the input values 
 */
int blas_mult(
    dim_t anrows,
    dim_t ancols,
    dim_t bncols,
    const real_t * arowval,
    const real_t * bcolval,
    real_t * crowval,
    const dim_t * blocks,
    dim_t nblocks);


#define blas_spmult __gooseberry_blas_spmult
/**
 * @brief Performs sparse matrix dense matrix/vector multiplication:
 *        A x B = C
 *
 * @param anrows The number of rows of the matrix A
 * @param ancols The number of columns of the matrix A (should also be the same
 * number of rows of matrix B).
 * @param arowptr The starting index of each row of A 
 * @param arowind The column index of each nonzero value of A
 * @param arowval The nonzero values of A
 * @param bcolval The values of matrix B in column major format (B^t)
 * @param crowval The values of the output matrix C in row major format
 * @param blocks The starting index of rows for each block of matrix A and C.
 * (leave NULL for automatic assignment).
 * @param nblocks The number of blocks of A and C
 *
 * @return GOOSEBERRY_SUCCESS unless tehre is a problem with the input values 
 */
int blas_spmult(
    dim_t anrows,
    dim_t ancols,
    dim_t bncols, 
    const ind_t * arowptr,
    const dim_t * arowind,
    const real_t * arowval,
    const real_t * bcolval,
    real_t * crowval,
    const dim_t * blocks,
    dim_t nblocks);


#define blas_spmultsp __gooseberry_blas_spmultsp
/**
 * @brief Performs sparse matrix sparse matrix/vector multiplication:
 *        A x B = C
 *
 * @param anrows The number of rows of the matrix A
 * @param ancols The number of columns of the matrix A (should also be the same
 * number of rows of matrix B).
 * @param arowptr The starting index of each row of A 
 * @param arowind The column index of each nonzero value of A
 * @param arowval The nonzero values of A
 * @param browptr The starting index of each row of B 
 * @param browind The column index of each nonzero value of B
 * @param browval The nonzero values of B
 * @param crowval The values of the output matrix C in row major format
 * @param blocks The starting index of rows for each block of matrix A and C
 * (leave NULL for automatic assignment).
 * @param nblocks The number of blocks of A and C
 *
 * @return GOOSEBERRY_SUCCESS unless tehre is a problem with the input values 
 */
int blas_spmultsp(
    dim_t anrows,
    dim_t ancols,
    dim_t bncols,
    const ind_t * arowptr,
    const dim_t * arowind,
    const real_t * arowval,
    const ind_t * browptr,
    const dim_t * browind,
    const real_t * browval,
    real_t * crowval,
    const dim_t * blocks,
    dim_t nblocks);


#define blas_spmultsp_sp __gooseberry_spmultsp_sp
/**
 * @brief Performs sparse matrix sparse matrix/vector multiplication:
 *        A x B = C
 *
 * @param anrows The number of rows of the matrix A
 * @param ancols The number of columns of the matrix A (should also be the same
 * number of rows of matrix B).
 * @param arowptr The starting index of each row of A 
 * @param arowind The column index of each nonzero value of A
 * @param arowval The nonzero values of A
 * @param browptr The starting index of each row of B 
 * @param browind The column index of each nonzero value of B
 * @param browval The nonzero values of B
 * @param crowptr The starting index of each row of C 
 * @param crowind The column index of each nonzero value of C
 * @param crowval The nonzero values of C
 * @param blocks The starting index of rows for each block of matrix A and C
 * (leave NULL for automatic assignment).
 * @param nblocks The number of blocks of A and C
 *
 * @return GOOSEBERRY_SUCCESS unless tehre is a problem with the input values 
 */
int blas_spmultsp_sp(
    dim_t anrows,
    dim_t ancols,
    dim_t bncols,
    const ind_t * arowptr,
    const dim_t * arowind,
    const real_t * arowval,
    const ind_t * browptr,
    const dim_t * browind,
    const real_t * browval,
    ind_t ** crowptr,
    dim_t ** crowind,
    real_t ** crowval,
    const dim_t * blocks,
    dim_t nblocks);




#endif
