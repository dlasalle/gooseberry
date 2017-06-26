/**
 * @file gooseberry.h
 * @brief The Gooseberry (Falls) library header
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-04-30
 */




#ifndef GOOSEBERRY_H
#define GOOSEBERRY_H




#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define GOOSEBERRY_VER_MAJOR 0
#define GOOSEBERRY_VER_MINOR 1
#define GOOSEBERRY_VER_SUBMINOR 6




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


#ifndef BOWSTRING_TYPES_DEFINED

#ifdef GOOSEBERRY_SIGNED_TYPES
typedef int64_t gooseberry_int64_t;
#define PF_BSINT64_T "%"PRIi64
typedef int32_t gooseberry_int32_t;
#define PF_BSINT32_T "%"PRIi32
#else
typedef uint64_t gooseberry_int64_t;
#define PF_BSINT64_T "%"PRIu64
typedef uint32_t gooseberry_int32_t;
#define PF_BSINT32_T "%"PRIu32
#endif /* GOOSEBERRY_SIGNED_TYPES */

#ifdef GOOSEBERRY_64BIT_DIMENSIONS
typedef gooseberry_int64_t dim_t;
#define PF_DIM_T PF_BSINT64_T
#else
typedef gooseberry_int32_t dim_t;
#define PF_DIM_T PF_BSINT32_T
#endif /* GOOSEBERRY_64BIT_VERTICES */

#ifdef GOOSEBERRY_64BIT_INDICES
typedef gooseberry_int64_t ind_t;
#define PF_IND_T PF_BSINT64_T
#else
typedef gooseberry_int32_t ind_t;
#define PF_IND_T PF_BSINT32_T
#endif /* GOOSEBERRY_64BIT_EDGES */

#ifdef GOOSEBERRY_DOUBLE_REALS
typedef double real_t;
#define PF_REAL_T "%lf"
#else
typedef float real_t;
#define PF_REAL_T "%f"
#endif /* GOOSEBERRY_DOUBLE_REALS */
#define PF_MINREAL_T "%g"

#endif /* GOOSEBERRY_TYPES_DEFINED */


typedef enum gooseberry_error_t {
  GOOSEBERRY_SUCCESS=1,
  GOOSEBERRY_ERROR_UNIMPLEMENTED,
  GOOSEBERRY_ERROR_INVALIDINPUT,
  GOOSEBERRY_ERROR_READ,
  GOOSEBERRY_ERROR_WRITE,
  GOOSEBERRY_ERROR_FILENOTFOUND,
  GOOSEBERRY_ERROR_PERMISSIONDENIED,
  GOOSEBERRY_ERROR_UNKNOWN
} gooseberry_error_t;


typedef enum gooseberry_format_t {
  GOOSEBERRY_FORMAT_RAW,
  GOOSEBERRY_FORMAT_GRID,
  GOOSEBERRY_FORMAT_CSR,
  GOOSEBERRY_FORMAT_POINT,
  GOOSEBERRY_FORMAT_SVM,
  GOOSEBERRY_FORMAT_GRAPH,
  GOOSEBERRY_FORMAT_CLU,
  GOOSEBERRY_FORMAT_MATRIXMARKET
} gooseberry_format_t;




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#ifdef __cplusplus
extern 'C' {
#endif


/* BLAS **********************************************************************/


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
int gooseberry_scale(
    ind_t n,
    real_t * val,
    real_t s);


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
int gooseberry_add_scalar(
    ind_t n,
    real_t * val,
    real_t s);


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
int gooseberry_spmult(
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


int gooseberry_spmultsp(
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


int gooseberry_spmultsp_sp(
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


/* TRANSFORM *****************************************************************/


/**
 * @brief Split a sparse matrix along its rows
 *
 * The resulting matrices are accessible via:
 *   for (p=0;p<nparts;++p) {
 *     for (i=dist[p];i<dist[p+1];++i) {
 *       row = i-dist[p];
 *       for (j=(*r_rowptr)[i];j<(*r_rowptr)[i+1];++j) {
 *         col = (*r_rowind)[j];
 *         val = (*r_rowval)[j];
 *       }
 *     }
 *   }
 *
 * @param nrows The number of rows
 * @param ncols The number of columns
 * @param rowptr The pointer to the row adjacency lists
 * @param rowind The column indexes for each row
 * @param rowval The values for each row
 * @param nparts The number of partitions
 * @param dist The row distribution
 * @param map The mapping of rows to partitions
 * @param r_rowptr The output segmented row adjacncy lists
 * @param r_rowind The output segmented column indexes
 * @param r_rowval The output segmented row values
 *
 * @return GOOSEBERRY_SUCCESS on success
 */
int gooseberry_rowsplit_sparse(
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval,
    dim_t nparts,
    const dim_t * dist,
    const dim_t * map,
    ind_t ** r_rowptr,
    dim_t ** r_rowind,
    real_t ** r_rowval);


/**
 * @brief Split a sparse matrix along its columns
 *
 * The resulting matrices are accessible via:
 *   for (p=0;p<nparts;++p) {
 *     for (i=dist[p];i<dist[p+1];++i) {
 *       row = i-dist[p];
 *       for (j=(*r_rowptr)[i];j<(*r_rowptr)[i+1];++j) {
 *         col = (*r_rowind)[j];
 *         val = (*r_rowval)[j];
 *       }
 *     }
 *   }
 *
 * @param nrows The number of rows
 * @param ncols The number of columns
 * @param rowptr The pointer to the row adjacency lists
 * @param rowind The column indexes for each row
 * @param rowval The values for each row
 * @param nparts The number of partitions
 * @param dist The row distribution
 * @param map The mapping of rows to partitions
 * @param r_rowptr The output segmented row adjacncy lists
 * @param r_rowind The output segmented column indexes
 * @param r_rowval The output segmented row values
 *
 * @return GOOSEBERRY_SUCCESS on success
 */
int gooseberry_colsplit_sparse(
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval,
    dim_t nparts,
    const dim_t * dist,
    const dim_t * map,
    ind_t ** r_rowptr,
    dim_t ** r_rowind,
    real_t ** r_rowval);


/**
 * @brief Create a symmetric matrix from a non-symmetric rectangular matrix via
 * doing an element-wise addition of the matrix and its transpose.
 * B = A + A^T
 *
 * @param nrows The number of rows in the matrix
 * @param ncols The number columns in the matrix
 * @param rowptr The pointer to the row adjacency lists
 * @param rowind The column indexes for each row
 * @param rowval The values for each row
 * @param r_rowptr A reference to the output pointer to the row adjacency lists
 * @param r_rowind A reference to the output pointer to the column indexes
 * @param r_rowval A reference to the output pointer to the row values
 *
 * @return GOOSEBERRY_SUCCESS on success
 */
int gooseberry_symmetrify_sparse(
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval,
    ind_t ** r_rowptr,
    dim_t ** r_rowind,
    real_t ** r_rowval);


/**
 * @brief Create a symmetric matrix from a non-symmetric rectangular matrix via
 * treating it as bipartite graph and converting it to an undirected generic
 * graph. 
 * B = [0, A; A^T, 0]
 *
 * @param nrows The number of rows in the matrix
 * @param ncols The number columns in the matrix
 * @param rowptr The pointer to the row adjacency lists
 * @param rowind The column indexes for each row
 * @param rowval The values for each row
 * @param r_rowptr The output pointer to the row adjacency lists
 * @param r_rowind The output pointer to the column indexes
 * @param r_rowval The output pointer to the row values
 *
 * @return GOOSEBERRY_SUCCESS on success
 */
int gooseberry_debipartify_sparse(
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval,
    ind_t ** r_rowptr,
    dim_t ** r_rowind,
    real_t ** r_rowval);


/**
 * @brief Transpose a sparse matrix: B = A^T
 *
 * @param nrows Number of rows
 * @param ncols Number of columns
 * @param rowptr The pointer to the row adjacency lists
 * @param rowind The column indexes for each row
 * @param rowval The values for each row
 * @param r_trowptr The output pointer to the row adjacency lists
 * @param r_trowind The output pointer to the column indexes
 * @param r_trowval The output pointer to the row values
 *
 * @return GOOSEBERRY_SUCCESS on success
 */
int gooseberry_transpose_sparse(
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval,
    ind_t ** r_trowptr,
    dim_t ** r_trowind,
    real_t ** r_trowval);


/* IO ************************************************************************/


int gooseberry_read_sparse_matrix(
    int type, 
    const char * file, 
    dim_t * nrows,
    dim_t * ncols,
    ind_t ** rowptr,
    dim_t ** rowind,
    real_t ** rowval);


int gooseberry_read_dense_matrix(
    int type, 
    const char * file, 
    dim_t * nrows,
    dim_t * ncols,
    real_t ** rowval);


int gooseberry_write_dense_matrix(
    int type, 
    const char * file, 
    dim_t nrows,
    dim_t ncols,
    const real_t * rowval);


int gooseberry_write_sparse_matrix(
    int type, 
    const char * file, 
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval);


int gooseberry_read_labels(
    const char * file,
    dim_t * nrows,
    dim_t ** labels);


int gooseberry_read_partition(
    const char * const file,
    dim_t * const r_nrows,
    dim_t * const r_nparts,
    dim_t ** const r_map,
    dim_t ** const r_perm,
    dim_t ** const r_dist);


#ifdef __cplusplus
}
#endif

#endif
