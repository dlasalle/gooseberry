/**
 * @file transform.h
 * @brief Matrix transformation functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-06-10
 */




#ifndef GOOSEBERRY_TRANSFORM_H
#define GOOSEBERRY_TRANSFORM_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define transform_symmetrify_sparse \
    __GOOSEBERRY_transform_symmetrify_sparse
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
int transform_symmetrify_sparse(
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval,
    ind_t ** r_rowptr,
    dim_t ** r_rowind,
    real_t ** r_rowval);


#define transform_debipartify_sparse \
    __GOOSEBERRY_transform_debipartify_sparse
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
int transform_debipartify_sparse(
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval,
    ind_t ** r_rowptr,
    dim_t ** r_rowind,
    real_t ** r_rowval);


#define transform_rowsplit_sparse \
    __GOOSEBERRY_transform_rowsplit_sparse
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
int transform_rowsplit_sparse(
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


#define transform_colsplit_sparse \
    __GOOSEBERRY_transform_colsplit_sparse
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
int transform_colsplit_sparse(
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


#define transform_transpose_sparse \
    __GOOSEBERRY_transform_transpose_sparse
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
int transform_transpose_sparse(
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval,
    ind_t ** r_trowptr,
    dim_t ** r_trowind,
    real_t ** r_trowval);


#endif
