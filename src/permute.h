/**
 * @file permute.h
 * @brief Functions for a re-ordering a matrix 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-06-05
 */




#ifndef GOOSEBERRY_PERMUTE_H
#define GOOSEBERRY_PERMUTE_H




#include "base.h"





/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


int permute_sparse(
    dim_t nrows,
    dim_t ncols,
    ind_t * rowptr,
    dim_t * rowind,
    real_t * rowval,
    const dim_t * rowperm,
    const dim_t * colperm);


int permute_cuthillmckee(
    dim_t nrows,
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval,
    const dim_t * blocks,
    dim_t nblocks,
    dim_t * perm);


int permute_dense(
    dim_t nrows,
    dim_t ncols,
    real_t * rowval,
    const dim_t * rowperm,
    const dim_t * colperm);


int permute_dense_rev(
    dim_t nrows,
    dim_t ncols,
    real_t * rowval,
    const dim_t * rowperm,
    const dim_t * colperm);




#endif
