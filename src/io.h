/**
 * @file io.h
 * @brief I/O function prototypes for matrices (and vectors).
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-05-01
 */





#ifndef GOOSEBERRY_IO_H
#define GOOSEBERRY_IO_H



#include "base.h"





/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


int read_sparse_matrix(
    int type, 
    const char * file,
    dim_t * nrows, 
    dim_t * ncols,
    ind_t ** rowptr,
    dim_t ** rowind,
    real_t ** rowval);


int read_dense_matrix(
    int type, 
    const char * file,
    dim_t * nrows, 
    dim_t * ncols,
    real_t ** rowval);


int read_labels(
    const char * filename,
    dim_t * nrows,
    dim_t ** labels);


int write_dense_matrix(
    int type, 
    const char * file,
    dim_t nrows, 
    dim_t ncols,
    const real_t * rowval);


int write_sparse_matrix(
    int type, 
    const char * file,
    dim_t nrows, 
    dim_t ncols,
    const ind_t * rowptr,
    const dim_t * rowind,
    const real_t * rowval);




#endif
