/**
 * @file sgd.h
 * @brief Function prototypes for SGD
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-04-26
 */




#ifndef GOOSEBERRY_SGD_H
#define GOOSEBERRY_SGD_H




#include "matrix.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


int sgd(
    matrix_t * mat, 
    matrix_t * v, 
    matrix_t * u,
    real_t rate,
    real_t lambda,
    size_t niter);





#endif
