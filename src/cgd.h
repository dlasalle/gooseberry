/**
 * @file cgd.h
 * @brief Conjugate gradient desecent function prototypes
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-05-06
 */




#ifndef GOOSEBERRY_CGD_H
#define GOOSEBERRY_CGD_H




#include "matrix.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


int cgd(
    matrix_t * a, 
    matrix_t * b,
    matrix_t * x,
    real_t error,
    size_t niter);




#endif
