/**
 * @file base.h
 * @brief Basic types for Gooseberry
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-05-01
 */




#ifndef GOOSEBERRY_BASE_H
#define GOOSEBERRY_BASE_H




#include <domlib.h>
#include <stdlib.h> 
#include <stdio.h>
#include <stdint.h>
#ifndef NO_OMP
#include <omp.h>
#endif




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


/* Define the types for internal use if they are externally defined.
 * Otherwise gooseberry.h will provide definitions for us. */
#ifdef BOWSTRING_TYPES_DEFINED

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


#include <gooseberry.h>




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


/* real_t */
#define DLMEM_PREFIX real
#define DLMEM_TYPE_T real_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLMATH_PREFIX real
#define DLMATH_TYPE_T real_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T

#define DLRAND_PREFIX real
#define DLRAND_TYPE_T real_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


/* dim_t */
#define DLMEM_PREFIX dim
#define DLMEM_TYPE_T dim_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLRAND_PREFIX dim
#define DLRAND_TYPE_T dim_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T

#define DLMATH_PREFIX dim
#define DLMATH_TYPE_T dim_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T

#define DLSTATS_PREFIX dim
#define DLSTATS_TYPE_T dim_t
#include "dlstats_headers.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T

#define DLSORTKV_PREFIX dd
#define DLSORTKV_KEY_T dim_t
#define DLSORTKV_VAL_T dim_t
#define DLSORTKV_SIZE_T dim_t
#include "dlsortkv_headers.h"
#undef DLSORTKV_SIZE_T
#undef DLSORTKV_PREFIX
#undef DLSORTKV_KEY_T
#undef DLSORTKV_VAL_T


/* ind_t */
#define DLMEM_PREFIX ind
#define DLMEM_TYPE_T ind_t
#include "dlmem_headers.h"
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLRAND_PREFIX ind
#define DLRAND_TYPE_T ind_t
#include "dlrand_headers.h"
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T

#define DLMATH_PREFIX ind
#define DLMATH_TYPE_T ind_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_headers.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T




/******************************************************************************
* CONTSANTS *******************************************************************
******************************************************************************/


static const ind_t NULL_IND = (ind_t)-1;
static const dim_t NULL_DIM = (dim_t)-1;
static const real_t NULL_REAL = (real_t)-INFINITY;

/* openmp scheduling */
static const size_t OMP_SMALL_BLOCK = 16;
static const size_t OMP_BIG_BLOCK = 256;




#endif
