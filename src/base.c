/**
 * @file base.c
 * @brief Basic types for Gooseberry (Falls)
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-05-01
 */





#ifndef GOOSEBERRY_BASE_C
#define GOOSEBERRY_BASE_C




#include "base.h"




/******************************************************************************
* DOMLIB FUNCTIONS ************************************************************
******************************************************************************/


#ifndef GOOSEBERRY_TYPES_DEFINED

/* real_t */
#define DLMEM_PREFIX real
#define DLMEM_TYPE_T real_t
#define DLMEM_DLTYPE DLTYPE_FLOAT
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLMATH_PREFIX real 
#define DLMATH_TYPE_T real_t
#define DLMATH_DLTYPE DLTYPE_FLOAT
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T

#define DLRAND_PREFIX real
#define DLRAND_TYPE_T real_t
#define DLRAND_DLTYPE DLTYPE_FLOAT
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T


/* dim_t */
#define DLMEM_PREFIX dim
#define DLMEM_TYPE_T dim_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLRAND_PREFIX dim
#define DLRAND_TYPE_T dim_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T

#define DLMATH_PREFIX dim
#define DLMATH_TYPE_T dim_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T

#define DLSTATS_PREFIX dim
#define DLSTATS_TYPE_T dim_t
#include "dlstats_funcs.h"
#undef DLSTATS_PREFIX
#undef DLSTATS_TYPE_T

#define DLSORTKV_PREFIX dd
#define DLSORTKV_KEY_T dim_t
#define DLSORTKV_VAL_T dim_t
#define DLSORTKV_SIZE_T dim_t
#include "dlsortkv_funcs.h"
#undef DLSORTKV_SIZE_T
#undef DLSORTKV_PREFIX
#undef DLSORTKV_KEY_T
#undef DLSORTKV_VAL_T


/* ind_t */
#define DLMEM_PREFIX ind
#define DLMEM_TYPE_T ind_t
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_funcs.h"
#undef DLMEM_DLTYPE
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T

#define DLRAND_PREFIX ind
#define DLRAND_TYPE_T ind_t
#define DLRAND_DLTYPE DLTYPE_INTEGRAL
#include "dlrand_funcs.h"
#undef DLRAND_DLTYPE
#undef DLRAND_PREFIX
#undef DLRAND_TYPE_T

#define DLMATH_PREFIX ind
#define DLMATH_TYPE_T ind_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#include "dlmath_funcs.h"
#undef DLMATH_DLTYPE
#undef DLMATH_PREFIX
#undef DLMATH_TYPE_T

#endif /* GOOSEBERRY_TYPES_DEFINED */





#endif
