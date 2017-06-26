/**
 * @file mpi_spmultsp.c
 * @brief Perform a sparse matrix matrix multiplication via MPI
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-05-13
 */




#ifndef GOOSEBERRY_MPI_SPMSPM_C_C
#define GOOSEBERRY_MPI_SPMSPM_C_C




#include <mpi.h>
#include <domlib.h>
#include <gooseberry.h>




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define MPI_DIM_T MPI_UNSIGNED
#define MPI_REAL_T MPI_FLOAT


#define MPI_Alltoallv(sb,sc,so,sd,rb,rc,ro,rd,c) \
  MPI_Alltoallv(sb,(int*)sc,(int*)so,sd,rb,(int*)rc,(int*)ro,rd,c)

#define MPI_Allgatherv(sb,sc,sd,rb,rc,ro,rd,c) \
  MPI_Allgatherv(sb,sc,sd,rb,(int*)rc,(int*)ro,rd,c)




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static const size_t NRUNS = 100;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX dim
#define DLMEM_TYPE_T dim_t
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX real
#define DLMEM_TYPE_T real_t
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMATH_PREFIX dim
#define DLMATH_TYPE_T dim_t
#define DLMATH_DLTYPE DLTYPE_INTEGRAL
#define DLMATH_STATIC
#include "dlmath_headers.h"
#undef DLMATH_STATIC
#undef DLMATH_DLTYPE
#undef DLMATH_TYPE_T
#undef DLMATH_PREFIX





/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __print_spm(
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const ptr,
    const dim_t * const ind,
    const real_t * const val)
{
  dim_t i, k;
  ind_t j;

  for (i=0;i<nrows;++i) {
    printf("[");
    j = ptr[i];
    for (k=0;k<ncols;++k) {
      if (ind[j] == k && j < ptr[i+1] && ptr[i] != ptr[i+1]) {
        printf(" "PF_REAL_T" ",val[j]);
        ++j;
      } else {
        printf(" "PF_REAL_T" ",0.0);
      }
    }
    printf("]\n");
  }
}



static dim_t * __generate_pre(
    const dim_t * const src, 
    const dim_t n)
{
  dim_t * dst = dim_alloc(n+1);
  dim_copy(dst,src,n); 
  dst[n] = 0;
  dim_prefixsum_exc(dst,n+1);
  return dst;
}



/**
 * This command is invoked with two arguments:
 *    mpi_spmultsp spmatA.csr spmatB.csr out.mat
 *
 * Each process then reads from spmatA.csr.P and spmatB.csr.P to read in its
 * portions of the matrices, and writes out.mat.P
 *
 */
int main(int argc, char ** argv)
{
  int rv, err;
  int rank, size; 
  size_t r;
  char aname[512], cname[512], * bname;
  dl_timer_t iotmr, cputmr, nettmr;
  dim_t anrows, ancols, bnrows, bncols, cnrows, cncols;
  ind_t * arowptr = NULL;
  dim_t * arowind = NULL;
  real_t * arowval = NULL, * browval = NULL, * crowval = NULL;

  if (sizeof(dim_t) != sizeof(int) || sizeof(ind_t) != sizeof(int)) {
    eprintf("Wrong type size. Aborting.\n");
    abort();
  }

  rv = 0;

  if (MPI_Init(&argc,&argv) != MPI_SUCCESS) {
    eprintf("Initialization failed, aborting...\n");
    abort();
  }

  if (argc != 4) {
    eprintf("Invalid number of arguments\n");
    eprintf("USAGE: %s a.csr b.vec c.vec\n",argv[0]);
    goto END;
  }

  dl_init_timer(&iotmr);
  dl_init_timer(&nettmr);
  dl_init_timer(&cputmr);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size); 

  /* setup file names */
  sprintf(aname,"%s.%d",argv[1],rank);
  bname = argv[2];
  sprintf(cname,"%s.%d",argv[3],rank);

  dl_start_timer(&iotmr);

  /* read in matrix chunks */ 
  err = gooseberry_read_sparse_matrix(GOOSEBERRY_FORMAT_CSR,aname,&anrows,
      &ancols,&arowptr,&arowind,&arowval);
  if (err != GOOSEBERRY_SUCCESS) {
    eprintf("Error reading first matrix: '%s'\n",aname);
    rv = 1;
    goto END;
  }

  err = gooseberry_read_dense_matrix(GOOSEBERRY_FORMAT_GRID,bname,&bnrows,
      &bncols,&browval);
  if (err != GOOSEBERRY_SUCCESS) {
    eprintf("Error reading first matrix: '%s'\n",bname);
    rv = 1;
    goto END;
  }

  /* just in case a's zero rows were ignored */
  if (ancols < bnrows) {
    ancols = bnrows;
  } else if (ancols > bnrows) {
    eprintf("B matrix has too few rows: found "PF_DIM_T" and expected atleast "
        PF_DIM_T"\n",bnrows,ancols);
    rv = 1;
    goto END;
  }

  cnrows = anrows;
  cncols = 1;

  dl_stop_timer(&iotmr);

  MPI_Barrier(MPI_COMM_WORLD);

  dl_start_timer(&cputmr);

  crowval = real_alloc(cnrows*cncols);

  for (r=0;r<NRUNS;++r) {
    gooseberry_spmult(anrows,ancols,bncols,arowptr,arowind,arowval,browval,
        crowval,NULL,0);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  dl_stop_timer(&cputmr);

  dl_start_timer(&iotmr);

  err = gooseberry_write_dense_matrix(GOOSEBERRY_FORMAT_GRID,cname,cnrows,
      cncols,crowval);

  MPI_Barrier(MPI_COMM_WORLD);
  dl_stop_timer(&iotmr);

  if (rank == 0) {
    dl_print_header("TIMING",'=');
    printf("IO: %lf\n",dl_poll_timer(&iotmr));
    printf("COMPUTATION: %lf\n",dl_poll_timer(&cputmr));
    dl_print_footer('=');
  }

  END:

  if (arowptr) {
    dl_free(arowptr);
  }
  if (arowind) {
    dl_free(arowind);
  }
  if (arowval) {
    dl_free(arowval);
  }
  if (browval) {
    dl_free(browval);
  }
  if (crowval) {
    dl_free(crowval);
  }

  MPI_Finalize();

  return rv;
}


#endif

