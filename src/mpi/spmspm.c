/**
 * @file mpi_spmultsp.c
 * @brief Perform a sparse matrix matrix multiplication via MPI
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-05-13
 */




#ifndef MPI_SPMULTSP_C
#define MPI_SPMULTSP_C




#include <mpi.h>
#include <domlib.h>
#include <gooseberry.h>





/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/





#define DLMEM_PREFIX real
#define DLMEM_TYPE_T real_t
#define DLMEM_STATIC
#define DLMEM_DLTYPE DLTYPE_FLOAT
#include "dlmem_headers.h"
#undef DLMEM_DLTYPE
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX dim
#define DLMEM_TYPE_T dim_t
#define DLMEM_STATIC
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_headers.h"
#undef DLMEM_DLTYPE
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX 


#define DLMEM_PREFIX ind 
#define DLMEM_TYPE_T ind_t
#define DLMEM_STATIC
#define DLMEM_DLTYPE DLTYPE_INTEGRAL
#include "dlmem_headers.h"
#undef DLMEM_DLTYPE
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX 


#define DLMATH_PREFIX dim
#define DLMATH_TYPE_T dim_t
#define DLMATH_STATIC
#include "dlmath_headers.h"
#undef DLMATH_STATIC
#undef DLMATH_TYPE_T
#undef DLMATH_PREFIX


#define DLMATH_PREFIX ind
#define DLMATH_TYPE_T ind_t
#define DLMATH_STATIC
#include "dlmath_headers.h"
#undef DLMATH_STATIC
#undef DLMATH_TYPE_T
#undef DLMATH_PREFIX




/******************************************************************************
* MPI MACROS ******************************************************************
******************************************************************************/


#define MPI_DIM_T MPI_UNSIGNED
#define MPI_REAL_T MPI_FLOAT


#define MPI_Alltoallv(sb,sc,so,sd,rb,rc,ro,rd,c) \
  MPI_Alltoallv(sb,(int*)sc,(int*)so,sd,rb,(int*)rc,(int*)ro,rd,c)

#define MPI_Allgatherv(sb,sc,sd,rb,rc,ro,rd,c) \
  MPI_Allgatherv(sb,sc,sd,rb,(int*)rc,(int*)ro,rd,c)





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




/******************************************************************************
* MAIN ************************************************************************
******************************************************************************/


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
  int rv, err, r;
  int rank, size; 
  int p;
  char aname[512], bname[512], cname[512];
  dl_timer_t iotmr, cputmr, nettmr, mpi_tmr;
  dim_t anrows, ancols, bnrows, mbncols,bncols, cnrows, cncols, i, k, nnr;
  ind_t j, l;
  ind_t * arowptr = NULL, * browptr = NULL, * crowptr = NULL, * fbrowptr = NULL;
  dim_t * arowind = NULL, * browind = NULL, * fbrowind = NULL, *crowind = NULL,
        * obrowind = NULL;
  real_t * arowval = NULL, * browval = NULL, * fbrowval = NULL, 
         * crowval = NULL;
  dim_t * mark, * rbc; 
  /* my buffers */
  ind_t nnz;
  dim_t * nr, * nrptr, * nrptr_ps, *ni, * xnnr, *xnnr_ps, *cnt, * ptr;
  real_t * nz;
  /* remote buffers */
  dim_t req_nnz;
  dim_t * req_nnr, * req_nnr_ps, * req_nr, * req_nrptr, * req_nrptr_ps,
        * req_ni, *req_cnt, *req_ptr;
  real_t * req_nz;

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
    eprintf("USAGE: %s a.csr b.csr c.mat\n",argv[0]);
    goto END;
  }

  dl_init_timer(&iotmr);
  dl_init_timer(&nettmr);
  dl_init_timer(&cputmr);
  dl_init_timer(&mpi_tmr);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size); 

  /* setup file names */
  sprintf(aname,"%s.%d",argv[1],rank);
  sprintf(bname,"%s.%d",argv[2],rank);
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

  err = gooseberry_read_sparse_matrix(GOOSEBERRY_FORMAT_CSR,bname,&bnrows,
      &bncols,&browptr,&obrowind,&browval);
  if (err != GOOSEBERRY_SUCCESS) {
    eprintf("Error reading first matrix: '%s'\n",bname);
    rv = 1;
    goto END;
  }
  browind = dim_alloc(browptr[bnrows]);

  /* just in case a's zero rows were ignored */
  if (ancols < bnrows) {
    ancols = bnrows;
  } else if (ancols > bnrows) {
    eprintf("B matrix has too few rows: found "PF_DIM_T" and expected atleast "
        PF_DIM_T"\n",bnrows,ancols);
    rv = 1;
    goto END;
  }

  dl_stop_timer(&iotmr);

  for (r=0;r<10;++r) {
    /* communicate needed portions of B ****************************************/
    MPI_Barrier(MPI_COMM_WORLD);
    dl_start_timer(&nettmr);

    /* build an ordered list of the rows needed of b */
    mark = dim_calloc(ancols);
    nnr = 0;
    for (i=0;i<anrows;++i) {
      for (j=arowptr[i];j<arowptr[i+1];++j) {
        k = arowind[j];
        if (!mark[k]) {
          mark[k] = 1;
          ++nnr;
        }
      }
    }
    nr = dim_alloc(nnr);
    nnr = 0;
    for (i=0;i<ancols;++i) {
      if (mark[i]) {
        nr[nnr++] = i;
      }
    }
    dl_free(mark);

    /* determine size of C and modify my B */
    cnrows = anrows;
    rbc = dim_alloc(size+1);
    rbc[size] = 0;
    dl_start_timer(&mpi_tmr);
    MPI_Allgather(&bncols,1,MPI_DIM_T,rbc,1,MPI_DIM_T,MPI_COMM_WORLD);
    dl_stop_timer(&mpi_tmr);
    dim_prefixsum_exc(rbc,size+1);
    for (i=0;i<browptr[bnrows];++i) {
      browind[i] = obrowind[i] + rbc[rank];
    }
    mbncols = cncols = rbc[size];
    dl_free(rbc);


    //printf("Output will be sparse "PF_DIM_T"x"PF_DIM_T" matrix\n",cnrows,cncols);

    #ifdef DEBUG
    for (p=0;p<size;++p) {
      if (rank == p) {
        printf("%d Multiplying "PF_DIM_T"x"PF_DIM_T" matrix by "PF_DIM_T"x"
            PF_DIM_T" matrix to get a "PF_DIM_T"x"PF_DIM_T" matrix.\n",rank,
            anrows,ancols,bnrows,mbncols,cnrows,cncols);
        printf("%d A\n",rank);
        __print_spm(anrows,ancols,arowptr,arowind,arowval);
        printf("%d B\n",rank);
        __print_spm(bnrows,mbncols,browptr,browind,browval);
        fflush(stdout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    #endif


    /* let everyone know how many rows I need -- and find out how many they 
     * need */
    req_nnr = dim_alloc(size);
    dl_start_timer(&mpi_tmr);
    MPI_Allgather(&nnr,1,MPI_DIM_T,req_nnr,1,MPI_DIM_T,MPI_COMM_WORLD);
    dl_stop_timer(&mpi_tmr);

    /* do a prefixsum on req_nnr so I can figure where to store requested rows */
    req_nnr_ps = __generate_pre(req_nnr,size);

    /* now actually send/recv row numbers */
    req_nr = dim_alloc(req_nnr_ps[size]);
    dl_start_timer(&mpi_tmr);
    MPI_Allgatherv(nr,nnr,MPI_DIM_T,req_nr,req_nnr,req_nnr_ps,MPI_DIM_T,
        MPI_COMM_WORLD);
    dl_stop_timer(&mpi_tmr);

    /* build a rowptr to send -- unprefix summed */
    req_nrptr = dim_alloc(req_nnr_ps[size]);
    for (p=0;p<size;++p) {
      for (i=req_nnr_ps[p];i<req_nnr_ps[p+1];++i) {
        req_nrptr[i] = browptr[req_nr[i]+1] - browptr[req_nr[i]];
      }
    }

    /* communicate row sizes */
    xnnr = dim_init_alloc(nnr,size);
    xnnr_ps = __generate_pre(xnnr,size);
    nrptr = dim_alloc(xnnr_ps[size]+1);
    dl_start_timer(&mpi_tmr);
    MPI_Alltoallv(req_nrptr,req_nnr,req_nnr_ps,MPI_DIM_T,nrptr,xnnr,xnnr_ps,
        MPI_DIM_T,MPI_COMM_WORLD);
    dl_stop_timer(&mpi_tmr);
    dl_free(xnnr);

    /* setup structures to recieve elements */
    nrptr_ps = __generate_pre(nrptr,xnnr_ps[size]);
    nnz = nrptr_ps[xnnr_ps[size]];

    cnt = dim_alloc(size);
    ptr = dim_alloc(size+1);
    ni = dim_alloc(nnz);
    nz = real_alloc(nnz);
    ptr[0] = 0;
    for (p=0;p<size;++p) {
      cnt[p] = 0;
      for (i=xnnr_ps[p];i<xnnr_ps[p+1];++i) {
        cnt[p] += nrptr[i];
      }
      ptr[p+1] = ptr[p] + cnt[p];
    }
    dl_free(xnnr_ps);
    dl_free(nrptr);

    /* setup and fill structure to send elements */
    req_nrptr_ps = __generate_pre(req_nrptr,req_nnr_ps[size]);
    dl_free(req_nrptr);
    req_nnz = req_nrptr_ps[req_nnr_ps[size]];


    req_ni = dim_alloc(req_nnz);
    req_nz = real_alloc(req_nnz);
    req_cnt = dim_alloc(size);
    req_ptr = dim_alloc(size+1);
    req_ptr[0] = 0;
    l = 0;
    for (p=0;p<size;++p) {
      for (i=req_nnr_ps[p];i<req_nnr_ps[p+1];++i) {
        for (j=browptr[req_nr[i]];j<browptr[req_nr[i]+1];++j) {
          req_ni[l] = browind[j];
          req_nz[l] = browval[j];
          ++l;
        }
      }
      req_cnt[p] = l - req_ptr[p];
      req_ptr[p+1] = l;
    }

    /* send elements and indices */
    dl_start_timer(&mpi_tmr);
    MPI_Alltoallv(req_ni,req_cnt,req_ptr,MPI_DIM_T,ni,cnt,ptr,MPI_DIM_T,
        MPI_COMM_WORLD);
    MPI_Alltoallv(req_nz,req_cnt,req_ptr,MPI_REAL_T,nz,cnt,ptr,MPI_REAL_T,
        MPI_COMM_WORLD);
    dl_stop_timer(&mpi_tmr);

    dl_free(req_ni);
    dl_free(req_nz);
    dl_free(req_cnt);
    dl_free(req_ptr);
    dl_free(req_nrptr_ps);

    /* ugly trakcing stuff */
    dim_t * rnnz = NULL, * snnz = NULL;
    dim_t new_nnz = nnz;
    dim_t max_snnz, max_rnnz, trnnz, tsnnz;
    if (r == 0) {
    if (rank == 0) {
      rnnz = dim_alloc(size);
      snnz = dim_alloc(size);
    }

    MPI_Gather(&new_nnz,1,MPI_DIM_T,rnnz,1,MPI_DIM_T,0,MPI_COMM_WORLD);
    MPI_Gather(&req_nnz,1,MPI_DIM_T,snnz,1,MPI_DIM_T,0,MPI_COMM_WORLD);

    trnnz = 0;
    tsnnz = 0;
    max_snnz = 0;
    max_rnnz = 0;
    if (rank == 0) {
      for (p=0;p<size;++p) {
        //printf("%d: recieving "PF_IND_T" non-zeros\n",p,rnnz[p]);
        //printf("%d: sending "PF_IND_T" non-zeros\n",p,snnz[p]);
        trnnz += rnnz[p];
        tsnnz += snnz[p];
        if (rnnz[p] > max_rnnz) {
          max_rnnz = rnnz[p];
        }
        if (snnz[p] > max_snnz) {
          max_snnz = snnz[p];
        }
      }
      printf("Recieve imbalance: %f\n",(max_rnnz*size)/(float)trnnz); 
      printf("Send imbalance: %f\n",(max_snnz*size)/(float)tsnnz); 
      dl_free(rnnz);
      dl_free(snnz);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    }
      
    /* assemble my b matrix */
    fbrowptr = ind_calloc(bnrows+1);
    fbrowind = dim_alloc(nnz);
    fbrowval = real_alloc(nnz);

    for (p=0;p<size;++p) {
      k = nnr*p;
      for (i=0;i<nnr;++i) {
        for (j=nrptr_ps[k+i];j<nrptr_ps[k+i+1];++j) {
          ++fbrowptr[nr[i]+1];
        }
      }
    }
    ind_prefixsum_exc(fbrowptr+1,bnrows);
    for (p=0;p<size;++p) {
      k = nnr*p;
      for (i=0;i<nnr;++i) {
        for (j=nrptr_ps[k+i];j<nrptr_ps[k+i+1];++j) {
          fbrowind[fbrowptr[nr[i]+1]] = ni[j];
          fbrowval[fbrowptr[nr[i]+1]] = nz[j];
          ++fbrowptr[nr[i]+1];
        }
      }
    }


    dl_free(ni);
    dl_free(nz);
    dl_free(nr);
    dl_free(nrptr_ps);


    dl_start_timer(&mpi_tmr);
    MPI_Barrier(MPI_COMM_WORLD);
    dl_stop_timer(&mpi_tmr);
    dl_stop_timer(&nettmr);



    #ifdef DEBUG
    for (p=0;p<size;++p) {
      if (rank == p) {
        printf("%d B'\n",rank);
        __print_spm(bnrows,mbncols,fbrowptr,fbrowind,fbrowval);
        fflush(stdout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    #endif


    /* compute section of C ****************************************************/

    dl_start_timer(&cputmr);


    gooseberry_spmultsp_sp(anrows,ancols,mbncols,arowptr,arowind,arowval,
        fbrowptr,fbrowind,fbrowval,&crowptr,&crowind,&crowval,NULL,0);

    dl_free(fbrowptr);
    dl_free(fbrowind);
    dl_free(fbrowval);

    MPI_Barrier(MPI_COMM_WORLD);
    dl_stop_timer(&cputmr);

    dl_start_timer(&iotmr);

    err = gooseberry_write_sparse_matrix(GOOSEBERRY_FORMAT_CSR,cname,cnrows,
        cncols,crowptr,crowind,crowval);

    MPI_Barrier(MPI_COMM_WORLD);
    dl_stop_timer(&iotmr);

    dl_free(crowptr);
    dl_free(crowind);
    dl_free(crowval);
  }

  if (rank == 0) {
    dl_print_header("TIMING",'=');
    printf("IO: %lf\n",dl_poll_timer(&iotmr));
    printf("COMMUNICATION: %lf\n",dl_poll_timer(&nettmr));
    printf("\tMPI_CALLS: %lf\n",dl_poll_timer(&mpi_tmr));
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
  if (browptr) {
    dl_free(browptr);
  }
  if (browind) {
    dl_free(browind);
  }
  if (browval) {
    dl_free(browval);
  }

  MPI_Finalize();

  return rv;
}


#endif

