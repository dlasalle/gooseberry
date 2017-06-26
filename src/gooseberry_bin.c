/**
 * @file bowstring_bin.c
 * @brief Command line interface for gooseberry
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-04-28
 */




#ifndef GOOSEBERRY_C
#define GOOSEBERRY_C




#include "base.h"
#include "matrix.h"
#include "io.h"
#include "blas.h"
#include "cgd.h"
#include "permute.h"
#include "analyze.h"




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define ARRAY_SIZE(a) \
  (sizeof(a) > 0 ?  (sizeof(a) / sizeof((a)[0])) : 0)

#ifndef NO_OMP
#define DEFAULT_NUMTHREADS omp_get_max_threads()
#else
#define DEFAULT_NUMTHREADS 1
#endif




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


/* COMMANDS ******************************************************************/

typedef enum command_t {
  COMMAND_HELP,
  COMMAND_ANALYSIS,
  COMMAND_PERMUTE,
  COMMAND_TRANSFORM,
  COMMAND_GENERATE,
  COMMAND_BLAS,
  COMMAND_CGD,
  COMMAND_SGD,
  COMMAND_PAGERANK
} command_t;



/* ANALYSIS ******************************************************************/

typedef enum analysis_t {
  ANALYSIS_MATRIXSTATS,
  ANALYSIS_CHOLESKY
} analysis_t;


typedef enum analysis_option_t {
  ANALYSIS_OPTION_HELP,
  ANALYSIS_OPTION_INFILE,
  ANALYSIS_OPTION_TIME,
  ANALYSIS_OPTION_TYPE,
  ANALYSIS_OPTION_PERMFILE 
} analysis_option_t;



/* PERMUTE *******************************************************************/

typedef enum permute_option_t {
  PERMUTE_OPTION_HELP,
  PERMUTE_OPTION_INFILE,
  PERMUTE_OPTION_OUTFILE,
  PERMUTE_OPTION_PERMUTATION,
  PERMUTE_OPTION_TIME,
  PERMUTE_OPTION_ROWPERM,
  PERMUTE_OPTION_COLPERM
} permute_option_t;


typedef enum permute_permutation_t {
  PERMUTE_PERMUTATION_FILE,
  PERMUTE_PERMUTATION_RANDOM,
  PERMUTE_PERMUTATION_ROWRANDOM,
  PERMUTE_PERMUTATION_COLRANDOM,
  PERMUTE_PERMUTATION_BANDWIDTH
} permute_permutation_t;


/* TRANSFORM *****************************************************************/

typedef enum transform_option_t {
  TRANSFORM_OPTION_HELP,
  TRANSFORM_OPTION_INFILE,
  TRANSFORM_OPTION_OUTFILE,
  TRANSFORM_OPTION_PARTFILE,
  TRANSFORM_OPTION_TIME,
  TRANSFORM_OPTION_OPERATION
} transform_option_t;


typedef enum transform_operation_t {
  TRANSFORM_OPERATION_CONVERT,
  TRANSFORM_OPERATION_SYMMETRIFY,
  TRANSFORM_OPERATION_DEBIPARTIFY,
  TRANSFORM_OPERATION_ROWSPLIT,
  TRANSFORM_OPERATION_COLSPLIT,
  TRANSFORM_OPERATION_ROWJOIN,
  TRANSFORM_OPERATION_COLJOIN,
  TRANSFORM_OPERATION_TRANSPOSE
} transform_operation_t;


/* GENERATE ******************************************************************/

typedef enum generate_option_t {
  GENERATE_OPTION_HELP,
  GENERATE_OPTION_OUTFILE,
  GENERATE_OPTION_TYPE,
  GENERATE_OPTION_SIZE,
  GENERATE_OPTION_TIME
} generate_option_t;


typedef enum generate_type_t {
  GENERATE_TYPE_NULL,
  GENERATE_TYPE_DENSE_VECTOR
} generate_type_t;


/* BLAS **********************************************************************/

typedef enum blas_option_t {
  BLAS_OPTION_HELP,
  BLAS_OPTION_OPERATION,
  BLAS_OPTION_INFILE,
  BLAS_OPTION_OUTFILE,
  BLAS_OPTION_TIME,
  BLAS_OPTION_RUNS,
  BLAS_OPTION_THREADS,
  BLAS_OPTION_ROWPERM,
  BLAS_OPTION_COLPERM,
  BLAS_OPTION_REDUCEBANDWIDTH
} blas_option_t;


typedef enum blas_operation_t {
  BLAS_OPERATION_NOOP,
  BLAS_OPERATION_MULTIPLY
} blas_operation_t;


/* CGD ***********************************************************************/

typedef enum cgd_option_t {
  CGD_OPTION_HELP,
  CGD_OPTION_INFILE,
  CGD_OPTION_OUTFILE,
  CGD_OPTION_ERROR,
  CGD_OPTION_NITER,
  CGD_OPTION_TIME,
  CGD_OPTION_RUNS,
  CGD_OPTION_THREADS,
  CGD_OPTION_ROWPERM,
  CGD_OPTION_COLPERM
} cgd_option_t;


/* PAGERANK ******************************************************************/

typedef enum pagerank_option_t {
  PAGERANK_OPTION_HELP,
  PAGERANK_OPTION_INFILE,
  PAGERANK_OPTION_OUTFILE,
  PAGERANK_OPTION_ERROR,
  PAGERANK_OPTION_DAMPING,
  PAGERANK_OPTION_NITER,
  PAGERANK_OPTION_TIME,
  PAGERANK_OPTION_RUNS,
  PAGERANK_OPTION_THREADS,
  PAGERANK_OPTION_PERM
} pagerank_option_t;





/******************************************************************************
* OPTION ARRAYS ***************************************************************
******************************************************************************/


/* COMMANDS ******************************************************************/

static const cmd_opt_pair_t COMMANDS[] = {
  [COMMAND_HELP] = {"help","Display list of available commands.", \
      COMMAND_HELP},
  [COMMAND_ANALYSIS] = {"analysis","Perform an analysis on a matrix/vector.", \
      COMMAND_ANALYSIS},
  [COMMAND_TRANSFORM] = {"transform","Transform a matrix/vector.", \
      COMMAND_TRANSFORM},
  [COMMAND_PERMUTE] = {"permute","Permute a matrix/vector.",COMMAND_PERMUTE},
  [COMMAND_GENERATE] = {"generate","Generate a matrix/vector.", \
      COMMAND_GENERATE},
  [COMMAND_BLAS] = {"blas","Perform a blas operation.",COMMAND_BLAS},
  [COMMAND_CGD] = {"cgd","Perform conjugate gradient descent.",COMMAND_CGD},
  [COMMAND_SGD] = {"sgd","Perform stocastic gradient descent.",COMMAND_SGD},
  [COMMAND_PAGERANK] = {"pagerank","Perform a pagerank on a square matrix.", \
      COMMAND_PAGERANK}
};
static const size_t NCOMMANDS = ARRAY_SIZE(COMMANDS);


/* ANALYSIS ******************************************************************/

static const cmd_opt_pair_t ANALYSIS[] = {
  {"matrixstats","Calculate statistics of a matrix/vector", \
      ANALYSIS_MATRIXSTATS},
  {"cholesky","Calculate the stats associated with a cholesky decomposition", \
      ANALYSIS_CHOLESKY}
};


static const cmd_opt_t ANALYSIS_OPTS[] = {
  {ANALYSIS_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG, \
      NULL,0},
  {ANALYSIS_OPTION_INFILE,'i',"infile","The input matrix.",CMD_OPT_STRING, \
      NULL,0},
  {ANALYSIS_OPTION_PERMFILE,'p',"permfile","The permutation vector.", \
      CMD_OPT_STRING,NULL,0},
  {ANALYSIS_OPTION_TYPE,'a',"type","The type of analysis to perform.", \
      CMD_OPT_CHOICE,ANALYSIS,ARRAY_SIZE(ANALYSIS)},
  {PERMUTE_OPTION_TIME,'t',"times","Print timing of the analysis.", \
      CMD_OPT_FLAG,NULL,0}
};
static const size_t NANALYSIS_OPTS = ARRAY_SIZE(ANALYSIS_OPTS);


/* PERMUTE *******************************************************************/

static const cmd_opt_pair_t PERMUTE_PERMUTATIONS[] = {
  {"random","Perform a random permutation on the rows and columns.", \
      PERMUTE_PERMUTATION_RANDOM},
  {"file","Perform a permuation based on input files (specified with -R " \
    "and/or -C, can be permutations or partitions).", \
      PERMUTE_PERMUTATION_FILE},
  {"bandwidth","Perform a bandwidth reducing permutation.", \
      PERMUTE_PERMUTATION_BANDWIDTH},
};

static const cmd_opt_t PERMUTE_OPTS[] = {
  {PERMUTE_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG, \
      NULL,0},
  {PERMUTE_OPTION_INFILE,'i',"infile","An input matrix/vector file.", \
      CMD_OPT_STRING,NULL,0},
  {PERMUTE_OPTION_OUTFILE,'o',"outfile","The output vector file.", \
      CMD_OPT_STRING,NULL,0},
  {PERMUTE_OPTION_PERMUTATION,'p',"permutation","The type of permutation to " \
      "perform.",CMD_OPT_CHOICE,PERMUTE_PERMUTATIONS, \
      ARRAY_SIZE(PERMUTE_PERMUTATIONS)},
  {PERMUTE_OPTION_TIME,'t',"times","Print timing of the permutation.", \
      CMD_OPT_FLAG,NULL,0},
  {PERMUTE_OPTION_ROWPERM,'R',"rowperm","Row permutation/partition file.", \
    CMD_OPT_STRING,NULL,0},
  {PERMUTE_OPTION_COLPERM,'C',"colperm","Column permutation/partition file.", \
      CMD_OPT_STRING,NULL,0}
};
static const size_t NPERMUTE_OPTS = ARRAY_SIZE(PERMUTE_OPTS);


/* TRANSFROM *****************************************************************/

static const cmd_opt_pair_t TRANSFORM_OPERATIONS[] = {
  {"convert","Convert from one matrix/vector format to another.", \
      TRANSFORM_OPERATION_CONVERT},
  {"symmetrify","Transform to a symmetric matrix: B = A + A^T.", \
      TRANSFORM_OPERATION_SYMMETRIFY},
  {"debipartify","Transform to a symmetric matrix: B = [ 0 , A ; A^T , 0].", \
      TRANSFORM_OPERATION_DEBIPARTIFY},
  {"rowsplit","Split a matrix row-wise into submatrices: " \
    "[ B.0 ; B.1 ; ... ] = A.",TRANSFORM_OPERATION_ROWSPLIT},
  {"colsplit","Split a matrix column-wise into submatrices: " \
    "[ B.0 , B.1 , ... ] = A.",TRANSFORM_OPERATION_COLSPLIT},
  #ifdef XXX
  {"rowjoin","Join submatrices row-wise into a single matrix: " \
    "B = [ A.0 ; A.1 ; ... ].",TRANSFORM_OPERATION_ROWJOIN},
  {"coljoin","Join submatrices column-wise into a single matrix: " \
    "B = [ A.0 , A.1 , ... ].",TRANSFORM_OPERATION_COLJOIN},
  #endif
  {"transpose","Transpose a matrix: " \
    "B = A^T.",TRANSFORM_OPERATION_TRANSPOSE}
};


static const cmd_opt_t TRANSFORM_OPTS[] = {
  {TRANSFORM_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG, \
      NULL,0},
  {TRANSFORM_OPTION_INFILE,'i',"infile","An input matrix/vector file.", \
      CMD_OPT_STRING,NULL,0},
  {TRANSFORM_OPTION_OUTFILE,'o',"outfile","The output vector file.", \
      CMD_OPT_STRING,NULL,0},
  {TRANSFORM_OPTION_PARTFILE,'p',"partfile","The partition vector file.", \
      CMD_OPT_STRING,NULL,0},
  {TRANSFORM_OPTION_OPERATION,'x',"operation","The type of permutation to " \
      "perform.",CMD_OPT_CHOICE,TRANSFORM_OPERATIONS, \
      ARRAY_SIZE(TRANSFORM_OPERATIONS)},
  {TRANSFORM_OPTION_TIME,'t',"times","Print timing of the permutation.", \
      CMD_OPT_FLAG,NULL,0}
};
static const size_t NTRANSFORM_OPTS = ARRAY_SIZE(TRANSFORM_OPTS);


/* TRANSFROM *****************************************************************/

static const cmd_opt_pair_t GENERATE_OPERATIONS[] = {
  {"vector","Generate a dense vector.",GENERATE_TYPE_DENSE_VECTOR}
};


static const cmd_opt_t GENERATE_OPTS[] = {
  {GENERATE_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG, \
      NULL,0},
  {GENERATE_OPTION_OUTFILE,'o',"outfile","The output vector file.", \
      CMD_OPT_STRING,NULL,0},
  {GENERATE_OPTION_TYPE,'g',"type","The type of permutation to " \
      "perform.",CMD_OPT_CHOICE,GENERATE_OPERATIONS, \
      ARRAY_SIZE(GENERATE_OPERATIONS)},
  {GENERATE_OPTION_SIZE,'s',"size","Size of the generate matrix/vector", \
      CMD_OPT_INT,NULL,0},
  {GENERATE_OPTION_TIME,'t',"times","Print timing of the permutation.", \
      CMD_OPT_FLAG,NULL,0}
};
static const size_t NGENERATE_OPTS = ARRAY_SIZE(GENERATE_OPTS);



/* BLAS **********************************************************************/

static const cmd_opt_pair_t BLAS_OPERATIONS[] = {
  {"noop","Perform no operation, just copy the input matrix/vector.", \
      BLAS_OPERATION_NOOP},
  {"multiply","Multiply a matrix/vector with a matrix/vector.", \
      BLAS_OPERATION_MULTIPLY}
};


static const cmd_opt_t BLAS_OPTS[] = {
  {BLAS_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG, \
    NULL,0},
  {BLAS_OPTION_OPERATION,'x',"operation","The type of operation to perform.", \
      CMD_OPT_CHOICE,BLAS_OPERATIONS,ARRAY_SIZE(BLAS_OPERATIONS)},
  {BLAS_OPTION_INFILE,'i',"infile","An input matrix/vector file.", \
      CMD_OPT_STRING,NULL,0},
  {BLAS_OPTION_OUTFILE,'o',"outfile","The output matrix/vector file.", \
      CMD_OPT_STRING,NULL,0},
  {BLAS_OPTION_TIME,'t',"times","Print timing of the blas routines.", \
      CMD_OPT_FLAG,NULL,0},
  {BLAS_OPTION_RUNS,'r',"runs","Number of repeated runs (only useful for " \
    "timing purposes).",CMD_OPT_INT,NULL,0},
#ifndef NO_OMP
  {BLAS_OPTION_THREADS,'T',"threads","Number of threads.",CMD_OPT_INT,NULL,0}, 
#endif
  {BLAS_OPTION_ROWPERM,'R',"rowperm","Row permutation file.",CMD_OPT_STRING, \
      NULL,0},
  {BLAS_OPTION_COLPERM,'C',"colperm","Column permutation file.", \
      CMD_OPT_STRING,NULL,0},
  {BLAS_OPTION_REDUCEBANDWIDTH,'b',"bandwidthreduce","Re-order the matrix " \
      "to reduce bandwidth.",CMD_OPT_FLAG, NULL,0}
};
static const size_t NBLAS_OPTS = ARRAY_SIZE(BLAS_OPTS);


/* CGD ****************************************************************/

static const cmd_opt_t CGD_OPTS[] = {
  {CGD_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG, \
      NULL,0},
  {CGD_OPTION_INFILE,'i',"infile","An input matrix/vector file.", \
      CMD_OPT_STRING,NULL,0},
  {CGD_OPTION_OUTFILE,'o',"outfile","The output vector file.", \
      CMD_OPT_STRING,NULL,0},
  {CGD_OPTION_ERROR,'e',"error","The RMSE to achieve before exiting.", \
      CMD_OPT_FLOAT,NULL,0},
  {CGD_OPTION_NITER,'I',"iter","The number of iterations to run before " \
      "exiting.",CMD_OPT_INT,NULL,0},
  {CGD_OPTION_TIME,'t',"times","Print timing of the cgd routines.", \
      CMD_OPT_FLAG,NULL,0},
  {CGD_OPTION_RUNS,'r',"runs","Number of repeated runs (only useful for " \
    "timing purposes).",CMD_OPT_INT,NULL,0},
#ifndef NO_OMP
  {CGD_OPTION_THREADS,'T',"threads","Number of threads.",CMD_OPT_INT,NULL,0},
#endif
  {CGD_OPTION_ROWPERM,'R',"rowperm","Row permutation file.",CMD_OPT_STRING, \
      NULL,0},
  {CGD_OPTION_COLPERM,'C',"colperm","Column permutation file.", \
      CMD_OPT_STRING,NULL,0}
};
static const size_t NCGD_OPTS = ARRAY_SIZE(CGD_OPTS);


/* PAGERANK ******************************************************************/

static const cmd_opt_t PAGERANK_OPTS[] = {
  {PAGERANK_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG, \
      NULL,0},
  {PAGERANK_OPTION_INFILE,'i',"infile","An input matrix file.", \
      CMD_OPT_STRING,NULL,0},
  {PAGERANK_OPTION_OUTFILE,'o',"outfile","The output matrix/vector file.", \
      CMD_OPT_STRING,NULL,0},
  {PAGERANK_OPTION_TIME,'t',"times","Print timing of the pagerank " \
      "calcuation.", CMD_OPT_FLAG,NULL,0},
  {PAGERANK_OPTION_RUNS,'r',"runs","Number of repeated runs (only useful " \
      "for timing purposes).",CMD_OPT_INT,NULL,0},
  {PAGERANK_OPTION_NITER,'I',"iterations","Maximum number of iterations.", \
      CMD_OPT_INT,NULL,0},
  {PAGERANK_OPTION_ERROR,'e',"error","Error threshold for stopping.", \
      CMD_OPT_FLOAT,NULL,0},
  {PAGERANK_OPTION_DAMPING,'d',"damping","Damping factor to use.", \
      CMD_OPT_FLOAT,NULL,0},
#ifndef NO_OMP
  {PAGERANK_OPTION_THREADS,'T',"threads","Number of threads.",CMD_OPT_INT, \
    NULL,0},
#endif
  {PAGERANK_OPTION_PERM,'p',"perm","Row and column permutation file.", \
    CMD_OPT_STRING,NULL,0}
};
static const size_t NPAGERANK_OPTS = ARRAY_SIZE(PAGERANK_OPTS);





/* FILE TYPES ****************************************************************/

static const char * RAW_EXTENSIONS[] = {"raw",NULL};
static const char * GRID_EXTENSIONS[] = {"mat","grid","vec","txt",NULL};
static const char * CSR_EXTENSIONS[] = {"csr",NULL};
static const char * SVM_EXTENSIONS[] = {"svm","libfm",NULL};
static const char * POINT_EXTENSIONS[] = {"ij","point",NULL};
static const char * GRAPH_EXTENSIONS[] = {"metis","chaco","graph",NULL};
static const char * CLU_EXTENSIONS[] = {"clu",NULL};
static const char * MATRIXMARKET_EXTENSIONS[] = {"mm","mtx",NULL};

static const char * const * const FILE_TYPES[] = {
  [GOOSEBERRY_FORMAT_RAW] = RAW_EXTENSIONS,
  [GOOSEBERRY_FORMAT_GRID] = GRID_EXTENSIONS,
  [GOOSEBERRY_FORMAT_CSR] = CSR_EXTENSIONS,
  [GOOSEBERRY_FORMAT_SVM] = SVM_EXTENSIONS,
  [GOOSEBERRY_FORMAT_POINT] = POINT_EXTENSIONS,
  [GOOSEBERRY_FORMAT_GRAPH] = GRAPH_EXTENSIONS,
  [GOOSEBERRY_FORMAT_CLU] = CLU_EXTENSIONS,
  [GOOSEBERRY_FORMAT_MATRIXMARKET] = MATRIXMARKET_EXTENSIONS
};




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static int __is_dense(
    int type)
{
  switch (type) {
    case GOOSEBERRY_FORMAT_RAW:
    case GOOSEBERRY_FORMAT_GRID:
      return 1;
    default:
      return 0;
  }
}


static int __get_file_type(
    const char * const name)
{
  size_t i,j;

  for (i=0;i<ARRAY_SIZE(FILE_TYPES);++i) {
    for (j=0;FILE_TYPES[i][j] != NULL;++j) {
      if (dl_string_endswith(name,FILE_TYPES[i][j])) {
        return i;
      }
    }
  }

  return -1;
}


static int __usage(
    const char * const name, 
    FILE * fout)
{
  size_t i;

  fprintf(fout,"USAGE:\n");
  fprintf(fout,"%s <command> [options]\n",name);
  fprintf(fout,"\n");
  fprintf(fout,"Commands:\n");

  for (i=0;i<NCOMMANDS;++i) {
    fprintf(fout,"\t%s : %s\n",COMMANDS[i].str,COMMANDS[i].desc);
  }

  return 1;
}


static int __command_usage(
    const char * const name, 
    const char * const cmd,
    const cmd_opt_t * const opts, 
    const size_t nopts, 
    FILE * fout)
{
  fprintf(stdout,"USAGE:\n");
  fprintf(stdout,"%s %s [options]\n",name,cmd);
  fprintf(stdout,"\n");
  fprint_cmd_opts(fout,opts,nopts);

  return 1;
}


/* COMMAND FUNCTIONS *********************************************************/


static int __help(
    int argc, 
    char ** argv)
{
  __usage(argv[0],stdout);
  return GOOSEBERRY_SUCCESS;
}


static int __analyze(
    int argc, 
    char ** argv)
{
  dl_timer_t io_tmr, op_tmr;
  size_t nargs, i;
  dim_t k, nec, ner, prows;
  ind_t j, nnz;
  real_t v, minvalue, maxvalue;
  dim_t * rowsize, * colsize, * perm, * order, * pk = NULL;
  int times, err, type, analysis;
  cmd_arg_t * args = NULL;
  matrix_t * mat = NULL;
  char const * matfile = NULL, * pfile = NULL;
  double nops;

  /* set defaults */
  times = 0;
  analysis = ANALYSIS_MATRIXSTATS;

  err = cmd_parse_args(argc-2,argv+2,ANALYSIS_OPTS,NANALYSIS_OPTS,&args, \
      &nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }

  err = GOOSEBERRY_SUCCESS;

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],ANALYSIS_OPTS,NANALYSIS_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case ANALYSIS_OPTION_HELP:
        __command_usage(argv[0],argv[1],ANALYSIS_OPTS,NANALYSIS_OPTS,stdout);
        goto END;
        break;
      case ANALYSIS_OPTION_TYPE:
        analysis = (analysis_t)args[i].val.o;
        break;
      case ANALYSIS_OPTION_INFILE:
        if (matfile == NULL) {
          matfile = args[i].val.s;
        } else {
          eprintf("Too many input files specified\n");
          err = GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        break;
      case ANALYSIS_OPTION_PERMFILE:
        pfile = args[i].val.s;
        break;
      case ANALYSIS_OPTION_TIME:
        times = 1;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (matfile == NULL) {
    eprintf("You must specify a matrix/vector input file.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_init_timer(&io_tmr);
    dl_init_timer(&op_tmr);
    dl_start_timer(&io_tmr);
  }

  /* read in input files */
  type = __get_file_type(matfile);
  if (type < 0) {
    eprintf("Unknown file format of '%s'\n",matfile);
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  } else {
    /* read in the matrix */
    mat = matrix_alloc(1);
    memset(mat,0,sizeof(matrix_t));
    if (__is_dense(type)) {
      err = gooseberry_read_dense_matrix(type,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_DENSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_DENSE_MATRIX;
      }
    } else {
      err = gooseberry_read_sparse_matrix(type,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowptr),&(mat->rowind),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_SPARSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_SPARSE_MATRIX;
      }
    }
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }

  if (pfile) {
    prows = mat->nrows;
    err = gooseberry_read_labels(pfile,&prows,&pk);
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
    if (prows != mat->nrows) {
      eprintf("Matrix is "PF_DIM_T"x"PF_DIM_T" but permutation file has " \
          PF_DIM_T" rows.\n",mat->nrows,mat->ncols,prows);
      err = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto END;
    }
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_start_timer(&op_tmr);
  }

  if (pk) {
    perm = dim_alloc(mat->nrows);
    order = dim_alloc(mat->nrows);
    dim_incset(order,0,1,mat->nrows);
    dd_countingsort_kv(pk,order,0,mat->nrows,mat->nrows,perm,NULL);
    matrix_permute(mat,perm,perm);
    dl_free(order);
    dl_free(perm);
  }

  nnz = mat->rowptr[mat->nrows];

  switch (analysis) {
    case ANALYSIS_MATRIXSTATS:
      if (nnz > 0) {
        minvalue = maxvalue = mat->rowval[0];
      } else {
        minvalue = maxvalue = 0;
      }
      rowsize = dim_alloc(mat->nrows);
      colsize = dim_init_alloc(0,mat->ncols);  
      ner = 0;
      for (i=0;i<mat->nrows;++i) {
        rowsize[i] = mat->rowptr[i+1] - mat->rowptr[i];
        if (rowsize[i] == 0) {
          ++ner;
        }
      }
      for (i=0;i<mat->nrows;++i) {
        for (j=mat->rowptr[i];j<mat->rowptr[i+1];++j) {
          k = mat->rowind[j];
          v = mat->rowval[j];
          ++colsize[k];
          if (v < minvalue) {
            minvalue = v;
          }
          if (v > maxvalue) {
            maxvalue = v;
          }
        }
      }
      nec = 0;
      for (i=0;i<mat->ncols;++i) {
        if (colsize[i] == 0) {
          ++nec;
        }
      }
      printf("Number of rows          = %16zu\n",(size_t)mat->nrows);
      printf("Number of columns       = %16zu\n",(size_t)mat->ncols);
      printf("Number of non-zeros     = %16zu\n",(size_t)mat->rowptr[mat->nrows]);
      printf("Median nnz / row        = %16zu\n", \
          (size_t)dim_median(rowsize,mat->nrows));
      printf("Mean nnz / row          = %16.3lf\n", \
          dim_arithmetic_mean(rowsize,mat->nrows));
      printf("Median nnz / column     = %16zu\n", \
          (size_t)dim_median(colsize,mat->ncols));
      printf("Mean nnz / column       = %16.3lf\n", \
          dim_arithmetic_mean(colsize,mat->ncols));
      printf("Maximum value           = %16.3lf\n",maxvalue);
      printf("Minimum value           = %16.3lf\n",minvalue);
      printf("Number of empty rows    = %16zu\n",(size_t)ner);
      printf("Number of empty columns = %16zu\n",(size_t)nec);
      dl_free(rowsize);
      dl_free(colsize);
      break;
    case ANALYSIS_CHOLESKY:
      analyze_cholesky(mat->nrows,mat->rowptr,mat->rowind,&nnz,&nops);
      printf("Number of non-zeroes = "PF_IND_T"\n",nnz);
      printf("Number of operations = %g\n",nops); 
      break;
    default:
      eprintf("Unknown analysis '%d'\n",analysis);
      err = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto END;
  }

  if (times) {
    dl_stop_timer(&op_tmr);
  } 

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_print_header("Times",'#');
    printf(" I/O: %0.04lf\n",dl_poll_timer(&io_tmr));
    printf(" Compute: %0.04lf\n",dl_poll_timer(&op_tmr));
    dl_print_footer('#');
  }

  END:

  if (pk) {
    dl_free(pk);
  }
  if (mat) {
    matrix_free(mat);
  } 
  if (args) {
    dl_free(args); 
  }

  return err;
}


static int __permute(
    int argc, 
    char ** argv)
{
  dl_timer_t io_tmr, op_tmr;
  size_t nargs,i;
  int times, j, err, permutation;
  dim_t prows;
  cmd_arg_t * args = NULL;
  dim_t * rpk = NULL, * cpk = NULL, *rperm = NULL, *cperm = NULL, *order;
  matrix_t * mat = NULL;
  const char * matfile = NULL, * outfile = NULL, * rpf = NULL, * cpf = NULL;

  /* set defaults */
  times = 0;
  permutation = PERMUTE_PERMUTATION_RANDOM;

  err = cmd_parse_args(argc-2,argv+2,PERMUTE_OPTS,NPERMUTE_OPTS,&args,&nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }

  err = GOOSEBERRY_SUCCESS;

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],PERMUTE_OPTS,NPERMUTE_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case PERMUTE_OPTION_HELP:
        __command_usage(argv[0],argv[1],PERMUTE_OPTS,NPERMUTE_OPTS,stdout);
        goto END;
        break;
      case PERMUTE_OPTION_INFILE:
        if (matfile == NULL) {
          matfile = args[i].val.s;
        } else {
          eprintf("Too many input files specified\n");
          err = GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        break;
      case PERMUTE_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case PERMUTE_OPTION_TIME:
        times = 1;
        break;
      case PERMUTE_OPTION_ROWPERM:
        rpf = args[i].val.s;
        break;
      case PERMUTE_OPTION_COLPERM:
        cpf = args[i].val.s;
        break;
      case PERMUTE_OPTION_PERMUTATION:
        permutation = (permute_permutation_t)args[i].val.o;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (matfile == NULL) {
    eprintf("You must specify a matrix/vector input file.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (outfile == NULL) {
    eprintf("You must specify an output matrix/vector file.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (permutation == PERMUTE_PERMUTATION_FILE) {
    if (rpf == NULL && cpf == NULL) {
      eprintf("You must specify a row permutation and/or a column permutation "
          "to permute from a file.\n");
      err = GOOSEBERRY_ERROR_INVALIDINPUT;
    }
  } else {
    if (rpf || cpf) {
      eprintf("Input row and column permutation files are only for use with "
          "file based permutations.\n");
      err = GOOSEBERRY_ERROR_INVALIDINPUT;
    }
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_init_timer(&io_tmr);
    dl_init_timer(&op_tmr);
    dl_start_timer(&io_tmr);
  }

  /* read in input files */
  j = __get_file_type(matfile);
  if (j < 0) {
    eprintf("Unknown file format of '%s'\n",matfile);
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  } else {
    /* read in the matrix */
    mat = matrix_calloc(1);
    if (__is_dense(j)) {
      err = gooseberry_read_dense_matrix(j,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_DENSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_DENSE_MATRIX;
      }
    } else {
      err = gooseberry_read_sparse_matrix(j,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowptr),&(mat->rowind),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_SPARSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_SPARSE_MATRIX;
      }
    }
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }

  if (permutation == PERMUTE_PERMUTATION_RANDOM) {
    if (mat->nrows != mat->ncols) {
      eprintf("Cannot apply a single permutation to columns and rows of a "
          "non-square matrix\n");
      err = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto END;
    }
  }

  /* read in permutation files if provided */
  if (permutation == PERMUTE_PERMUTATION_FILE) {
    if (rpf) {
      prows = mat->nrows;
      err = gooseberry_read_labels(rpf,&prows,&rpk);
      if (err != GOOSEBERRY_SUCCESS) {
        goto END;
      }
    }
    if (cpf) {
      prows = mat->ncols;
      err = gooseberry_read_labels(cpf,&prows,&cpk);
      if (err != GOOSEBERRY_SUCCESS) {
        goto END;
      }
    }
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_start_timer(&op_tmr);
  }

  /* permute the input matrices */
  switch (permutation) {
    case PERMUTE_PERMUTATION_FILE:
      if (rpk) {
        rperm = dim_alloc(mat->nrows);
        order = dim_alloc(mat->nrows);
        dim_incset(order,0,1,mat->nrows);
        dd_countingsort_kv(rpk,order,0,mat->nrows,mat->nrows,rperm,
            &mat->rdist);
        dl_free(order);
        dl_free(rpk);
        rpk = NULL;
      } 
      if (cpk) {
        cperm = dim_alloc(mat->ncols);
        order = dim_alloc(mat->ncols);
        dim_incset(order,0,1,mat->ncols);
        dd_countingsort_kv(cpk,order,0,mat->ncols,mat->ncols,cperm,NULL);
        dl_free(order);
        dl_free(cpk);
        cpk = NULL;
      } 
      break;
    case PERMUTE_PERMUTATION_RANDOM:
      rperm = dim_alloc(mat->nrows);
      dim_incset(rperm,0,1,mat->nrows);
      dim_pseudo_shuffle(rperm,mat->nrows/4,mat->nrows);
      cperm = dim_duplicate(rperm,mat->nrows);
      break;
    case PERMUTE_PERMUTATION_ROWRANDOM:
      rperm = dim_alloc(mat->nrows);
      dim_incset(rperm,0,1,mat->nrows);
      dim_pseudo_shuffle(rperm,mat->nrows/4,mat->nrows);
      break;
    case PERMUTE_PERMUTATION_COLRANDOM:
      cperm = dim_alloc(mat->ncols);
      dim_incset(cperm,0,1,mat->ncols);
      dim_pseudo_shuffle(cperm,mat->ncols/4,mat->ncols);
      break;
    case PERMUTE_PERMUTATION_BANDWIDTH:
      rperm = dim_alloc(mat->nrows);
      if ((err = permute_cuthillmckee(mat->nrows,mat->ncols,mat->rowptr, \
          mat->rowind,mat->rowval,NULL,0,rperm)) != GOOSEBERRY_SUCCESS) {
        goto END;
      }
      cperm = dim_duplicate(rperm,mat->nrows);
      break;
  }
  matrix_permute(mat,rperm,cperm);

  if (times) {
    dl_stop_timer(&op_tmr);
    dl_start_timer(&io_tmr);
  } 

  /* save the output */
  j = __get_file_type(outfile);
  if (__is_dense(j)) {
    if (mat->type != MATRIX_TYPE_DENSE_VECTOR &&
        mat->type != MATRIX_TYPE_DENSE_MATRIX) {
      matrix_densify(mat);
    }
    err = gooseberry_write_dense_matrix(j,outfile,mat->nrows,mat->ncols,
        mat->rowval);
  } else {
    if (mat->type != MATRIX_TYPE_SPARSE_VECTOR &&
        mat->type != MATRIX_TYPE_SPARSE_MATRIX) {
      matrix_sparsify(mat);
    }
    err = gooseberry_write_sparse_matrix(j,outfile,mat->nrows,mat->ncols,
        mat->rowptr,mat->rowind,mat->rowval);
  }

  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_print_header("Times",'#');
    printf(" I/O: %0.04lf\n",dl_poll_timer(&io_tmr));
    printf(" Compute: %0.04lf\n",dl_poll_timer(&op_tmr));
    dl_print_footer('#');
  }

  END:

  if (mat) {
    matrix_free(mat);
  } 
  if (rperm) {
    dl_free(rperm);
  }
  if (rpk) {
    dl_free(rpk);
  }
  if (cperm) {
    dl_free(cperm);
  }
  if (cpk) {
    dl_free(cpk);
  }
  if (args) {
    dl_free(args); 
  }

  return err;
}



static int __transform(
    int argc, 
    char ** argv)
{
  dl_timer_t io_tmr, op_tmr;
  size_t nargs, i;
  int times, j, err, operation;
  ind_t offset;
  dim_t prows, pcols, nout, p, nparts;
  cmd_arg_t * args = NULL;
  matrix_t * mat = NULL, * out = NULL;
  dim_t * dist = NULL, * map = NULL;
  char * sfile;
  char const * matfile = NULL, * outfile = NULL, * partfile = NULL;

  /* set defaults */
  nparts = 0;
  nout = 0;
  times = 0;
  pcols = prows = 0;
  operation = -1;

  err = cmd_parse_args(argc-2,argv+2,TRANSFORM_OPTS,NTRANSFORM_OPTS,&args, \
      &nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }

  err = GOOSEBERRY_SUCCESS;

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],TRANSFORM_OPTS,NTRANSFORM_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case TRANSFORM_OPTION_HELP:
        __command_usage(argv[0],argv[1],TRANSFORM_OPTS,NTRANSFORM_OPTS,stdout);
        goto END;
        break;
      case TRANSFORM_OPTION_INFILE:
        if (matfile == NULL) {
          matfile = args[i].val.s;
        } else {
          eprintf("Extra input file specified: '%s'\n",args[i].val.s);
          err = GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        break;
      case TRANSFORM_OPTION_OUTFILE:
        if (outfile == NULL) {
          outfile = args[i].val.s;
        } else {
          eprintf("Extra output file specified: '%s'\n",args[i].val.s);
          err = GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        break;
      case TRANSFORM_OPTION_TIME:
        times = 1;
        break;
      case TRANSFORM_OPTION_OPERATION:
        operation = (transform_operation_t)args[i].val.o;
        break;
      case TRANSFORM_OPTION_PARTFILE:
        if (partfile == NULL) {
          partfile = args[i].val.s;
        } else {
          eprintf("Extra part file specified: '%s'\n",args[i].val.s);
          err = GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (matfile == NULL) {
    eprintf("You must specify a matrix/vector input file.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (outfile == NULL) {
    eprintf("You must specify an output matrix/vector file.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if ((operation == TRANSFORM_OPERATION_ROWSPLIT || \
      operation == TRANSFORM_OPERATION_COLSPLIT) && \
      (partfile == NULL && nparts == 0)) {
    eprintf("You must specify a part file for splitting on or a number of "
        "partitions.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_init_timer(&io_tmr);
    dl_init_timer(&op_tmr);
    dl_start_timer(&io_tmr);
  }

  /* read in input file */
  j = __get_file_type(matfile);
  if (j < 0) {
    eprintf("Unknown file format of '%s'\n",matfile);
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  } else {
    /* read in the matrix */
    mat = matrix_calloc(1);
    if (__is_dense(j)) {
      err = gooseberry_read_dense_matrix(j,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_DENSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_DENSE_MATRIX;
      }
    } else {
      err = gooseberry_read_sparse_matrix(j,matfile,&(mat->nrows), \
          &(mat->ncols),&(mat->rowptr),&(mat->rowind),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_SPARSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_SPARSE_MATRIX;
      }
    }
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }
  if (operation == TRANSFORM_OPERATION_ROWSPLIT) {
    if (partfile != NULL) {
      err = gooseberry_read_partition(partfile,&prows,&nout,&map,NULL,&dist);
      if (err != GOOSEBERRY_SUCCESS) {
        goto END;
      } else if (prows != mat->nrows) {
        eprintf("Invalid number of rows in partition file: found '"PF_DIM_T"' "
            "but matrix has '"PF_DIM_T"'\n",prows,mat->nrows);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
      }
    } else {
      dl_error("Unimplemented\n"); 
    }
  } else if (operation == TRANSFORM_OPERATION_COLSPLIT) {
    if (partfile) {
      err = gooseberry_read_partition(partfile,&pcols,&nout,&map,NULL,&dist);
      if (err != GOOSEBERRY_SUCCESS) {
        goto END;
      } else if (pcols != mat->ncols) {
        eprintf("Invalid number of columns in partition file: found '"PF_DIM_T
            "' but matrix has '"PF_DIM_T"'\n",pcols,mat->ncols);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
      }
    } else {
      dl_error("Unimplemented\n"); 
    }
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_start_timer(&op_tmr);
  }

  out = matrix_calloc(1);
  out->type = MATRIX_TYPE_SPARSE_MATRIX;
  switch (operation) {
    case TRANSFORM_OPERATION_CONVERT:
      out->nrows = mat->nrows;
      out->ncols = mat->ncols;
      if (mat->type == out->type) {
        out->rowptr = ind_duplicate(mat->rowptr,mat->nrows+1);
        out->rowind = dim_duplicate(mat->rowind,mat->rowptr[mat->nrows]);
        out->rowval = real_duplicate(mat->rowval,mat->rowptr[mat->nrows]);
      } else {
        dl_error("Not finished yet");
      }
      break;
    case TRANSFORM_OPERATION_SYMMETRIFY:
      out->ncols = out->nrows = dl_max(mat->nrows,mat->ncols);
      out->type = MATRIX_TYPE_SPARSE_MATRIX;
      err = gooseberry_symmetrify_sparse(mat->nrows,mat->ncols,mat->rowptr,
          mat->rowind,mat->rowval,&out->rowptr,&out->rowind,&out->rowval);
      break;
    case TRANSFORM_OPERATION_DEBIPARTIFY:
      out->ncols = out->nrows = mat->nrows+mat->ncols;
      out->type = MATRIX_TYPE_SPARSE_MATRIX;
      err = gooseberry_debipartify_sparse(mat->nrows,mat->ncols,mat->rowptr,
          mat->rowind,mat->rowval,&out->rowptr,&out->rowind,&out->rowval);
      break;
    case TRANSFORM_OPERATION_TRANSPOSE:
      out->ncols = mat->nrows;
      out->nrows = mat->ncols;
      err = gooseberry_transpose_sparse(mat->nrows,mat->ncols,mat->rowptr,
          mat->rowind,mat->rowval,&out->rowptr,&out->rowind,&out->rowval);
      break;
    case TRANSFORM_OPERATION_ROWSPLIT:
      out->ncols = mat->ncols;
      out->nrows = mat->nrows;
      gooseberry_rowsplit_sparse(mat->nrows,mat->ncols,mat->rowptr,mat->rowind,
          mat->rowval,nout,dist,map,&out->rowptr,&out->rowind,&out->rowval);
      break;
    case TRANSFORM_OPERATION_COLSPLIT:
      out->ncols = 0;
      for (p=0;p<nout;++p) {
        if (dist[p+1] - dist[p] > out->ncols)  {
          out->ncols = dist[p+1] - dist[p];
        }
      }
      out->nrows = mat->nrows;
      gooseberry_colsplit_sparse(mat->nrows,mat->ncols,mat->rowptr,mat->rowind,
          mat->rowval,nout,dist,map,&out->rowptr,&out->rowind,&out->rowval);
      break;
    default:
      eprintf("Unknown transform operation: '%d'\n",operation);
      err = GOOSEBERRY_ERROR_INVALIDINPUT;
      break;
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_stop_timer(&op_tmr);
    dl_start_timer(&io_tmr);
  } 

  j = __get_file_type(outfile);
  if (operation == TRANSFORM_OPERATION_COLSPLIT || 
      operation == TRANSFORM_OPERATION_ROWSPLIT) {
    /* long enough to append a 64 bit number and null character */
    sfile = malloc(strlen(outfile)+22);
    for (p=0;p<nout;++p) {
      sprintf(sfile,"%s."PF_DIM_T,outfile,p);
      if (operation == TRANSFORM_OPERATION_COLSPLIT) {
        prows = mat->nrows;
        pcols = dist[p+1] - dist[p];
        offset = p*mat->nrows;
      } else {
        prows = dist[p+1] - dist[p];
        pcols = mat->ncols;
        offset = dist[p];
      }
      if (__is_dense(j)) {
        offset = mat->nrows*mat->ncols*p;
        if (out->type != MATRIX_TYPE_DENSE_VECTOR && 
            out->type != MATRIX_TYPE_DENSE_MATRIX) {
          matrix_densify(out);
        }
        err = gooseberry_write_dense_matrix(j,sfile,prows,pcols,
            out->rowval+offset);
      } else {
        if (out->type != MATRIX_TYPE_SPARSE_VECTOR && 
            out->type != MATRIX_TYPE_SPARSE_MATRIX) {
          matrix_sparsify(out);
        }
        err = gooseberry_write_sparse_matrix(j,sfile,prows,pcols,
            out->rowptr+offset,out->rowind,out->rowval);
      }
    }
    dl_free(sfile);
  } else {
    /* save the output */
    if (__is_dense(j)) {
      if (out->type != MATRIX_TYPE_DENSE_VECTOR && 
          out->type != MATRIX_TYPE_DENSE_MATRIX) {
        matrix_densify(out);
      }
      err = gooseberry_write_dense_matrix(j,outfile,out->nrows,out->ncols,
          out->rowval);
    } else {
      if (out->type != MATRIX_TYPE_SPARSE_VECTOR && 
          out->type != MATRIX_TYPE_SPARSE_MATRIX) {
        matrix_sparsify(out);
      }
      err = gooseberry_write_sparse_matrix(j,outfile,out->nrows,out->ncols,
          out->rowptr,out->rowind,out->rowval);
    }
  }

  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_print_header("Times",'#');
    printf(" I/O: %0.04lf\n",dl_poll_timer(&io_tmr));
    printf(" Compute: %0.04lf\n",dl_poll_timer(&op_tmr));
    dl_print_footer('#');
  }

  END:

  if (mat) {
    matrix_free(mat);
  } 
  if (out) {
    matrix_free(out);
  }
  if (args) {
    dl_free(args); 
  }

  return err;
}


static int __generate(
    int argc, 
    char ** argv)
{
  dl_timer_t io_tmr, op_tmr;
  size_t nargs, i;
  int times, j, err, type;
  dim_t size;
  cmd_arg_t * args = NULL;
  matrix_t * out = NULL;
  const char * outfile = NULL;

  /* set defaults */
  size = 0;
  times = 0;
  type = GENERATE_TYPE_NULL;

  err = cmd_parse_args(argc-2,argv+2,GENERATE_OPTS,NGENERATE_OPTS,&args,
      &nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }

  err = GOOSEBERRY_SUCCESS;

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],GENERATE_OPTS,NGENERATE_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case GENERATE_OPTION_HELP:
        __command_usage(argv[0],argv[1],GENERATE_OPTS,NGENERATE_OPTS,stdout);
        goto END;
        break;
      case GENERATE_OPTION_OUTFILE:
        if (outfile == NULL) {
          outfile = args[i].val.s;
        } else {
          eprintf("Extra output file specified: '%s'\n",args[i].val.s);
          err = GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        break;
      case GENERATE_OPTION_TIME:
        times = 1;
        break;
      case GENERATE_OPTION_TYPE:
        type = (generate_type_t)args[i].val.o;
        break;
      case GENERATE_OPTION_SIZE:
        size = (dim_t)args[i].val.i;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (size == 0) {
    eprintf("You must specify a size greater than zero.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (type == GENERATE_TYPE_NULL) {
    eprintf("You must specify a type to generate.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (outfile == NULL) {
    eprintf("You must specify an output matrix/vector file.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_init_timer(&io_tmr);
    dl_init_timer(&op_tmr);
    dl_start_timer(&op_tmr);
  }

  out = matrix_calloc(1);
  switch (type) {
    case GENERATE_TYPE_DENSE_VECTOR:
      out->type = MATRIX_TYPE_DENSE_VECTOR;
      out->nrows = size;
      out->ncols = 1;
      out->rowval = real_alloc(out->nrows);
      real_fill_rand(-1.0,1.0,out->rowval,out->nrows);
      break;
    default:
      eprintf("Unknown generate type: '%d'\n",type);
      err = GOOSEBERRY_ERROR_INVALIDINPUT;
      break;
  }

  if (times) {
    dl_stop_timer(&op_tmr);
    dl_start_timer(&io_tmr);
  } 

  j = __get_file_type(outfile);
  if (__is_dense(j)) {
    if (out->type != MATRIX_TYPE_DENSE_VECTOR && 
        out->type != MATRIX_TYPE_DENSE_MATRIX) {
      matrix_densify(out);
    }
    err = gooseberry_write_dense_matrix(j,outfile,out->nrows,out->ncols,
        out->rowval);
  } else {
    if (out->type != MATRIX_TYPE_SPARSE_VECTOR && 
        out->type != MATRIX_TYPE_SPARSE_MATRIX) {
      matrix_sparsify(out);
    }
    err = gooseberry_write_sparse_matrix(j,outfile,out->nrows,out->ncols,
        out->rowptr,out->rowind,out->rowval);
  }

  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_print_header("Times",'#');
    printf(" I/O: %0.04lf\n",dl_poll_timer(&io_tmr));
    printf(" Compute: %0.04lf\n",dl_poll_timer(&op_tmr));
    dl_print_footer('#');
  }

  END:

  if (out) {
    matrix_free(out);
  }
  if (args) {
    dl_free(args); 
  }

  return err;
}


static int __blas(
    int argc, 
    char ** argv)
{
  dl_timer_t io_tmr, op_tmr, aux_tmr;
  size_t nargs, runs, r, i, ninfiles = 0,nthreads;
  int times, j, err, oper, redband;
  dim_t outrows, outcols, prows;
  cmd_arg_t * args = NULL;
  dim_t *rperm = NULL, *cperm = NULL, *bperm = NULL; 
  matrix_t * in[256], * out = NULL;
  const char * infiles[256], * outfile = NULL, * rpf = NULL, * cpf = NULL;

  /* set defaults */
  redband = 0;
  times = 0;
  runs = 1;
  oper = BLAS_OPERATION_NOOP;

  err = cmd_parse_args(argc-2,argv+2,BLAS_OPTS,NBLAS_OPTS,&args,&nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }

  err = GOOSEBERRY_SUCCESS;

  nthreads = DEFAULT_NUMTHREADS;

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],BLAS_OPTS,NBLAS_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case BLAS_OPTION_HELP:
        __command_usage(argv[0],argv[1],BLAS_OPTS,NBLAS_OPTS,stdout);
        goto END;
        break;
      case BLAS_OPTION_OPERATION:
        oper = (blas_operation_t)args[i].val.o;
        break;
      case BLAS_OPTION_INFILE:
        infiles[ninfiles++] = args[i].val.s;
        break;
      case BLAS_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case BLAS_OPTION_TIME:
        times = 1;
        break;
      case BLAS_OPTION_RUNS:
        runs = (size_t)args[i].val.i; 
        break;
#ifndef NO_OMP
      case BLAS_OPTION_THREADS:
        nthreads = (size_t)args[i].val.i;
        omp_set_num_threads(nthreads);
        break;
#endif
      case BLAS_OPTION_ROWPERM:
        rpf = args[i].val.s;
        break;
      case BLAS_OPTION_COLPERM:
        cpf = args[i].val.s;
        break;
      case BLAS_OPTION_REDUCEBANDWIDTH:
        redband = 1;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_init_timer(&io_tmr);
    dl_init_timer(&op_tmr);
    dl_init_timer(&aux_tmr);
    dl_start_timer(&io_tmr);
  }

  /* read in input files */
  for (i=0;i<ninfiles;++i) {
    j = __get_file_type(infiles[i]);
    if (j < 0) {
      eprintf("Unknown file format of '%s'\n",infiles[i]);
      err = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto END;
    } else {
      /* read in the matrix/vector */
      in[i] = matrix_calloc(1);
      if (__is_dense(j)) {
        err = gooseberry_read_dense_matrix(j,infiles[i],&(in[i]->nrows),
            &(in[i]->ncols),&(in[i]->rowval));
        if (in[i]->ncols == 1 || in[i]->nrows == 1) {
          in[i]->type = MATRIX_TYPE_DENSE_VECTOR;
        } else {
          in[i]->type = MATRIX_TYPE_DENSE_MATRIX;
        }
      } else {
        err = gooseberry_read_sparse_matrix(j,infiles[i],&(in[i]->nrows),
            &(in[i]->ncols),&(in[i]->rowptr),&(in[i]->rowind),
            &(in[i]->rowval));
        if (in[i]->ncols == 1 || in[i]->nrows == 1) {
          in[i]->type = MATRIX_TYPE_SPARSE_VECTOR;
        } else {
          in[i]->type = MATRIX_TYPE_SPARSE_MATRIX;
        }
      }
      if (err != GOOSEBERRY_SUCCESS) {
        goto END;
      }
    }
  }

  /* read in permutation files if provided */
  if (rpf) {
    prows = in[0]->nrows;
    err = gooseberry_read_partition(rpf,&prows,&in[0]->nrdist,NULL,&rperm,
        &in[0]->rdist);
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
    printf("Read "PF_DIM_T"-way row partition\n",in[0]->nrdist);
  } else {
    in[0]->nrdist = 1;
  }
  if (cpf) {
    prows = in[0]->ncols;
    err = gooseberry_read_partition(cpf,&prows,&in[0]->ncdist,NULL,&cperm,
        &in[0]->cdist);
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
    printf("Read "PF_DIM_T"-way column partition\n",in[0]->ncdist);
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_start_timer(&aux_tmr);
  }

  matrix_permute(in[0],rperm,cperm);
  matrix_permute(in[1],cperm,NULL);

  if (redband) {
    bperm = dim_alloc(in[0]->nrows); 
    if ((err = permute_cuthillmckee(in[0]->nrows,in[0]->ncols,in[0]->rowptr,
          in[0]->rowind,in[0]->rowval,in[0]->rdist,in[0]->nrdist,bperm)) 
        != GOOSEBERRY_SUCCESS) {
      goto END;
    }
    if (rperm && cperm) {
      for (i=0;i<in[0]->nrows;++i) {
        rperm[i] = bperm[rperm[i]];
      }
      for (i=0;i<in[0]->nrows;++i) {
        cperm[i] = bperm[cperm[i]];
      }
    }
    matrix_permute(in[0],bperm,bperm);
    matrix_permute(in[1],bperm,NULL);
    dl_free(bperm);
  }

  /* allocate the output matrix */
  outrows = in[0]->nrows;
  outcols = in[ninfiles-1]->ncols;
  out = matrix_alloc(1);
  j = __get_file_type(outfile);
  if (__is_dense(j)) {
    if (outrows == 1|| outcols == 1) {
      matrix_init(MATRIX_TYPE_DENSE_VECTOR,outrows,outcols,0,out);
    } else {
      matrix_init(MATRIX_TYPE_DENSE_MATRIX,outrows,outcols,0,out);
    }
  } else {
    if (outrows == 1|| outcols == 1) {
      matrix_init(MATRIX_TYPE_SPARSE_VECTOR,outrows,outcols,NULL_IND,out);
    } else {
      matrix_init(MATRIX_TYPE_SPARSE_MATRIX,outrows,outcols,NULL_IND,out);
    }
  }

  if (times) {
    dl_stop_timer(&aux_tmr);
    dl_start_timer(&op_tmr);
  }

  for (r=0;r<runs;++r) {
    /* perform operation */
    switch (oper) {
      case BLAS_OPERATION_MULTIPLY:
        if (in[0]->ncols > in[1]->nrows) {
          eprintf("Matrix dimensions do not match for multiplication: "
              PF_DIM_T"x"PF_DIM_T" and "PF_DIM_T"x"PF_DIM_T"\n",in[0]->nrows,
              in[0]->ncols,in[1]->nrows,in[1]->ncols);
          err = GOOSEBERRY_ERROR_INVALIDINPUT;
          goto END;
        }
        switch(in[0]->type) {
          case MATRIX_TYPE_SPARSE_VECTOR:
          case MATRIX_TYPE_SPARSE_MATRIX:
            switch (in[1]->type) {
              case MATRIX_TYPE_DENSE_VECTOR:
              case MATRIX_TYPE_DENSE_MATRIX:
                matrix_buildindex(in[1]);
                if ((err = blas_spmult(in[0]->nrows,in[0]->ncols,in[1]->ncols,
                        in[0]->rowptr,in[0]->rowind,in[0]->rowval,
                        in[1]->colval,out->rowval,in[0]->rdist,in[0]->nrdist)) 
                    != GOOSEBERRY_SUCCESS) {
                  goto END; 
                }
                break;
              case MATRIX_TYPE_SPARSE_VECTOR:
              case MATRIX_TYPE_SPARSE_MATRIX:
                switch (out->type) {
                  case MATRIX_TYPE_DENSE_VECTOR:
                  case MATRIX_TYPE_DENSE_MATRIX:
                    if ((err = blas_spmultsp(in[0]->nrows,in[0]->ncols,
                            in[1]->ncols,in[0]->rowptr,in[0]->rowind,
                            in[0]->rowval,in[1]->rowptr,in[1]->rowind,
                            in[1]->rowval,out->rowval,in[0]->rdist,
                            in[0]->nrdist)) != GOOSEBERRY_SUCCESS) {
                      goto END; 
                    }
                    break;
                  case MATRIX_TYPE_SPARSE_VECTOR:
                  case MATRIX_TYPE_SPARSE_MATRIX:
                    if ((err = blas_spmultsp_sp(in[0]->nrows,in[0]->ncols,
                            in[1]->ncols,in[0]->rowptr,in[0]->rowind,
                            in[0]->rowval,in[1]->rowptr,in[1]->rowind,
                            in[1]->rowval,&out->rowptr,&out->rowind,
                            &out->rowval,in[0]->rdist,in[0]->nrdist)) != 
                        GOOSEBERRY_SUCCESS) {
                      goto END; 
                    }
                    break;
                  default:
                    eprintf("Unsupported output matrix type: %d\n",out->type);
                    err = GOOSEBERRY_ERROR_INVALIDINPUT;
                    break;
                }
                break;
              default:
                eprintf("Unsupported matrix combinations\n");
                err = GOOSEBERRY_ERROR_INVALIDINPUT;
                break;
            }
            break;
          case MATRIX_TYPE_DENSE_VECTOR:
          case MATRIX_TYPE_DENSE_MATRIX:
            switch (in[1]->type) {
              case MATRIX_TYPE_DENSE_VECTOR:
              case MATRIX_TYPE_DENSE_MATRIX:
                matrix_buildindex(in[1]);
                if ((err = blas_mult(in[0]->nrows,in[0]->ncols,in[1]->ncols,
                        in[0]->rowval,in[1]->colval,out->rowval,in[0]->rdist,
                        in[0]->nrdist)) != GOOSEBERRY_SUCCESS) {
                  goto END; 
                }
                break;
              case MATRIX_TYPE_SPARSE_VECTOR:
              case MATRIX_TYPE_SPARSE_MATRIX:
                eprintf("The operation mmsp is unsupported at the moment.\n");
                err = GOOSEBERRY_ERROR_UNIMPLEMENTED;
                break;
              default:
                eprintf("Unsupported matrix combinations\n");
                err = GOOSEBERRY_ERROR_INVALIDINPUT;
                break;
            }
            break;
          default:
            eprintf("Unknown matrix type '%d'\n",in[0]->type);
            err = GOOSEBERRY_ERROR_INVALIDINPUT;
            break;
        }
        break;
      default:
        eprintf("Unknown operation '%d'\n",oper);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
        break;
    }
  }

  if (times) {
    dl_stop_timer(&op_tmr);
    dl_start_timer(&aux_tmr);
  } 

  if (rperm) {
    matrix_unpermute(out,rperm,NULL);
  }

  if (times) {
    dl_stop_timer(&aux_tmr);
    dl_start_timer(&io_tmr);
  }

  /* save the output */
  j = __get_file_type(outfile);
  if (__is_dense(j)) {
    err = gooseberry_write_dense_matrix(j,outfile,out->nrows,out->ncols,
        out->rowval);
  } else {
    err = gooseberry_write_sparse_matrix(j,outfile,out->nrows,out->ncols,
        out->rowptr,out->rowind,out->rowval);
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_print_header("Times",'#');
    printf(" I/O: %0.04lf\n",dl_poll_timer(&io_tmr));
    printf(" Auxillary: %0.04lf\n",dl_poll_timer(&aux_tmr));
    printf(" Compute: %0.04lf\n",dl_poll_timer(&op_tmr));
    dl_print_footer('#');
  }

  END:

  if (out) {
    matrix_free(out);
  }
  if (rperm) {
    dl_free(rperm);
  }
  if (cperm) {
    dl_free(cperm);
  }
  for (i=0;i<ninfiles;++i) {
    if (in[i]) {
      matrix_free(in[i]);
      in[i] = NULL;
    }
  }
  if (args) {
    dl_free(args); 
  }

  return err;
}


static int __cgd(
    int argc, 
    char ** argv)
{
  dl_timer_t io_tmr, op_tmr;
  real_t error = 0;
  size_t nargs, runs, r, i, niter = 0;
  int times, j, err;
  dim_t outrows, outcols, prows;
  cmd_arg_t * args = NULL;
  dim_t * rpk = NULL, * cpk = NULL, *rperm = NULL, *cperm = NULL, *order;
  matrix_t * mat = NULL, * vec = NULL, * out = NULL;
  const char * matfile = NULL, * vecfile = NULL, * outfile = NULL, 
             * rpf = NULL, * cpf = NULL;

  /* set defaults */
  times = 0;
  runs = 1;

  err = cmd_parse_args(argc-2,argv+2,CGD_OPTS,NCGD_OPTS,&args,&nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }

  err = GOOSEBERRY_SUCCESS;

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],CGD_OPTS,NCGD_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case CGD_OPTION_HELP:
        __command_usage(argv[0],argv[1],CGD_OPTS,NCGD_OPTS,stdout);
        goto END;
        break;
      case CGD_OPTION_ERROR:
        error = (real_t)args[i].val.f;
        break;
      case CGD_OPTION_NITER:
        niter = (size_t)args[i].val.i;
        break;
      case CGD_OPTION_INFILE:
        if (matfile == NULL) {
          matfile = args[i].val.s;
        } else if (vecfile == NULL) {
          vecfile = args[i].val.s;
        } else {
          eprintf("Too many input files specified\n");
          err = GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        break;
      case CGD_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case CGD_OPTION_TIME:
        times = 1;
        break;
      case CGD_OPTION_RUNS:
        runs = (size_t)args[i].val.i; 
        break;
#ifndef NO_OMP
      case CGD_OPTION_THREADS:
        omp_set_num_threads(args[i].val.i);
        break;
#endif
      case CGD_OPTION_ROWPERM:
        rpf = args[i].val.s;
        break;
      case CGD_OPTION_COLPERM:
        cpf = args[i].val.s;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (matfile == NULL || vecfile == NULL) {
    eprintf("You must specify both a matrix input file and a vector input "
        "file (in that order).\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_init_timer(&io_tmr);
    dl_init_timer(&op_tmr);
    dl_start_timer(&io_tmr);
  }

  /* read in input files */
  j = __get_file_type(matfile);
  if (j < 0) {
    eprintf("Unknown file format of '%s'\n",matfile);
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  } else {
    /* read in the matrix */
    mat = matrix_calloc(1);
    if (__is_dense(j)) {
      err = gooseberry_read_dense_matrix(j,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_DENSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_DENSE_MATRIX;
      }
    } else {
      err = gooseberry_read_sparse_matrix(j,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowptr),&(mat->rowind),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_SPARSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_SPARSE_MATRIX;
      }
    }
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }
  j = __get_file_type(vecfile);
  if (j < 0) {
    eprintf("Unknown file format of '%s'\n",vecfile);
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  } else {
    /* read in the vector */
    vec = matrix_calloc(1);
    if (__is_dense(j)) {
      err = gooseberry_read_dense_matrix(j,vecfile,&(vec->nrows),
          &(vec->ncols),&(vec->rowval));
      if (vec->ncols == 1 || vec->nrows == 1) {
        vec->type = MATRIX_TYPE_DENSE_VECTOR;
      } else {
        vec->type = MATRIX_TYPE_DENSE_MATRIX;
      }
    } else {
      err = gooseberry_read_sparse_matrix(j,vecfile,&(vec->nrows),
          &(vec->ncols),&(vec->rowptr),&(vec->rowind),&(vec->rowval));
      if (vec->ncols == 1 || vec->nrows == 1) {
        vec->type = MATRIX_TYPE_SPARSE_VECTOR;
      } else {
        vec->type = MATRIX_TYPE_SPARSE_MATRIX;
      }
    }
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }

  printf("matrix = "PF_DIM_T"x"PF_DIM_T" vector = "PF_DIM_T"x"PF_DIM_T"\n",
      mat->nrows,mat->ncols,vec->nrows,vec->ncols);

  /* read in permutation files if provided */
  if (rpf) {
    prows = mat->nrows;
    err = gooseberry_read_labels(rpf,&prows,&rpk);
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }
  if (cpf) {
    prows = mat->ncols;
    err = gooseberry_read_labels(cpf,&prows,&cpk);
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }

  if (times) {
    dl_stop_timer(&io_tmr);
  }

  /* permute the input matrices */
  if (rpk) {
    rperm = dim_alloc(mat->nrows);
    order = dim_alloc(mat->nrows);
    dim_incset(order,0,1,mat->nrows);
    dd_countingsort_kv(rpk,order,0,mat->nrows,mat->nrows,rperm,
        &mat->rdist);
    dl_free(order);
    dl_free(rpk);
    rpk = NULL;
  } 
  if (cpk) {
    cperm = dim_alloc(mat->ncols);
    order = dim_alloc(mat->ncols);
    dim_incset(order,0,1,mat->ncols);
    dd_countingsort_kv(cpk,order,0,mat->ncols,mat->ncols,cperm,NULL);
    dl_free(order);
    dl_free(cpk);
    cpk = NULL;
  } 
  matrix_permute(mat,rperm,cperm);
  matrix_permute(vec,rperm,NULL);

  /* allocate the output matrix */
  outrows = mat->ncols;
  outcols = 1;
  out = matrix_alloc(1);
  j = __get_file_type(outfile);
  if (__is_dense(j)) {
    if (outrows == 1|| outcols == 1) {
      matrix_init(MATRIX_TYPE_DENSE_VECTOR,outrows,outcols,0,out);
    } else {
      matrix_init(MATRIX_TYPE_DENSE_MATRIX,outrows,outcols,0,out);
    }
  } else {
    if (outrows == 1|| outcols == 1) {
      matrix_init(MATRIX_TYPE_SPARSE_VECTOR,outrows,outcols,NULL_IND,out);
    } else {
      matrix_init(MATRIX_TYPE_SPARSE_MATRIX,outrows,outcols,NULL_IND,out);
    }
  }

  if (times) {
    dl_start_timer(&op_tmr);
  }

  for (r=0;r<runs;++r) {
    /* perform operation */
    if ((err = cgd(mat,vec,out,error,niter)) != GOOSEBERRY_SUCCESS) {
      goto END; 
    }
  }

  if (times) {
    dl_stop_timer(&op_tmr);
  } 

  if (cperm) {
    matrix_unpermute(out,cperm,NULL);
  }

  if (times) {
    dl_start_timer(&io_tmr);
  }

  /* save the output */
  err = gooseberry_write_dense_matrix(GOOSEBERRY_FORMAT_GRID,outfile,
      out->nrows,out->ncols,out->rowval);
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_print_header("Times",'#');
    printf(" I/O: %0.04lf\n",dl_poll_timer(&io_tmr));
    printf(" Compute: %0.04lf\n",dl_poll_timer(&op_tmr));
    dl_print_footer('#');
  }

  END:

  if (mat) {
    matrix_free(mat);
  } 
  if (vec) {
    matrix_free(vec);
  }
  if (out) {
    matrix_free(out);
  }
  if (rperm) {
    dl_free(rperm);
  }
  if (rpk) {
    dl_free(rpk);
  }
  if (cperm) {
    dl_free(cperm);
  }
  if (cpk) {
    dl_free(cpk);
  }
  if (args) {
    dl_free(args); 
  }

  return err;
}


static int __sgd(
    int argc, 
    char ** argv)
{
  int err = GOOSEBERRY_SUCCESS;
#ifdef XXX
  int times, j, err;
  dl_timer_t io_tmr, op_tmr;
  real_t error = 0;
  size_t nargs, runs, r, i, niter = 0;
  dim_t outrows, outcols, prows;
  cmd_arg_t * args = NULL;
  dim_t * rpk = NULL, * cpk = NULL, *rperm = NULL, *cperm = NULL, *order;
  matrix_t * mat = NULL, * vec = NULL, * out = NULL;
  const char * matfile = NULL, * ufile = NULL, * vfile, * rpf = NULL, 
             * cpf = NULL;

  /* set defaults */
  times = 0;
  runs = 1;

  err = cmd_parse_args(argc-2,argv+2,SGD_OPTS,NSGD_OPTS,&args,&nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],SGD_OPTS,NSGD_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case SGD_OPTION_HELP:
        __command_usage(argv[0],argv[1],SGD_OPTS,NSGD_OPTS,stdout);
        goto END;
        break;
      case SGD_OPTION_ERROR:
        error = (real_t)args[i].val.f;
        break;
      case SGD_OPTION_NITER:
        niter = (size_t)args[i].val.i;
        break;
      case SGD_OPTION_INFILE:
        if (matfile == NULL) {
          matfile = args[i].val.s;
        } else if (vecfile == NULL) {
          vecfile = args[i].val.s;
        } else {
          eprintf("Too many input files specified\n");
          err = GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        break;
      case SGD_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case SGD_OPTION_TIME:
        times = 1;
        break;
      case SGD_OPTION_RUNS:
        runs = (size_t)args[i].val.i; 
        break;
#ifndef NO_OMP
      case SGD_OPTION_THREADS:
        omp_set_num_threads(args[i].val.i);
        break;
#endif
      case SGD_OPTION_ROWPERM:
        rpf = args[i].val.s;
        break;
      case SGD_OPTION_COLPERM:
        cpf = args[i].val.s;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (matfile == NULL || vecfile == NULL) {
    eprintf("You must specify both a matrix input file and a vector input "
        "file (in that order).\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_init_timer(&io_tmr);
    dl_init_timer(&op_tmr);
    dl_start_timer(&io_tmr);
  }

  /* read in input files */
  j = __get_file_type(matfile);
  if (j < 0) {
    eprintf("Unknown file format of '%s'\n",matfile);
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  } else {
    /* read in the matrix */
    mat = matrix_calloc(1);
    if (__is_dense(j)) {
      err = gooseberry_read_dense_matrix(j,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_DENSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_DENSE_MATRIX;
      }
    } else {
      err = gooseberry_read_sparse_matrix(j,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowptr),&(mat->rowind),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_SPARSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_SPARSE_MATRIX;
      }
    }
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }
  j = __get_file_type(vecfile);
  if (j < 0) {
    eprintf("Unknown file format of '%s'\n",vecfile);
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  } else {
    /* read in the vector */
    vec = matrix_calloc(1);
    if (__is_dense(j)) {
      err = gooseberry_read_dense_matrix(j,vecfile,&(vec->nrows),
          &(vec->ncols),&(vec->rowval));
      if (vec->ncols == 1 || vec->nrows == 1) {
        vec->type = MATRIX_TYPE_DENSE_VECTOR;
      } else {
        vec->type = MATRIX_TYPE_DENSE_MATRIX;
      }
    } else {
      err = gooseberry_read_sparse_matrix(j,vecfile,&(vec->nrows),
          &(vec->ncols),&(vec->rowptr),&(vec->rowind),&(vec->rowval));
      if (vec->ncols == 1 || vec->nrows == 1) {
        vec->type = MATRIX_TYPE_SPARSE_VECTOR;
      } else {
        vec->type = MATRIX_TYPE_SPARSE_MATRIX;
      }
    }
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }

  printf("matrix = "PF_DIM_T"x"PF_DIM_T" vector = "PF_DIM_T"x"PF_DIM_T"\n",
      mat->nrows,mat->ncols,vec->nrows,vec->ncols);

  /* read in permutation files if provided */
  if (rpf) {
    prows = mat->nrows;
    err = gooseberry_read_labels(rpf,&prows,&rpk);
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }
  if (cpf) {
    prows = mat->ncols;
    err = gooseberry_read_labels(cpf,&prows,&cpk);
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }

  if (times) {
    dl_stop_timer(&io_tmr);
  }

  /* permute the input matrices */
  if (rpk) {
    rperm = dim_alloc(mat->nrows);
    order = dim_alloc(mat->nrows);
    dim_incset(order,0,1,mat->nrows);
    dim_countingsort_v(rpk,order,rperm,0,mat->nrows,mat->nrows);
    dl_free(order);
    dl_free(rpk);
    rpk = NULL;
  } 
  if (cpk) {
    cperm = dim_alloc(mat->ncols);
    order = dim_alloc(mat->ncols);
    dim_incset(order,0,1,mat->ncols);
    dim_countingsort_v(cpk,order,cperm,0,mat->ncols,mat->ncols);
    dl_free(order);
    dl_free(cpk);
    cpk = NULL;
  } 
  matrix_permute(mat,rperm,cperm);
  matrix_permute(vec,rperm,NULL);

  /* allocate the output matrix */
  outrows = mat->ncols;
  outcols = 1;
  out = matrix_alloc(1);
  j = __get_file_type(outfile);
  if (__is_dense(j)) {
    if (outrows == 1|| outcols == 1) {
      matrix_init(MATRIX_TYPE_DENSE_VECTOR,outrows,outcols,0,out);
    } else {
      matrix_init(MATRIX_TYPE_DENSE_MATRIX,outrows,outcols,0,out);
    }
  } else {
    if (outrows == 1|| outcols == 1) {
      matrix_init(MATRIX_TYPE_SPARSE_VECTOR,outrows,outcols,NULL_IND,out);
    } else {
      matrix_init(MATRIX_TYPE_SPARSE_MATRIX,outrows,outcols,NULL_IND,out);
    }
  }

  if (times) {
    dl_start_timer(&op_tmr);
  }

  for (r=0;r<runs;++r) {
    /* perform operation */
    if ((err = cgd(mat,vec,out,error,niter)) != GOOSEBERRY_SUCCESS) {
      goto END; 
    }
  }

  if (times) {
    dl_stop_timer(&op_tmr);
  } 

  if (cperm) {
    matrix_unpermute(out,cperm,NULL);
  }

  if (times) {
    dl_start_timer(&io_tmr);
  }

  /* save the output */
  err = gooseberry_write_dense_matrix(GOOSEBERRY_FORMAT_GRID,outfile,
      out->nrows,out->ncols,out->rowval);
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_print_header("Times",'#');
    printf(" I/O: %0.04lf\n",dl_poll_timer(&io_tmr));
    printf(" Compute: %0.04lf\n",dl_poll_timer(&op_tmr));
    dl_print_footer('#');
  }

  END:

  if (mat) {
    matrix_free(mat);
  } 
  if (vec) {
    matrix_free(vec);
  }
  if (out) {
    matrix_free(out);
  }
  if (rperm) {
    dl_free(rperm);
  }
  if (rpk) {
    dl_free(rpk);
  }
  if (cperm) {
    dl_free(cperm);
  }
  if (cpk) {
    dl_free(cpk);
  }
  if (args) {
    dl_free(args); 
  }

#endif
  return err;
}


static int __pagerank(
    int argc, 
    char ** argv)
{
  dl_timer_t io_tmr, op_tmr, pre_tmr, mul_tmr;
  size_t nargs,i,iter,runs,r;
  int times, j, err;
  ind_t l;
  dim_t prows, k, nsinks, m;
  real_t minerror, error, damping, diff, dist, deg, wgt;
  cmd_arg_t * args = NULL;
  real_t * rank = NULL, * indeg = NULL;
  dim_t * pk = NULL, * perm = NULL, * order = NULL, * sinks = NULL;
  matrix_t * mat = NULL, * out = NULL;
  const char * matfile = NULL, * outfile = NULL, * pf = NULL;

  /* set defaults */
  times = 0;
  runs = 1;
  minerror = 0.0;
  iter = 0;
  damping = 0.85;

  err = cmd_parse_args(argc-2,argv+2,PAGERANK_OPTS,NPAGERANK_OPTS,&args,
      &nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return GOOSEBERRY_ERROR_INVALIDINPUT;
  }

  err = GOOSEBERRY_SUCCESS;

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],PAGERANK_OPTS,NPAGERANK_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case PAGERANK_OPTION_HELP:
        __command_usage(argv[0],argv[1],PAGERANK_OPTS,NPAGERANK_OPTS,stdout);
        goto END;
        break;
      case PAGERANK_OPTION_INFILE:
        if (matfile == NULL) {
          matfile = args[i].val.s;
        } else {
          eprintf("Too many input files specified\n");
          err = GOOSEBERRY_ERROR_INVALIDINPUT;
        }
        break;
      case PAGERANK_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case PAGERANK_OPTION_TIME:
        times = 1;
        break;
      case PAGERANK_OPTION_PERM:
        pf = args[i].val.s;
        break;
      case PAGERANK_OPTION_RUNS:
        runs = (size_t)args[i].val.i;
        break;
      case PAGERANK_OPTION_NITER:
        iter = (size_t)args[i].val.i;
        break;
      case PAGERANK_OPTION_DAMPING:
        damping = (real_t)args[i].val.f;
        break;
      case PAGERANK_OPTION_ERROR:
        minerror = (real_t)args[i].val.f;
        break;
#ifndef NO_OMP
      case PAGERANK_OPTION_THREADS:
        omp_set_num_threads(args[i].val.i);
        break;
#endif
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (matfile == NULL) {
    eprintf("You must specify a matrix input file.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (outfile == NULL) {
    eprintf("You must specify an output vector file.\n");
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
  }
  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_init_timer(&io_tmr);
    dl_init_timer(&op_tmr);
    dl_init_timer(&pre_tmr);
    dl_init_timer(&mul_tmr);
    dl_start_timer(&io_tmr);
  }

  /* read in input files */
  j = __get_file_type(matfile);
  if (j < 0) {
    eprintf("Unknown file format of '%s'\n",matfile);
    err = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto END;
  } else {
    /* read in the matrix */
    mat = matrix_calloc(1);
    if (__is_dense(j)) {
      err = gooseberry_read_dense_matrix(j,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_DENSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_DENSE_MATRIX;
      }
    } else {
      err = gooseberry_read_sparse_matrix(j,matfile,&(mat->nrows),
          &(mat->ncols),&(mat->rowptr),&(mat->rowind),&(mat->rowval));
      if (mat->ncols == 1 || mat->nrows == 1) {
        mat->type = MATRIX_TYPE_SPARSE_VECTOR;
      } else {
        mat->type = MATRIX_TYPE_SPARSE_MATRIX;
      }
    }
    if (mat->nrows != mat->ncols) {
      if (mat->type == MATRIX_TYPE_DENSE_MATRIX) {
        eprintf("PageRank requires a square matrix, but input matrix '%s' is "
            PF_DIM_T"x"PF_DIM_T".\n",matfile,mat->nrows,mat->ncols);
        err = GOOSEBERRY_ERROR_INVALIDINPUT;
      } else {
        if (mat->nrows > mat->ncols) {
          mat->ncols = mat->nrows;
        } else if (mat->ncols > mat->nrows) {
          /* stretch the matrix */
          mat->rowptr = ind_realloc(mat->rowptr,mat->ncols+1);
          ind_set(mat->rowptr+mat->nrows+1,mat->rowptr[mat->nrows],
              mat->ncols-mat->nrows);
          mat->nrows = mat->ncols;
        }
      }
    }
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }

  /* read in permutation file if provided */
  if (pf) {
    prows = mat->nrows;
    err = gooseberry_read_labels(pf,&prows,&pk);
    if (err != GOOSEBERRY_SUCCESS) {
      goto END;
    }
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_start_timer(&pre_tmr);
  }

  /* permute the input matrix */
  if (pk) {
    perm = dim_alloc(mat->nrows);
    order = dim_alloc(mat->nrows);
    dim_incset(order,0,1,mat->nrows);
    dd_countingsort_kv(pk,order,0,mat->nrows,mat->nrows,perm,
        &mat->rdist);
    dl_free(order);
    dl_free(pk);
    pk = NULL;
    matrix_permute(mat,perm,perm);
  }

  /* normalize input matrix and find sinks */
  sinks = dim_alloc(mat->nrows);
  nsinks = 0;
  indeg = real_calloc(mat->ncols);
  switch (mat->type) {
    case MATRIX_TYPE_SPARSE_MATRIX:
      for (k=0;k<mat->nrows;++k) {
        if (mat->rowptr[k] == mat->rowptr[k+1]) {
          sinks[nsinks++] = k;
        } else {
          for (l=mat->rowptr[k];l<mat->rowptr[k+1];++l) {
            indeg[mat->rowind[l]] += mat->rowval[l];
          }
        }
      }
      for (k=0;k<mat->nrows;++k) {
        for (l=mat->rowptr[k];l<mat->rowptr[k+1];++l) {
          mat->rowval[l] /= indeg[mat->rowind[l]];
        }
      }
      break;
    case MATRIX_TYPE_DENSE_MATRIX:
      for (k=0;k<mat->nrows;++k) {
        deg = 0;
        for (m=0;m<mat->ncols;++m) {
          wgt = mat->rowval[(k*mat->ncols)+m];
          if (wgt != 0) {
            deg += wgt;
            indeg[m] += wgt;
          }
        }
        if (deg == 0) {
          sinks[nsinks++] = k;
        }
      }
      for (k=0;k<mat->nrows;++k) {
        for (m=0;m<mat->ncols;++m) {
          mat->rowval[(k*mat->ncols)+m] /= indeg[m];
        }
      }
      break;
    default:
      eprintf("Unknown matrix type: %d\n",mat->type);
      err = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto END;
  }
  dl_free(indeg);
  indeg = NULL;

  /* allocate output matrix */
  out = matrix_calloc(1);
  matrix_init(MATRIX_TYPE_DENSE_VECTOR,mat->ncols,1,0,out);

  if (times) {
    dl_stop_timer(&pre_tmr);
    dl_start_timer(&op_tmr);
  }

  /* peform pagerank */
  rank = real_alloc(mat->ncols);
  for (r=0;r<runs;++r) {
    real_set(rank,1.0/out->nrows,out->nrows);
    for (i=0;iter==0||i<iter;++i) {
      if (times) {
        dl_start_timer(&mul_tmr);
      }
      gooseberry_spmult(mat->nrows,mat->ncols,1,mat->rowptr,mat->rowind,
          mat->rowval,rank,out->rowval,NULL,0);
      if (times) {
        dl_stop_timer(&mul_tmr);
      }
      gooseberry_scale(out->nrows,out->rowval,damping);
      /* redistrubite sunk ranks */
      dist = 0;
      for (k=0;k<nsinks;++k) {
        dist += rank[sinks[k]]; 
      }
      gooseberry_add_scalar(out->nrows,out->rowval,
          ((1.0-damping)/out->nrows) + dist);
      /* only check RMSE every 10 iterations */
      if (i%10 == 0) {
        #ifndef NO_OMP
        #pragma omp parallel default(none) shared(out,rank,error) private(diff)
        { 
          error = 0.0;
          #pragma omp for schedule(static,OMP_BIG_BLOCK) \
            reduction(+:error)
          for (k=0;k<out->nrows;++k) {
            diff = rank[k] - out->rowval[k];
            error += diff*diff;
          }
        }
        #else
        error = 0.0;
        for (k=0;k<out->nrows;++k) {
          diff = rank[k] - out->rowval[k];
          error += diff*diff;
        }
        #endif
        if (error <= minerror*minerror) {
          ++i;
          /* skip recalculating RMSE */
          goto FINISH;
        }
      }
      dl_swap(rank,out->rowval);
    }
    /* check error when finished */
    #ifndef NO_OMP
    #pragma omp parallel default(none) shared(out,rank,error) private(diff)
    { 
      error = 0.0;
      #pragma omp for schedule(static,OMP_BIG_BLOCK) \
        reduction(+:error)
      for (k=0;k<out->nrows;++k) {
        diff = rank[k] - out->rowval[k];
        error += diff*diff;
      }
    }
    #else
    error = 0.0;
    for (k=0;k<out->nrows;++k) {
      diff = rank[k] - out->rowval[k];
      error += diff*diff;
    }
    #endif
    FINISH:
    error = sqrt(error);
  }

  printf("PageRank finished in %zu iterations with an RMSE of "PF_REAL_T"\n",
      i,error);

  if (times) {
    dl_stop_timer(&op_tmr);
  } 

  if (perm) {
    matrix_unpermute(out,perm,NULL);
  }

  if (times) {
    dl_start_timer(&io_tmr);
  } 

  /* save the output */
  j = __get_file_type(outfile);
  if (__is_dense(j)) {
    if (out->type != MATRIX_TYPE_DENSE_VECTOR &&
        out->type != MATRIX_TYPE_DENSE_MATRIX) {
      matrix_densify(out);
    }
    err = gooseberry_write_dense_matrix(j,outfile,out->nrows,out->ncols,
        out->rowval);
  } else {
    if (out->type != MATRIX_TYPE_SPARSE_VECTOR &&
        out->type != MATRIX_TYPE_SPARSE_MATRIX) {
      matrix_sparsify(out);
    }
    err = gooseberry_write_sparse_matrix(j,outfile,out->nrows,out->ncols,
        out->rowptr,out->rowind,out->rowval);
  }

  if (err != GOOSEBERRY_SUCCESS) {
    goto END;
  }

  if (times) {
    dl_stop_timer(&io_tmr);
    dl_print_header("Times",'#');
    printf(" I/O: %0.04lf\n",dl_poll_timer(&io_tmr));
    printf(" Preprocessing: %0.04lf\n",dl_poll_timer(&pre_tmr));
    printf(" Compute: %0.04lf\n",dl_poll_timer(&op_tmr));
    printf("   SpMV: %0.04lf\n",dl_poll_timer(&mul_tmr));
    dl_print_footer('#');
  }

  END:

  if (indeg) {
    dl_free(indeg);
  }
  if (sinks) {
    dl_free(sinks);
  }
  if (rank) {
    dl_free(rank);
  }
  if (mat) {
    matrix_free(mat);
  } 
  if (perm) {
    dl_free(perm);
  }
  if (pk) {
    dl_free(pk);
  }
  if (out) {
    matrix_free(out);
  }
  if (args) {
    dl_free(args); 
  }

  return err;
}


typedef int (*__cmdfuncptr_t)(int,char**); 
static const __cmdfuncptr_t COMMAND_FUNCS[] = {
  [COMMAND_HELP] = __help,
  [COMMAND_ANALYSIS] = __analyze,
  [COMMAND_PERMUTE] = __permute,
  [COMMAND_TRANSFORM] = __transform,
  [COMMAND_GENERATE] = __generate,
  [COMMAND_BLAS] = __blas,
  [COMMAND_CGD] = __cgd,
  [COMMAND_SGD] = __sgd,
  [COMMAND_PAGERANK] = __pagerank
};

/* don't get burned */
DL_STATIC_ASSERT(ARRAY_SIZE(COMMANDS) == ARRAY_SIZE(COMMAND_FUNCS));




/******************************************************************************
* MAIN ************************************************************************
******************************************************************************/


int main(
    int argc, 
    char ** argv)
{
  int err;
  char * cmdstr;
  size_t i;

  dl_init_rand();

  if (argc < 2) {
    eprintf("Must supply a command.\n");
    __usage(argv[0],stderr);
    return 1;
  }
  cmdstr = argv[1];

  for (i=0;i<NCOMMANDS;++i) {
    if (COMMANDS[i].str != NULL && strcmp(cmdstr,COMMANDS[i].str) == 0) {
      err = COMMAND_FUNCS[i](argc,argv);
      break;
    }
  }
  if (i == NCOMMANDS) {
    eprintf("Unrecognized command '%s'.\n",cmdstr);
    __usage(argv[0],stderr);
    return 1;
  }

  if (err == GOOSEBERRY_SUCCESS) {
    return 0;
  } else {
    eprintf("Operation failed.\n");
    return 2;
  }
}

  


#endif
