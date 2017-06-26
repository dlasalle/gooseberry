/**
 * @file io.c
 * @brief I/O functions for matrices (and vectors)
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-05-01
 */




#ifndef GOOSEBERRY_IO_C
#define GOOSEBERRY_IO_C




#include "io.h"




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define IDX(i,j,ncols) ((((ind_t)i)*((ind_t)ncols))+((ind_t)j))




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef enum matrix_market_types {
  MATRIX_MARKET_NULL,
  MATRIX_MARKET_MATRIX,
  MATRIX_MARKET_VECTOR
} matrix_market_types;


typedef struct matrow_t {
  dim_t ncol;
  dim_t * ind;
  real_t * val;
  struct matrow_t * next;
} matrow_t;




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static const ind_t SVM_OFFSET = 1;
static const size_t BUFFERSIZE = 0x1000;
static const char COMMENT_CHARS[256] = {
  ['#']=1,
  ['%']=1,
  ['\'']=1,
  ['"']=1,
  ['/']=1
};
static const char * const MATRIX_MARKET_HEADER = "%%MatrixMarket";
static const char * const MATRIX_MARKET_WHITESPACE = " \t";




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static inline matrow_t * __create_row(
    dim_t const ncols)
{
  matrow_t * row;

  row = malloc(sizeof(matrow_t)+((sizeof(dim_t)+sizeof(real_t))*ncols));

  row->ncol = ncols;
  row->next = NULL;
  row->ind = (dim_t*)(row+1);
  row->val = (real_t*)(row->ind+ncols);

  return row;
} 


/* UTILITY FUNCTIONS *********************************************************/

/**
 * @brief Wrapper around the dl_file routines to handle error codes
 *
 * @param filename The file to open
 * @param mode The mode to open the file in
 * @param r_file The reference to the file pointer
 *
 * @return GOOSEBERRY_SUCCESS, or an error code
 */
static inline int __open_file(
    const char * const filename, 
    const char * const mode, 
    file_t ** const r_file)
{
  int rv;
  if ((rv = dl_open_file(filename,mode,r_file)) != DL_FILE_SUCCESS) {
    switch (rv) {
      case DL_FILE_BAD_PARAMETERS:
      case DL_FILE_PATH_PARSE_FAILURE:
        eprintf("Bad filename '%s'\n",filename);
        rv = GOOSEBERRY_ERROR_INVALIDINPUT;
        break;
      case DL_FILE_PATH_BAD:
        eprintf("File not found '%s'\n",filename);
        rv = GOOSEBERRY_ERROR_FILENOTFOUND;
        break;
      case DL_FILE_PATH_ACCESS_DENIED:
      case DL_FILE_READ_ACCESS_DENIED:
      case DL_FILE_WRITE_ACCESS_DENIED:
        eprintf("Permission denied '%s'\n",filename);
        rv = GOOSEBERRY_ERROR_PERMISSIONDENIED;
        break;
      default:
        eprintf("Unknown failure: %d opening '%s'\n",rv,filename);
        rv = GOOSEBERRY_ERROR_UNKNOWN;
        break;
    }
  } else {
    rv = GOOSEBERRY_SUCCESS;
  }
  return rv;
}


/* READ MATRIX FUNCTIONS *****************************************************/

static int __read_grid(
    const char * const file,
    dim_t * const r_nrows,
    dim_t * const r_ncols,
    real_t ** const r_rowval)
{
  int rv;
  ssize_t ll;
  size_t bufsize;
  dim_t ncols, nrows, j;
  ind_t k;
  real_t f;
  file_t * fin = NULL;
  char * line = NULL, * eptr, * sptr;
  real_t * rowval = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(file,"r",&fin)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  nrows = 0;
  ncols = 0;

  /* read the whole matrix to get the size */
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (COMMENT_CHARS[(unsigned int)line[0]]) {
      /* skip comments */
      continue;
    }
    j = 0;
    sptr = line;
    f = strtod(sptr,&eptr);
    while (eptr != sptr) {
      ++j;
      sptr = eptr;
      f = strtod(sptr,&eptr);
    }
    if (ncols == 0) {
      ncols = j;
    } else {
      if (j == 0) {
        eprintf("Invalid line in grid matrix file: "PF_DIM_T"\n",nrows+1);
        rv = GOOSEBERRY_ERROR_INVALIDINPUT;
        goto ERROR;
      } else if (j != ncols) {
        eprintf("Different number of columns ("PF_DIM_T") for line "PF_DIM_T
            " when previous had "PF_DIM_T" columns.\n",j,nrows,ncols);
        rv = GOOSEBERRY_ERROR_INVALIDINPUT;
        goto ERROR;
      }
    }
    ++nrows; 
  }

  dl_reset_file(fin);

  *r_rowval = rowval = real_alloc(nrows*ncols);

  k = 0;
  /* read the whole matrix to get the size */
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (COMMENT_CHARS[(unsigned int)line[0]]) {
      /* skip comments */
      continue;
    }
    sptr = line;
    f = (real_t)strtod(sptr,&eptr);
    while (eptr != sptr) {
      rowval[k++] = f;
      sptr = eptr;
      f = (real_t)strtod(sptr,&eptr);
    }
  }

  dl_free(line);
  dl_close_file(fin);

  *r_nrows = nrows;
  *r_ncols = ncols;

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (line) {
    dl_free(line);
  }
  if (rowval) {
    dl_free(rowval);
  }

  return rv;
}


static int __read_graph(
    const char * const file,
    dim_t * const nrows,
    dim_t * const ncols,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  int rv, do_wgt;
  ssize_t ll;
  size_t bufsize;
  ind_t j;
  dim_t i, n;
  real_t f;
  file_t * fin = NULL;
  char * line = NULL, * eptr, * sptr;

  ind_t * rowptr = NULL;
  dim_t * rowind = NULL;
  real_t * rowval = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(file,"r",&fin)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  n = 0;
  j = 0;
  do_wgt = 0;
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* badness */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    break;
  }
  sptr = line;
  n = (dim_t)strtoull(sptr,&eptr,10);
  sptr = eptr;
  j = (ind_t)strtoull(sptr,&eptr,10) * 2;
  sptr = eptr;
  do_wgt = (int)strtoull(sptr,&eptr,10);
  if (sptr == eptr) {
    do_wgt = 0;
  }

  if (n == 0 && j == 0) {
    eprintf("Unable to read graph file or empty graph specified\n");
    rv = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto ERROR;
  }

  *r_rowptr = rowptr = ind_alloc(n+1); 
  *r_rowind = rowind = dim_alloc(j); 
  *r_rowval = rowval = real_alloc(j); 

  n = 0;
  j = 0;
  rowptr[0] = 0;
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* assume island vertex */
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = (dim_t)strtoull(sptr,&eptr,10);

    if (do_wgt) {
      sptr = eptr;
      f = (real_t)strtod(sptr,&eptr);
    }
    while (eptr != sptr) { 
      rowind[j] = i-1;
      if (do_wgt) {
        rowval[j] = f;
      } else {
        rowval[j] = 1.0;
      }
      ++j;
      sptr = eptr;
      i = (dim_t)strtoull(sptr,&eptr,10);
      if (do_wgt) {
        sptr = eptr;
        f = (real_t)strtod(sptr,&eptr);
      }
    }
    rowptr[++n] = j;
  }

  dl_free(line);
  dl_close_file(fin);

  *nrows = n;
  *ncols = n;

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (line) {
    dl_free(line);
  }
  if (rowptr) {
    dl_free(rowptr);
  }
  if (rowind) {
    dl_free(rowind);
  }
  if (rowval) {
    dl_free(rowval);
  }

  return rv;
}


static int __read_csr(
    const char * const file,
    dim_t * const nrows,
    dim_t * const ncols,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  int rv;
  ssize_t ll;
  size_t bufsize;
  ind_t j;
  dim_t i, maxcol, mincol, n, rowcol;
  real_t f;
  file_t * fin = NULL;
  char * line = NULL, * eptr, * sptr;

  ind_t * rowptr = NULL;
  dim_t * rowind = NULL;
  real_t * rowval = NULL;

  matrow_t * oldrow, * row;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(file,"r",&fin)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  oldrow = row = NULL;
  n = 0;
  j = 0;
  mincol = 1;
  maxcol = 0;
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* a blank line we'll assume means an empty row */
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }

    /* count the row */
    rowcol = (dim_t)(dl_get_ne_str(line) / 2);

    if (row == NULL) {
      oldrow = row = __create_row(rowcol);
    } else {
      row->next = __create_row(rowcol);
      row = row->next;
    }

    /* save the row */
    rowcol = 0;
    sptr = line;
    i = (dim_t)strtoull(sptr,&eptr,10);
    sptr = eptr;
    f = strtod(sptr,&eptr);
    while (eptr != sptr) {
      if (i > maxcol) {
        maxcol = i;
      }
      if (i < mincol) {
        mincol = i;
      }

      row->ind[rowcol] = i;
      row->val[rowcol] = f;
      ++rowcol;

      sptr = eptr;
      i = (dim_t)strtoull(sptr,&eptr,10);
      sptr = eptr;
      f = strtod(sptr,&eptr);
    }
    j += rowcol;

    ++n;
  }

  dl_reset_file(fin);

  *r_rowptr = rowptr = ind_alloc(n+1); 
  *r_rowind = rowind = dim_alloc(j); 
  *r_rowval = rowval = real_alloc(j); 

  n = 0;
  j = 0;
  rowptr[0] = 0;
  row = oldrow;
  while(row != NULL) {
    for (i=0;i<row->ncol;++i) {
      /* mincol is at most 1 */
      rowind[j] = row->ind[i]-mincol;
      rowval[j] = row->val[i];
      ++j;
    }
    oldrow = row;
    row = oldrow->next;
    dl_free(oldrow);
    rowptr[++n] = j;
  }

  dl_free(line);
  dl_close_file(fin);

  *nrows = n;
  *ncols = maxcol-mincol+1;

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (line) {
    dl_free(line);
  }
  if (rowptr) {
    dl_free(rowptr);
  }
  if (rowind) {
    dl_free(rowind);
  }
  if (rowval) {
    dl_free(rowval);
  }

  return rv;
}


static int __read_clu(
    const char * const file,
    dim_t * const r_nrows,
    dim_t * const r_ncols,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  int rv;
  ssize_t ll;
  size_t bufsize;
  ind_t j;
  dim_t i, maxcol, mincol, n;
  real_t f;
  file_t * fin = NULL;
  char * line = NULL, * eptr, * sptr;

  ind_t * rowptr = NULL;
  dim_t * rowind = NULL;
  real_t * rowval = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(file,"r",&fin)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  n = 0;
  j = 0;
  mincol = 1;
  maxcol = 0;
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* badness */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    break;
  }
  mincol = 1; /* specified by cluto manual */
  sscanf(line,PF_DIM_T" "PF_DIM_T" "PF_IND_T,&n,&maxcol,&j);

  *r_rowptr = rowptr = ind_alloc(n+1); 
  *r_rowind = rowind = dim_alloc(j); 
  *r_rowval = rowval = real_alloc(j); 

  n = 0;
  j = 0;
  rowptr[0] = 0;
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* badness */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = (dim_t)strtoull(sptr,&eptr,10);
    sptr = eptr;
    f = (real_t)strtod(sptr,&eptr);
    while (eptr != sptr) { 
      /* mincol is at most 1 */
      rowind[j] = i-mincol;
      rowval[j] = f;
      ++j;
      sptr = eptr;
      i = (dim_t)strtoull(sptr,&eptr,10);
      sptr = eptr;
      f = (real_t)strtod(sptr,&eptr);
    }
    rowptr[++n] = j;
  }

  dl_free(line);
  dl_close_file(fin);

  if (r_nrows) {
    *r_nrows = n;
  }
  if (r_ncols) {
    *r_ncols = maxcol-mincol+1;
  }

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (line) {
    dl_free(line);
  }
  if (rowptr) {
    dl_free(rowptr);
  }
  if (rowind) {
    dl_free(rowind);
  }
  if (rowval) {
    dl_free(rowval);
  }

  return rv;
}


static int __read_svm(
    const char * const file,
    dim_t * const r_nrows,
    dim_t * const r_ncols,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  int rv;
  ssize_t ll;
  size_t bufsize;
  ind_t j;
  dim_t i, maxcol,n;
  real_t f;
  file_t * fin = NULL;
  char * line = NULL, * eptr, * sptr;

  ind_t * rowptr = NULL;
  dim_t * rowind = NULL;
  real_t * rowval = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(file,"r",&fin)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  n = 0;
  j = 0;
  maxcol = 0;
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* badness */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    /* read leading number */
    sptr = line;
    i = (dim_t)strtoll(sptr,&eptr,10);

    sptr = eptr;
    i = (dim_t)strtoull(sptr,&eptr,10);
    /* skip colon */
    sptr = eptr+1;
    f = (real_t)strtod(sptr,&eptr);
    while (eptr != sptr) {
      if (i > maxcol) {
        maxcol = i;
      }
      ++j;
      sptr = eptr;
      i = (dim_t)strtoull(sptr,&eptr,10);
      sptr = eptr+1;
      f = (real_t)strtod(sptr,&eptr);
    }
    ++n;
  }

  dl_reset_file(fin);

  *r_rowptr = rowptr = ind_alloc(n+1); 
  *r_rowind = rowind = dim_alloc(j); 
  *r_rowval = rowval = real_alloc(j); 

  n = 0;
  j = 0;
  rowptr[0] = 0;
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* badness */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }

    /* read leading number */
    sptr = line;
    i = (dim_t)strtoll(sptr,&eptr,10);

    sptr = eptr;
    i = (dim_t)strtoull(sptr,&eptr,10);
    /* skip colon */
    sptr = eptr+1;
    f = (real_t)strtod(sptr,&eptr);
    while (eptr != sptr) {
      rowind[j] = i-SVM_OFFSET;
      rowval[j] = f;
      ++j;
      sptr = eptr;
      i = (dim_t)strtoull(sptr,&eptr,10);
      sptr = eptr+1;
      f = (real_t)strtod(sptr,&eptr);
    }
    rowptr[++n] = j;
  }

  dl_free(line);
  dl_close_file(fin);

  if (r_nrows) {
    *r_nrows = n;
  } 
  if (r_ncols) {
    *r_ncols = maxcol;
  }

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (line) {
    dl_free(line);
  }
  if (rowptr) {
    dl_free(rowptr);
  }
  if (rowind) {
    dl_free(rowind);
  }
  if (rowval) {
    dl_free(rowval);
  }

  return rv;
}


static int __read_point(
    const char * const filename,
    dim_t * const r_nrows,
    dim_t * const r_ncols,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  file_t * file;
  int rv;
  dim_t i, j, ncols, nrows, mincol, minrow;
  ind_t nnz;
  ssize_t ll;
  size_t bufsize;
  char * line = NULL,*eptr, * sptr;
  real_t w;

  ind_t * rowptr = NULL;
  dim_t * rowind = NULL;
  real_t * rowval = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(filename,"r",&file)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }
  line = char_alloc(bufsize);

  minrow = 1;
  mincol = 1;
  nrows = 0;
  ncols = 0;

  /* first pass to count edges and vertices */
  nnz = 0;
  while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* ignore blank lines */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = (dim_t)strtoull(sptr,&eptr,10);
    sptr = eptr;
    j = (dim_t)strtoull(sptr,&eptr,10);
    sptr = eptr;
    w = (real_t)strtod(sptr,&eptr);
    if (eptr == sptr) {
      eprintf("Error in point list at line "PF_IND_T"\n",nnz);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    }
    if (i > nrows) {
      nrows = i;
    } else {
      minrow = i;
    }
    if (j > ncols) {
      ncols = j;
    } else {
      mincol = j;
    }
    ++nnz;
  }

  /* adjust for 1 based indexing */
  nrows += 1-minrow;
  ncols += 1-mincol;

  rowptr = ind_calloc(nrows+1);
  rowind = dim_alloc(nnz); 
  rowval = real_alloc(nnz);

  dl_reset_file(file);

  /* second pass to count row widths */
  nnz = 0;
  while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* ignore blank lines */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = (dim_t)strtoull(sptr,&eptr,10) - minrow;
    sptr = eptr;
    j = (dim_t)strtoull(sptr,&eptr,10) - mincol;
    sptr = eptr;
    w = (real_t)strtod(sptr,&eptr);

    /* increment rowptr */
    ++rowptr[i+1]; 
  }

  ind_prefixsum_exc(rowptr+1,nrows);
  dl_reset_file(file);

  /* third pass to fill arrays */
  while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* ignore blank lines */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = (dim_t)strtoull(sptr,&eptr,10) - minrow;
    sptr = eptr;
    j = (dim_t)strtoull(sptr,&eptr,10) - mincol;
    sptr = eptr;
    w = (real_t)strtod(sptr,&eptr);

    /* fill arrays */
    nnz = rowptr[i+1]++; 
    rowind[nnz] = j;
    rowval[nnz] = w;
  }

  dl_free(line);
  dl_close_file(file);

  if (r_nrows) {
    *r_nrows = nrows;
  }
  if (r_ncols) {
    *r_ncols = ncols;
  }
  if (r_rowptr) {
    *r_rowptr = rowptr;
  } else {
    dl_free(rowptr);
  }
  if (r_rowind) {
    *r_rowind = rowind;
  } else if (rowind) {
    dl_free(rowind);
  }
  if (r_rowval) {
    *r_rowval = rowval;
  } else if (rowval) {
    dl_free(rowval);
  }

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (line) {
    dl_free(line);
  }
  if (rowptr) {
    dl_free(rowptr);
  }
  if (rowind) {
    dl_free(rowind);
  }
  if (rowval) {
    dl_free(rowval);
  }

  return rv;
}


static int __read_matrixmarket(
    const char * const filename,
    dim_t * const r_nrows,
    dim_t * const r_ncols,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    real_t ** const r_rowval)
{
  file_t * file;
  int rv, type, sym, wgts;
  dim_t i, j, ncols, nrows;
  ind_t nnz;
  ssize_t ll;
  size_t bufsize;
  char * line = NULL,*eptr, * sptr;
  real_t w;

  ind_t * rowptr = NULL;
  dim_t * rowind = NULL;
  real_t * rowval = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(filename,"r",&file)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }
  line = char_alloc(bufsize);

  /* read the type code line */
  if ((ll = dl_get_next_line(file,&line,&bufsize)) < 1) {
    eprintf("Missing header for matrix market file '%s'\n",filename);
    rv = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto ERROR;
  }
  if (sscanf(line,PF_DIM_T" "PF_DIM_T" "PF_IND_T,&nrows,&ncols,&nnz) == 3) {
    /* missing header -- assume sparse real unsymmetric */
    sym = 0;
  } else {
    /* figure out what this is */
    sptr = strtok(line,MATRIX_MARKET_WHITESPACE);
    if (strcmp(sptr,MATRIX_MARKET_HEADER) != 0) {
      eprintf("Bad header '%s', for matrix market file '%s'\n",sptr,filename);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    }
    sptr = strtok(NULL,MATRIX_MARKET_WHITESPACE);
    type = MATRIX_MARKET_NULL;
    if (strcmp(sptr,"matrix") == 0) {
      type = MATRIX_MARKET_MATRIX; 
    } else if (strcmp(sptr,"vector") == 0) {
      type = MATRIX_MARKET_VECTOR;
    } else {
      eprintf("Unknown type '%s' in matrix market header\n",sptr);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    }

    if (type != MATRIX_MARKET_MATRIX) {
      eprintf("Only matrices are supported currently\n");
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    }

    sptr = strtok(NULL,MATRIX_MARKET_WHITESPACE);
    if (strcmp(sptr,"coordinate") == 0) {
      /* keep going */
    } else if (strcmp(sptr,"array") == 0) {
      eprintf("Not a sparse matrix '%s'\n",filename);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    } else {
      eprintf("Unsupported matrix format '%s'\n",sptr);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    }
    sptr = strtok(NULL,MATRIX_MARKET_WHITESPACE);
    if (strcmp(sptr,"real") != 0) {
      eprintf("Unsupported matrix format '%s'\n",sptr);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    }
    sptr = strtok(NULL,MATRIX_MARKET_WHITESPACE);
    if (strcmp(sptr,"general") == 0) {
      sym = 0; 
    } else if (strcmp(sptr,"symmetric") == 0) {
      sym = 1;
    } else {
      eprintf("Unsupported matrix format '%s'\n",sptr);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    }

    /* find header line */
    while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
      if (ll == 0) {
        /* ignore blank lines */
        continue;
      } else {
        if (COMMENT_CHARS[(unsigned int)line[0]]) {
          /* skip comments */
          continue;
        }
      }
      break;
    }

    /* read dimensions */
    sptr = line;
    nrows = (dim_t)strtoull(sptr,&eptr,10);
    sptr = eptr;
    ncols = (dim_t)strtoull(sptr,&eptr,10);
    sptr = eptr;
    nnz = (ind_t)strtoull(sptr,&eptr,10);
    if (sptr == eptr) {
      eprintf("Failed to read nrows, ncols, and nnz from header line: " \
          "'%s'\n",line);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    }

    if (sym) {
      nnz *= 2;
    }
  }

  rowptr = ind_calloc(nrows+1);
  rowind = dim_alloc(nnz); 
  rowval = real_alloc(nnz);

  dl_mark_file(file);

  /* read in coordinates */
  wgts = 0;
  nnz = 0;
  while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* ignore blank lines */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = (dim_t)strtoull(sptr,&eptr,10) - 1;
    if (eptr == sptr) {
      eprintf("Bad 'i' in point list at line "PF_IND_T"\n",nnz);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    }

    sptr = eptr;
    j = (dim_t)strtoull(sptr,&eptr,10) - 1;
    if (eptr == sptr) {
      eprintf("Bad 'j' in point list at line "PF_IND_T"\n",nnz);
      rv = GOOSEBERRY_ERROR_INVALIDINPUT;
      goto ERROR;
    }

    sptr = eptr;
    w = (real_t)strtod(sptr,&eptr);
    if (eptr != sptr) {
      if (nnz == 0) {
        /* there are values here */
        wgts = 1;
      } else if (wgts == 0) {
        /* there are some values here */
        eprintf("Found first value on line "PF_IND_T", but not on preceding " \
            "lines\n",nnz);
        rv = GOOSEBERRY_ERROR_INVALIDINPUT;
        goto ERROR;
      }
    } else { 
      if (wgts) {
        eprintf("Missing value at line "PF_IND_T", but previous lines had " \
            "values\n",nnz);
        rv = GOOSEBERRY_ERROR_INVALIDINPUT;
        goto ERROR;
      }
    }

    /* increment rowptr */
    if (sym) {
      ++rowptr[j+1];
      ++nnz;
    }
    ++rowptr[i+1]; 
    ++nnz;
  }

  ind_prefixsum_exc(rowptr+1,nrows);
  dl_restore_file(file);

  /* second pass to fill arrays */
  w = 1;
  while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* ignore blank lines */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = (dim_t)strtoull(sptr,&eptr,10) - 1;
    sptr = eptr;
    j = (dim_t)strtoull(sptr,&eptr,10) - 1;
    if (wgts) {
      sptr = eptr;
      w = (real_t)strtod(sptr,&eptr);
    } 

    /* fill arrays */
    if (sym) {
      nnz = rowptr[j+1]++; 
      rowind[nnz] = i;
      rowval[nnz] = w;
    }
    nnz = rowptr[i+1]++; 
    rowind[nnz] = j;
    rowval[nnz] = w;
  }

  dl_free(line);
  dl_close_file(file);

  if (r_nrows) {
    *r_nrows = nrows;
  }
  if (r_ncols) {
    *r_ncols = ncols;
  }
  if (r_rowptr) {
    *r_rowptr = rowptr;
  } else {
    dl_free(rowptr);
  }
  if (r_rowind) {
    *r_rowind = rowind;
  } else if (rowind) {
    dl_free(rowind);
  }
  if (r_rowval) {
    *r_rowval = rowval;
  } else if (rowval) {
    dl_free(rowval);
  }

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (line) {
    dl_free(line);
  }
  if (rowptr) {
    dl_free(rowptr);
  }
  if (rowind) {
    dl_free(rowind);
  }
  if (rowval) {
    dl_free(rowval);
  }

  return rv;
}


/* WRITE MATRIX FUNCTIONS ****************************************************/

static int __write_grid(
    const char * const file,
    const dim_t nrows, 
    const dim_t ncols,
    const real_t * rowval) 
{
  int rv;
  dim_t i,j;
  file_t * fout = NULL;

  if ((rv = __open_file(file,"w",&fout)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  for (i=0;i<nrows;++i) {
    dl_fprintf(fout,PF_MINREAL_T,rowval[IDX(i,0,ncols)]);
    for (j=1;j<ncols;++j) {
      dl_fprintf(fout," "PF_MINREAL_T,rowval[IDX(i,j,ncols)]);
    }
    dl_fprintf(fout,"\n");
  }

  dl_close_file(fout);
  fout = NULL;

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (fout) {
    dl_close_file(fout);
  }

  return rv;
}



static int __write_csr(
    const char * const file,
    const dim_t nrows,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval)
{
  int rv;
  dim_t i;
  ind_t j;

  file_t * fout = NULL;

  if ((rv = __open_file(file,"w",&fout)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  for (i=0;i<nrows;++i) {
    if (rowptr[i] < rowptr[i+1]) {
      dl_fprintf(fout,PF_DIM_T" "PF_MINREAL_T,rowind[rowptr[i]]+1,
          rowval[rowptr[i]]);
      for (j=rowptr[i]+1;j<rowptr[i+1];++j) {
        dl_fprintf(fout," "PF_DIM_T" "PF_MINREAL_T,rowind[j]+1,rowval[j]);
      }
    }
    dl_fprintf(fout,"\n");
  }

  dl_close_file(fout);
  fout = NULL;

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (fout) {
    dl_close_file(fout);
  }

  return rv;

}


static int __write_graph(
    const char * const file,
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval)
{
  int rv, ewgts;
  dim_t i;
  ind_t j;
  file_t * fout = NULL;

  ind_t nnz = rowptr[nrows];

  if (nrows != ncols) {
    eprintf("Graphs must have an equal number of rows and columns: "PF_DIM_T \
        "x"PF_DIM_T"\n",nrows,ncols);
    rv = GOOSEBERRY_ERROR_INVALIDINPUT;
    goto ERROR;
  }

  if ((rv = __open_file(file,"w",&fout)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  /* filter loops */
  for (i=0;i<nrows;++i) {
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      if (i == rowind[j]) {
        --nnz;
        break;
      }
    }
  }

  /* see if we should write edge weights */
  ewgts = 0;
  if (rowval) {
    for (j=1;j<nnz;++j) {
      if (rowval[j] != rowval[j-1]) {
        ewgts = 1;
        break;
      }
    }
  }

  if (rowval && ewgts) {
    dl_fprintf(fout,PF_DIM_T" "PF_IND_T" 01\n",ncols,nnz/2);

    for (i=0;i<nrows;++i) {
      for (j=rowptr[i];j<rowptr[i+1];++j) {
        if (i != rowind[j]) {
          if (rowval[j] == (real_t)((long long)(rowval[j]))) {
            dl_fprintf(fout," "PF_DIM_T" %lld",rowind[j]+1,(long long)rowval[j]);
          } else {
            dl_fprintf(fout," "PF_DIM_T" "PF_MINREAL_T,rowind[j]+1,rowval[j]);
          }
        }
      }
      dl_fprintf(fout,"\n");
    }
  } else {
    dl_fprintf(fout,PF_DIM_T" "PF_IND_T"\n",ncols,nnz/2);

    for (i=0;i<nrows;++i) {
      for (j=rowptr[i];j<rowptr[i+1];++j) {
        if (i != rowind[j]) {
          dl_fprintf(fout," "PF_DIM_T,rowind[j]+1);
        }
      }
      dl_fprintf(fout,"\n");
    }
  }

  dl_close_file(fout);
  fout = NULL;

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (fout) {
    dl_close_file(fout);
  }

  return rv;
}


static int __write_matrixmarket(
    const char * const file,
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval)
{
  int rv;
  dim_t i;
  ind_t j;
  file_t * fout = NULL;

  ind_t nnz = rowptr[nrows];

  if ((rv = __open_file(file,"w",&fout)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  dl_fprintf(fout,"%%%%MatrixMarket matrix coordinate real general\n");
  dl_fprintf(fout,"%%\n");
  dl_fprintf(fout,"%% Written by Gooseberry version %d.%d.%d\n", \
      GOOSEBERRY_VER_MAJOR,GOOSEBERRY_VER_MINOR,GOOSEBERRY_VER_SUBMINOR);
  dl_fprintf(fout,"%%   http://cs.umn.edu/~lasalle/gooseberry\n");
  dl_fprintf(fout,"%%\n");

  /* write the header */
  dl_fprintf(fout,PF_DIM_T" "PF_DIM_T" "PF_IND_T"\n",nrows,ncols,nnz);

  /* write the matrix loops */
  for (i=0;i<nrows;++i) {
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      dl_fprintf(fout,PF_DIM_T"\t"PF_DIM_T"\t"PF_REAL_T"\n",i+1,rowind[j]+1, \
          rowval[j]);
    }
  }

  dl_close_file(fout);
  fout = NULL;

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (fout) {
    dl_close_file(fout);
  }

  return rv;
}


static int __write_point(
    const char * const file,
    const dim_t nrows,
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval)
{
  int rv;
  dim_t i;
  ind_t j;
  file_t * fout = NULL;

  if ((rv = __open_file(file,"w",&fout)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  /* write the matrix loops */
  for (i=0;i<nrows;++i) {
    for (j=rowptr[i];j<rowptr[i+1];++j) {
      dl_fprintf(fout,PF_DIM_T"\t"PF_DIM_T"\t"PF_REAL_T"\n",i+1,rowind[j]+1, \
          rowval[j]);
    }
  }

  dl_close_file(fout);
  fout = NULL;

  return GOOSEBERRY_SUCCESS;

  ERROR:

  if (fout) {
    dl_close_file(fout);
  }

  return rv;
}


/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int read_sparse_matrix(
    int type, 
    const char * const file,
    dim_t * nrows, 
    dim_t * ncols,
    ind_t ** rowptr,
    dim_t ** rowind,
    real_t ** rowval)
{
  switch(type) {
    case GOOSEBERRY_FORMAT_CSR:
      return __read_csr(file,nrows,ncols,rowptr,rowind,rowval);
    case GOOSEBERRY_FORMAT_GRAPH:
      return __read_graph(file,nrows,ncols,rowptr,rowind,rowval);
    case GOOSEBERRY_FORMAT_SVM:
      return __read_svm(file,nrows,ncols,rowptr,rowind,rowval);
    case GOOSEBERRY_FORMAT_CLU:
      return __read_clu(file,nrows,ncols,rowptr,rowind,rowval);
    case GOOSEBERRY_FORMAT_POINT:
      return __read_point(file,nrows,ncols,rowptr,rowind,rowval);
    case GOOSEBERRY_FORMAT_MATRIXMARKET:
      return __read_matrixmarket(file,nrows,ncols,rowptr,rowind,rowval);
    default:
      eprintf("Unrecognized sparse matrix format '%d'\n",type);
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
}


int read_dense_matrix(
    int type, 
    const char * const file,
    dim_t * nrows, 
    dim_t * ncols,
    real_t ** rowval)
{
  switch (type) {
    case GOOSEBERRY_FORMAT_GRID:
      return __read_grid(file,nrows,ncols,rowval);
    default:
      eprintf("Unrecognized dense matrix format '%d'\n",type);
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
}


int write_dense_matrix(
    const int type, 
    const char * const file,
    const dim_t nrows, 
    const dim_t ncols,
    const real_t * rowval)
{
  switch(type) {
    case GOOSEBERRY_FORMAT_GRID:
      return __write_grid(file,nrows,ncols,rowval);
    default:
      eprintf("Unrecognized dense matrix format '%d'\n",type);
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
}


int write_sparse_matrix(
    const int type, 
    const char * file,
    const dim_t nrows, 
    const dim_t ncols,
    const ind_t * const rowptr,
    const dim_t * const rowind,
    const real_t * const rowval)
{
  switch(type) {
    case GOOSEBERRY_FORMAT_CSR:
      return __write_csr(file,nrows,rowptr,rowind,rowval);
    case GOOSEBERRY_FORMAT_GRAPH:
      return __write_graph(file,nrows,ncols,rowptr,rowind,rowval);
    case GOOSEBERRY_FORMAT_MATRIXMARKET:
      return __write_matrixmarket(file,nrows,ncols,rowptr,rowind,rowval);
    case GOOSEBERRY_FORMAT_POINT:
      return __write_point(file,nrows,ncols,rowptr,rowind,rowval);
    default:
      eprintf("Unrecognized sparse matrix format '%d'\n",type);
      return GOOSEBERRY_ERROR_INVALIDINPUT;
  }
}


int read_labels(
    const char * const filename,
    dim_t * const r_nrows,
    dim_t ** const r_labels)
{
  file_t * file = NULL;
  dim_t nl, rowcap;
  dim_t l;
  int rv;
  ssize_t ll;
  size_t bufsize;
  char * line = NULL;
  dim_t * labels = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(filename,"r",&file)) != GOOSEBERRY_SUCCESS) {
    goto ERROR;
  }

  line = char_alloc(bufsize);
  ll = dl_get_next_line(file,&line,&bufsize);
  /* skip comments */
  while (ll > 0 && COMMENT_CHARS[(unsigned int)line[0]]) {
    ll = dl_get_next_line(file,&line,&bufsize);
  }
  if (!r_nrows || *r_nrows == 0) {
    /* determine the size of the file */
    nl = 0;
    while (ll > 0 && sscanf(line,PF_DIM_T,&l) == 1) {
      ++nl; 
      ll = dl_get_next_line(file,&line,&bufsize);
    }
    dl_reset_file(file);
    ll = dl_get_next_line(file,&line,&bufsize);
    /* skip comments */
    while (ll > 0 && COMMENT_CHARS[(unsigned int)line[0]]) {
      ll = dl_get_next_line(file,&line,&bufsize);
    }
  } else {
    /* the size of the file was specified */
    nl = *r_nrows;
  }

  labels = dim_alloc(nl);
  rowcap = nl;
  nl = 0;
  /* read in the file */
  while (ll > 0 && sscanf(line,PF_DIM_T,&l) == 1) {
    if (nl == rowcap) {
      wprintf("Found excess of rows in permutation file\n");
      break;
    }
    labels[nl++] = l;
    ll = dl_get_next_line(file,&line,&bufsize);
  }

  dl_close_file(file);
  file = NULL;

  if (r_nrows) {
    *r_nrows = nl;
  }

  dl_free(line);

  if (r_labels) {
    *r_labels = labels;
  } else {
    dl_free(labels);
  }

  return GOOSEBERRY_SUCCESS; 

  ERROR:

  if (file) {
    dl_close_file(file);
  }
  if (line) {
    dl_free(line);
  } 
  if (labels) {
    dl_free(labels);
  }

  return rv;
}


#endif
