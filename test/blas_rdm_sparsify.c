#include "test.h"

#define N 100
#define M 200




sint_t test(void) {
  sint_t i,j;

  real_dense_matrix_t * mat = real_dm_init(N,M);

  for (i=0;i<N;++i) {
    for (j=0;j<M;++j) {
      if (i != j) {
        mat->val[i][j] = 0.0;
      } else {
        mat->val[i][j] = 1.0;
      }
    }
  }

  real_sparse_matrix_t * spm = real_dm_sparsify(mat);
  real_dm_free(mat);

  for (i=0;i<N;++i) {
    TESTEQUALS(spm->val[i],1.0,PF_REAL_T);
  }

  real_sm_free(spm);

  return 0;
}
