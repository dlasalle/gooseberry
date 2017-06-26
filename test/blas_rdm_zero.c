#include "test.h"

#define N 100
#define M 200




sint_t test(void) {
  sint_t i,j;

  real_dense_matrix_t * mat = real_dm_init(N,M);

  real_dm_zero(mat);

  for (i=0;i<N;++i) {
    for (j=0;j<M;++j) {
      TESTEQUALS(mat->val[i][j],0.0,PF_REAL_T);
    }
  }

  real_dm_free(mat);

  return 0;
}
