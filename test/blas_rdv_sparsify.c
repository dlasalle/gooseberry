#include "test.h"

#define N 100
#define M 200






sint_t test(void) {
  sint_t i;

  real_dense_vector_t * vec = real_dv_init(N);

  for (i=0;i<N;++i) {
    if (i%5 != 0) {
      vec->val[i] = 0.0;
    } else {
      vec->val[i] = 1.0;
    }
  }

  real_sparse_vector_t * spv = real_dv_sparsify(vec);
  real_dv_free(vec);

  TESTEQUALS(spv->nnz,20L,PF_SSIZE_T);
  for (i=0;i<spv->nnz;++i) {
    TESTEQUALS(spv->val[i],1.0,PF_REAL_T);
  }

  real_sv_free(spv);

  return 0;
}
