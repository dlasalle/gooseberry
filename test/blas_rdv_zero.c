#include "test.h"

#define N 100




sint_t test(void) {
  sint_t i;

  real_dense_vector_t * vec = real_dv_init(N);

  real_dv_zero(vec);

  for (i=0;i<N;++i) {
    TESTEQUALS(vec->val[i],0.0,PF_REAL_T);
  }

  real_dv_free(vec); 

  return 0;
}
