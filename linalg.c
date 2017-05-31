#include <math.h>
#include <pmmintrin.h>
#include "linalg.h"

#define SIMD 0

void do_mat_vec_mul(const int n, const double * __restrict matrix, const double * __restrict x, double * __restrict y) { 
  int i, j;
  for(i = 0; i < n; i++) {
    double sum = 0;
    for(j = 0; j < n; j++)
      sum += matrix[j*n+i] * x[j];
    y[i] = sum;
  }
}

double get_dot_prod(const int N, const double * a, const double * b) { // simd dot product
  double sum = 0;
  int i;
#if SIMD
  const __m128d vec_zero = _mm_setzero_pd();
  for(i = 0; i < N; i += 2)
  {
  	const __m128d a_vec = _mm_load_pd(a + i);
  	const __m128d b_vec = _mm_load_pd(b + i);
  	const __m128d mult = _mm_mul_pd(a_vec, b_vec);
  	const __m128d reduc = _mm_hadd_pd(mult, vec_zero);
  	sum += _mm_cvtsd_f64(reduc);
  }
#else
  for(i = 0; i < N; i++)		 
    sum += a[i] * b[i];
#endif
  return sum;
}

/* Note: since we assume the matrix is symmetric, the one-norm and the
   inf-norm are the same, so it does not matter if we take row sums or
   column sums. */
double get_one_norm(const int N, const double * A) {
  double maxColSum = 0;
  int i, k;
  for(i = 0; i < N; i++) {
    double currColSum = 0;
    for(k = 0; k < N; k++)
      currColSum += fabs(A[i*N+k]); // better locality with RowSum
    if(currColSum > maxColSum)
      maxColSum = currColSum;
  }
  return maxColSum;
}
