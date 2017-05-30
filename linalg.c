#include <math.h>
#include "linalg.h"

void do_mat_vec_mul(const int n, const double * matrix, const double * x, double* y) {
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
  for(i = 0; i < N; i++)
    sum += a[i] * b[i];
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
