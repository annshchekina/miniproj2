#include "powermethod.h"
#include "linalg.h"
#include <stdio.h>
#include <math.h>

#define time_test 1

void do_power_method(int N, // add restrict
		     const double* guess,
		     const double* matrix,
		     double* resultEigVal,
		     double* resultVec,
		     int noOfIterations) 
{
  if(!time_test)
  	printf("This is the do_power_method() function, N = %d, noOfIterations = %d\n", N, noOfIterations);
  double v[N];
  int i, k;
  for(i = 0; i < N; i++)
    v[i] = guess[i];
  double eigValApprox = 0;
  for(k = 0; k < noOfIterations; k++) 
  {
    double y[N];
    do_mat_vec_mul(N, matrix, v, y);
    double dotprod = get_dot_prod(N, v, y); // why not moving it inside printing loop?
    double vv_prod = get_dot_prod(N, v, v); 
    eigValApprox = dotprod / vv_prod;
    if(!time_test && k % 2000 == 0)
      printf("k = %6d  eigValApprox = %17.10f\n", k, eigValApprox);
    for(i = 0; i < N; i++) // simd, denominator once
      v[i] = y[i] / sqrt(vv_prod);
  }
  double norm_of_v = sqrt(get_dot_prod(N, v, v));
  for(i = 0; i < N; i++)
    resultVec[i] = v[i] / norm_of_v;
  *resultEigVal = eigValApprox;
}
