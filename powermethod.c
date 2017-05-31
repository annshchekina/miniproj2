#include "powermethod.h"
#include "linalg.h"
#include <stdio.h>
#include <math.h>
#include <pmmintrin.h>

#define time_test 0

void do_power_method(const int N, // add restrict
		     const double * __restrict guess,
		     const double * __restrict matrix,
		     double * __restrict resultEigVal,
		     double * __restrict resultVec,
		     const int noOfIterations) 
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
    const double vv_prod = get_dot_prod(N, v, v); 
    if(k % 2000 == 0)
	{
		const double dotprod = get_dot_prod(N, v, y);
		eigValApprox = dotprod / vv_prod;
		if(!time_test)
			printf("k = %6d  eigValApprox = %17.10f\n", k, eigValApprox);
	}
	if(k == noOfIterations - 1)
	{
		const double dotprod = get_dot_prod(N, v, y);
		eigValApprox = dotprod / vv_prod;
	}
	const __m128d vv_prod_sqrt_vec = _mm_set1_pd(sqrt(vv_prod));
    for(i = 0; i < N; i += 2)
	{
		const __m128d y_vec = _mm_load_pd(y + i);
		const __m128d div_vec = _mm_div_pd(y_vec, vv_prod_sqrt_vec);
		_mm_store_pd(v + i, div_vec);
	}
  }
  double norm_of_v = sqrt(get_dot_prod(N, v, v));
  for(i = 0; i < N; i++)
    resultVec[i] = v[i] / norm_of_v;
  *resultEigVal = eigValApprox;
}
