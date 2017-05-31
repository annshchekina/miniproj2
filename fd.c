#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "plotfunc.h"
#include "linalg.h"
#include "hamiltonian.h"
#include "powermethod.h"
#include "time_mes.h"

#define pi 3.1415926535897932384626433832795028L
// time measurement
#define time_test 0
const int time_test_iter = 1; 

void fd(const double mass, const double omega, const int N, const int M, const int nThreads)
{
    const double x_max = 1.0;
    const double h = x_max / N; /* h is the spacing between grid points */
	if(!time_test)
    	printf("h = %f\n", h);

    const double omega_factor = 0.5 * mass * omega * omega; // bitwise
    const double factor = -1 / (2 * mass * h * h); // bitwise
	
    /* Construct Hamiltonian matrix H */
    /* Create H_shifted as H shifted by one_norm_of_H in order to make
    the smallest eigenvalue the one with largest magnitude. */
    double* H = (double*)malloc(N*N*sizeof(double));
	double* H_shifted = (double*)malloc(N*N*sizeof(double));
	double potential_vec[N];
	const double factor_two = 2 * factor; // bitwise
    int i, j;
    for(i = 0; i < N; i++)
    	for(j = 0; j < N; j++)
		{
			const int offset = i*N+j; // keep in mind blas has different align
			double h_elem = 0;
			if (i == j)
			{
				const double potential = get_potential_value(h*i, omega_factor);
				h_elem += potential - factor_two;	
				potential_vec[i] = potential;
			}		    	
			if(i == j-1 || i == j+1 || (i == N-1 && j == 0) || (i == 0 && j == N-1))
				h_elem += factor;		
			H[offset] = h_elem; 
	        H_shifted[offset] = h_elem;
		}

	double one_norm_of_H = get_one_norm(N, H);
	if(!time_test)
	    printf("one_norm_of_H = %f\n", one_norm_of_H);
	for(i = 0; i < N; i++)
		H_shifted[i*N+i] -= one_norm_of_H; 
	
    /* Set up guess vector */
    double v_guess[N]; 
    const double N_factor = N * 0.001;
	const double N_half = N / 2; // bitwise
    for(i = 0; i < N_half + 1; i++) 
        v_guess[i] = i * 0.001;
    for(i = N_half + 1; i < N; i++)
        v_guess[i] = N_factor - i * 0.001;
    
    double eigValApprox = 0;
    double v_from_power_method[N];
    int noOfIterations = M;
    do_power_method(N,
  		  v_guess,
  		  H_shifted,
  		  &eigValApprox,
  		  v_from_power_method,
  		  noOfIterations);
    const double eigValApprox_shiftedBack = eigValApprox + one_norm_of_H;

    /* Get wavefunction squared */
	/* Normalize wavefunction */
	/* Get wavefunction_squared_sum */
	/* Get wavefunction_squared_max */
    double wavefunction_squared[N];
	double wavefunction_normalized_squared[N];
	double wavefunction_squared_sum = 0;
	double wavefunction_squared_max = 0;
	const double h_sqrt = sqrt(h);
    for(i = 0; i < N; i++)
    {
  	  const double v_pm = v_from_power_method[i];
	  const double wf_sqr = v_pm * v_pm; 
  	  wavefunction_squared[i] = wf_sqr;
	  wavefunction_squared_sum += wf_sqr; 
	  const double wf_norm = v_pm / h_sqrt;
	  const double wf_norm_sqr = wf_norm * wf_norm; 
	  wavefunction_normalized_squared[i] = wf_norm_sqr; 
      if(wf_norm_sqr > wavefunction_squared_max)
        wavefunction_squared_max = wf_norm_sqr;
    }
	
    /* Get analytical solution */
    double wavefunction_analytical_squared[N];
	double maxAbsDiff = 0;
    const double pow_factor = pow(mass*omega/pi, 0.25);
    const double mass_omega_factor = -0.5 * mass * omega; // bitwise
    for(i = 0; i < N; i++) {
      double x = h*i - 0.5;
      double psi = pow_factor * exp(mass_omega_factor * x * x); 
	  const double wf_anal_sqr = psi*psi;
      wavefunction_analytical_squared[i] = wf_anal_sqr;
      const double absDiff = fabs(wavefunction_normalized_squared[i] - wf_anal_sqr);
      if(absDiff > maxAbsDiff)
        maxAbsDiff = absDiff;
    }
    const double energy_analytical = 0.5 * omega; // bitwise
	const double maxRelDiff = maxAbsDiff / wavefunction_squared_max;
	const double energy_diff = eigValApprox_shiftedBack - energy_analytical;
		
    if(!time_test)
	{
		printf("eigValApprox_shiftedBack = %f\n", eigValApprox_shiftedBack);
		printf("energy_analytical = %f\n", energy_analytical);
		plot_function(potential_vec, N, "potential");
	    plot_function(wavefunction_squared, N, "wavefunction_squared");
	    plot_function(wavefunction_analytical_squared, N, "wavefunction_analytical_squared");
    	printf("wavefunction_squared_sum = %f\n", wavefunction_squared_sum);
    	printf("wavefunction_squared_max = %f\n", wavefunction_squared_max);
		printf("Max rel difference between finite-difference computed and analytical wavefunction_squared value: %g\n", maxRelDiff);
		printf("energy_diff = %13.9f\n", energy_diff);
	}
	
    free(H);
    free(H_shifted);
}

int main(int argc, char** argv) 
{
  if(argc != 6) 
  {
    printf("Please give 5 arguments:\n");
    printf("   mass\n");
    printf("   omega\n");
    printf("   N (number of grid points)\n");
    printf("   M (number of power-method iterations)\n");
    printf("   nThreads (number of threads to use)\n");
    return -1;
  }
  const double mass  = atof(argv[1]);
  const double omega = atof(argv[2]);
  const int N        = atoi(argv[3]);
  const int M        = atoi(argv[4]);
  const int nThreads = atoi(argv[5]);
  printf("mass     = %f\n", mass);
  printf("omega    = %f\n", omega);
  printf("N        = %d\n", N);
  printf("M        = %d\n", M);
  printf("nThreads = %d\n", nThreads);

  int i;
  double t_end_min = DBL_MAX;
  for(i = 0; i < time_test_iter; i++)
  {
	  const double t_start = get_wall_seconds();
	  fd(mass, omega, N, M, nThreads);
	  const double t_end = get_wall_seconds() - t_start;
	  if(t_end < t_end_min)
		  t_end_min = t_end;
	  printf("Iteration %d finished\n", i);
  }
  
  printf("Execution took %f sec\n", t_end_min);
  
  return 0;
}
