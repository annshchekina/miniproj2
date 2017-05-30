#include "hamiltonian.h"

double get_potential_value(const double x, const double omega_factor) {
  /* Harmonic oscillator potential */
  const double dx = x - 0.5; 
  return omega_factor * dx * dx; 
}

double get_kinetic_energy_term(const int i, const int j, const int N, const double factor) { // implement inside i j loop
  /* Use 2nd order central difference scheme: (1,-2,1) */
  if(i == j-1 || (i == N-1 && j == 0))
    return factor;
  if(i == j)
    return -2 * factor;
  if(i == j+1 || (i == 0 && j == N-1)) // rule out some opportunities
    return factor;
  /* In all other cases, return 0 */
  return 0;
}

double get_potential_energy_term(const int i, const int j, const double h, const double omega_factor) {
  if(i != j) // i==j is a thing anyway, so handle it inside i j loop
    return 0;
  const double x = h*i; 
  return get_potential_value(x, omega_factor);
}

double get_H_matrix_element(const int i, const int j, const int N, const double h, const double factor, const double omega_factor) { // inline?
  /* Two terms: kinetic energy term and potential energy term */
  double sum = 0;
  sum += get_kinetic_energy_term(i, j, N, factor);
  sum += get_potential_energy_term(i, j, h, omega_factor);
  return sum;
}
