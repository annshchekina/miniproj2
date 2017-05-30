#include "hamiltonian.h"

double get_potential_value(double x, double mass, double omega) {
  /* Harmonic oscillator potential */
  double dx = x - 0.5;
  return 0.5*mass*omega*omega*dx*dx; // pass 0.5*mass*omega*omega
}

double get_kinetic_energy_term(int i, int j, int N, double h, double mass) {
  double factor = -1 / (2 * mass) / (h*h);
  /* Use 2nd order central difference scheme: (1,-2,1) */
  if(i == j-1 || (i == N-1 && j == 0))
    return 1 * factor;
  if(i == j)
    return -2 * factor;
  if(i == j+1 || (i == 0 && j == N-1))
    return 1 * factor;
  /* In all other cases, return 0 */
  return 0;
}

double get_potential_energy_term(int i, int j, int N, double h, double mass, double omega) {
  if(i != j)
    return 0;
  double x = h*i; 
  return get_potential_value(x, mass, omega);
}

double get_H_matrix_element(int i, int j, int N, double h, double mass, double omega) {
  /* Two terms: kinetic energy term and potential energy term */
  double sum = 0;
  sum += get_kinetic_energy_term(i, j, N, h, mass);
  sum += get_potential_energy_term(i, j, N, h, mass, omega);
  return sum;
}
