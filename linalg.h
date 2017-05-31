void do_mat_vec_mul(const int n, const double* __restrict matrix, const double* __restrict x, double* __restrict y);
double get_dot_prod(const int N, const double* a, const double* b) __attribute__((const));
double get_one_norm(const int N, const double* A) __attribute__((const));
