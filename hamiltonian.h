static inline double get_potential_value(const double x, const double omega_factor) 
	{
	    const double dx = x - 0.5; 
	    return omega_factor * dx * dx; 
	};
