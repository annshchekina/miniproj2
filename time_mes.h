#include <sys/time.h>
#include <stdlib.h>

static inline double get_wall_seconds()
	{
	    struct timeval tv;
	    gettimeofday(&tv, NULL);
	    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
	    return seconds;
	};