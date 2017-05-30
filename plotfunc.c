#include <stdio.h>
#include "plotfunc.h"

void plot_function(const double* vec, const int N, const char* name) {
  printf("plot_function for '%s':\n", name);
  const int nPlotPts = 120;
  if(N < nPlotPts) {
    printf("Cannot plot: N is too small.\n");
    return;
  }
  double minValue = vec[0];
  double maxValue = vec[0]; 
  int i, k;
  for(i = 0; i < N; i++) {
    if(vec[i] < minValue)
      minValue = vec[i];
    if(vec[i] > maxValue)
      maxValue = vec[i];
  }
  minValue -= (maxValue-minValue)/10;
  maxValue += (maxValue-minValue)/10;
  const int ny = 15;
  char plotMatrix[nPlotPts*ny];
  const double stepLen = (double)N / nPlotPts;
  for(i = 0; i < nPlotPts; i++) {
    const double desiredIndex = i*stepLen;
    int idx = (int)desiredIndex;
    if(idx < 0)
      idx = 0;
    if(idx > N-1)
      idx = N-1;
    const double value = vec[idx];
    int y = (int)(ny*(value-minValue) / (maxValue-minValue));
    if(y < 0)
      y = 0;
    if(y > ny-1)
      y = ny-1;
    for(k = 0; k < ny; k++) 
	{
      char c = ' ';
      if(y >= k)
		  c = 'X';
      plotMatrix[i*ny+k] = c;
    }
  }
  for(k = 0; k < ny; k++) 
  {
    const int y = ny-1-k;
    for(i = 0; i < nPlotPts; i++)
      printf("%c", plotMatrix[i*ny+y]);
    printf("\n");
  }
}
