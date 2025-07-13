#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../mfdfa.h"	// common stuff
#include "mfdfa3d.h"	// only for 3D MFDFA version, header file in the current directory

#define XDIM	790		// adjust these to your data
#define YDIM	790
#define ZDIM	790

double *data;

int main() {
  FILE *h;
  char fname[256], fbase[] = "A1001_790_790_790";
  DFA_CONFIG cfg = { 0 };										// configuration structure
  double H[MAXQ], H2[MAXQ], tau[MAXQ], f[MAXQ], alpha[MAXQ];	// results go here
  int i, w, m, n, total, nq, integrate;
  double q, qmin, qmax, dq;
  time_t start, end;

  start = time(0);	// timer start
  qmin = -4.000001; qmax = 4.0; dq = 0.1;	// skip q=0

  data = malloc(XDIM * YDIM * ZDIM * sizeof(double));
  sprintf(fname, "data/%s.bin", fbase);
  h = fopen(fname, "rb+");
  total=fread(data, sizeof(double), XDIM*YDIM*ZDIM, h);	// read the data
  fclose(h);
  printf("%d data read from %s\n",total,fname);

  cfg.minbox = 4;									// fill in the cofiguration structure
  cfg.maxbox = min(min(XDIM,YDIM),ZDIM) / 4;
  cfg.boxratio = pow(2.0, 1.0 / 8.0);				// standard choice for multiplicative factor
  cfg.goback = 1;									// if non-overlapping windows, go backwards as well
  //cfg.sw = 1;	// sliding window - more statistics, but MUCH, MUCH, MUCH slower, be careful with this option...

  rscale(&cfg);		// prepare scale
  mfdfa_3d(&cfg, ZDIM, YDIM, XDIM, data, qmin, qmax, dq, eps, H, H2);	// call 3D MFDFA
  free(data);

  end = time(0);	// timer end
  printf("\nElapsed %d seconds\n", end-start);

  // print out the results to screen and file
  sprintf(fname, "data/%s_mfdfa%s.txt", fbase, cfg.sw ? "_SW":"");
  h = fopen(fname, "w+");

  fprintf(h,"Runtime %d seconds\n", end - start);

  printf("q\tH\ttau\talpha\tf\n");
  fprintf(h, "q\tH\ttau\talpha\tf\n");
  nq = 0;
  for (q = qmin; q < qmax; q += dq) {
	  tau[nq] = q * H[nq] - 3.0;
	  alpha[nq] = H[nq] + q * (H2[nq] - H[nq]) / eps;
	  f[nq] = q * alpha[nq] - tau[nq];
	  printf("%f\t%f\t%f\t%f\t%f\n", q, H[nq], tau[nq], alpha[nq], f[nq]);
	  fprintf(h, "%f\t%f\t%f\t%f\t%f\n", q, H[nq], tau[nq], alpha[nq], f[nq]);
	  nq++;
  }
  fclose(h);

  return 1;
}