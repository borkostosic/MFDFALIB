#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../mfdfa.h"

double data[10000000];

/*********************************************************************************/
int load(char* filename, double* buf)
{
	/* local variables */
	long size = 0L;
	double val;
	FILE* fp;

	fp = fopen(filename, "r+");						/* open file */

	while (fscanf(fp, "%lf", &val) == 1) {			/* loop until end of file */
		size++;
		buf[size - 1] = val;							/* store data */
	}

	fclose(fp);										/* close file */
	return (int)(size);
}

/*********************************************************************************/

void main() {
  char fname[256], fnom[] = "ser16";
  int i, w, m, n;
  double r;
  double q, qmin, qmax, dq;
  int window = 250, jump = 1, total, first, offset = 0, nlin=0;
  FILE *h;
  DFA_CONFIG cfg = { 0 };		// configuration structure
  double H[MAXQ], tau[MAXQ], f[MAXQ], alpha[MAXQ];
  int nq, integrate;

  qmin = -10.000001; qmax = 10.0; dq = 0.1;

  sprintf(fname, "../data/%s.txt", fnom);
  total = load(fname, data);

  printf("%d data read from %s\n",total,fname);

  cfg.minbox = 4;
  cfg.maxbox = total / 4;
  cfg.boxratio = pow(2.0, 1.0 / 8.0);
  rscale(&cfg);		// prepares scale, allocates and fills in x values
  cfg.nfit = 1;		// nfit==1 means linear fit, nfit==2 means 2nd degree polynomial, and nfit==3 is 3rd degree
  cfg.goback = 1;	// if non-overlapping windows, go backwards as well
  //cfg.sw = 1;	// sliding window (more statistics, but slower...)

  //call the MFDFA algorithm
  mfdfa(&cfg, total, data, qmin, qmax, dq, H, tau, alpha, f);

  // print out the results to screen and file
  sprintf(fname, "%s_mfdfa%s.txt", fnom, cfg.sw ? "_SW":"");
  h = fopen(fname, "w+");
  printf("q\tH\ttau\talpha\tf\n");
  fprintf(h, "q\tH\ttau\talpha\tf\n");
  nq = 0;
  for (q = qmin; q < qmax; q += dq) {
	  printf("%f\t%f\t%f\t%f\t%f\n", q, H[nq], tau[nq], alpha[nq], f[nq]);
	  fprintf(h, "%f\t%f\t%f\t%f\t%f\n", q, H[nq], tau[nq], alpha[nq], f[nq]);
	  fflush(h);
	  nq++;
  }
  fclose(h);

}