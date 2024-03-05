#include "..\mfdfa.h"

/********* wrapper function; only pointers can be passed from R **************/

void mfdfa_calc(double *data, int* total, 
	double *H, double* tau, double* f, double* alpha,
	double* dmse, int *rs, int *nrs, 
	double *qmin, double *qmax, double *dq,
	int *minbox, int *maxbox, double *boxratio,
	int *nfit, int *sw)
{
	int i, nq, nr;
	double q;
	DFA_CONFIG cfg = {0};

	cfg.minbox = *minbox ? *minbox : 4;
	cfg.maxbox = *maxbox ? *maxbox : *total / 4;
	cfg.boxratio = *boxratio ? *boxratio : pow(2.0, 1.0 / 8.0);
	cfg.nfit = *nfit ? *nfit : 1;		// nfit==1 linear fit, nfit==2 2nd degree poly, etc.
	cfg.sw = *sw ? 1 : 0;
	cfg.goback = 1;	// if non-overlapping windows, go backwards as well

	
	if (*nrs) {				// user supplied scale
		int tot, n;
		cfg.nr = n= *nrs;
		for (nr = 0; nr < cfg.nr; nr++) {
			cfg.rs[nr] = rs[nr];			// copy to cfg structure
		}
		tot = cfg.rs[n - 1];
		cfg.x = (double*)malloc(tot * sizeof(double));
		if (!cfg.x)
			return 0;
		for (i = 0; i < tot; i++)	// fill in x values for fitting
			cfg.x[i] = i + 1;
	}
	else
		rscale(&cfg);	// prepares scale, allocates and fills in x value
	
	mfdfa(&cfg, *total, data, 
		*qmin, *qmax, *dq, H, tau, alpha, f);	//call the MFDFA algorithm

	if (!*nrs) {
		*nrs = cfg.nr;
		for (nr = 0; nr < cfg.nr; nr++) {
			rs[nr] = cfg.rs[nr];
		}
	}

	nq = 0;
	for (q = *qmin; q < *qmax; q += *dq) {
		for (nr = 0; nr < cfg.nr; nr++) {
			dmse[nq * MAX_BOX + nr] = cfg.dmse[nq][nr];
		}
		nq++;
	}

	if (cfg.x)
		free(cfg.x);

	return cfg.nr;
}

