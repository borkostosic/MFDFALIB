#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#define MAX_BOX	200			// maximum 200 points on logarithmic scale...
#define MAXQ	201			// max q resolution -10,...,10 dq=0.1

typedef struct {			// configurateion structure, holds all parameters
	int minbox;				// minimum box size
	int maxbox;				// maximum box size
	double boxratio;		// multiplicative factor for box size
	int rs[MAX_BOX];		// box size array 
	double *x;					// absicssa for fitting
	double mse[MAX_BOX];	// fluctuation array
	double dmse[MAXQ][MAX_BOX];	//detailed fluctuation info
	double dmse2[MAXQ][MAX_BOX];	//detailed fluctuation info
	int nfit;				// order of the regression fit
	int nr;					// number of box sizes 
	int sw;					// sliding window flag
	int goback;				// go backwards
}DFA_CONFIG;

// function declarations
int rscale(DFA_CONFIG *cfg);
double polyfit(double *x, double *y, int boxsize, int nfit);
double fit_log(int *x, double *y, int n, double *a, double *chisq);
double fit_poly1(double *x, double *y, int n, double *a, double *chisq);
double fit_poly2(double *x, double *y, int n, double *a, double *chisq);
double fit_poly3(double *x, double *y, int n, double *a, double *chisq);


double mfdfa(DFA_CONFIG* cfg, int nx, double* seq, 
	double qmin, double qmax, double dq, 
	double* H, double* tau, double* alpha, double* f);

double eps = 0.001;
int integrate=1;
/************************************************************************************************/

/* modified from C.K.Peng's original code to use the DFA_CONFIG structure and 0 offset for rs array:
rscale() allocates and fills rs[], the array of box sizes used by mfdfa()
below.  The box sizes range from (exactly) minbox to (approximately) maxbox,
and are arranged in a geometric series such that the ratio between
consecutive box sizes is (approximately) boxratio.  The return value is
the number of box sizes in rs[].
*/
int rscale(DFA_CONFIG *cfg)
{
	int i, ir, n, tot;
	long rw, rslen;

	/* Determine how many scales are needed. */
	rslen = (int)(log10(cfg->maxbox / (double)cfg->minbox) / log10(cfg->boxratio) + 1.5);
	/* Thanks to Peter Domitrovich for pointing out that a previous version
	of the above calculation undercounted the number of scales in some situations. */
	for (ir = 1, n = 1, cfg->rs[0] = cfg->minbox; n <= rslen && cfg->rs[n - 1] < cfg->maxbox; ir++)
		if ((rw = cfg->minbox * pow(cfg->boxratio, ir) + 0.5) > cfg->rs[n - 1])
			cfg->rs[n++] = rw;
	if (cfg->rs[n-1] > cfg->maxbox) --n;
	cfg->nr = n;
	tot = cfg->rs[n - 1];
	cfg->x = (double *)malloc(tot * sizeof(double));
	if (!cfg->x)
		return 0;
	for (i = 0; i < tot; i++)
		cfg->x[i] = i+1;
	return (n);
}

/*************** wrapper function for polynomial fit *****************************/
double polyfit(double *x, double *y, int boxsize, int nfit)
{
	double mse, chisq, a[5];
	switch (nfit)
	{
	case 1:
		mse=fit_poly1(x, y, boxsize, a, &chisq);	// linear fit
		break;
	case 2:
		mse=fit_poly2(x, y, boxsize, a, &chisq);	// second order polynomil
		break;
	case 3:
		mse=fit_poly3(x, y, boxsize, a, &chisq);	// third order polynomial
		break;
	default:
		return -1.0;
	}
	return mse;
}

/*************** linear fit *************************************/
double fit_poly1(double *x, double *y, int n, double *a, double *chisq)
{
double sx=0, sx2=0, sy=0, sxy=0, temp, mse;
int i;

for(i=0;i<n;i++)
	{
	sx+=x[i];
	sy+=y[i];
	sx2+=x[i]*x[i];
	sxy+=x[i]*y[i];
	}
a[1]=(n*sxy-sx*sy)/(n*sx2-sx*sx);
a[0]=(sy-(a[1])*sx)/n;

*chisq=0;
for(i=0;i<n;i++)
	{
	temp=a[1]*x[i]+a[0]-y[i];
	*chisq+=temp*temp;
	}
mse=*chisq;
*chisq /= n;
return mse;
}

/*************** fit log log scale *************************************/
double fit_log(int *x, double *y, int n, double *a, double *chisq)
{
	double sx = 0, sx2 = 0, sy = 0, sxy = 0, temp, mse;
	double logx, logy;
	int i;

	for (i = 0; i<n; i++)
	{
		logx = log10(x[i]);
		logy = log10(y[i])/2;
#ifdef DUMP
		fprintf(h, "%f\t%f\n", logx, logy);
#endif
		sx += logx;
		sy += logy;
		sx2 += logx*logx;
		sxy += logx*logy;
	}
	a[1] = (n*sxy - sx*sy) / (n*sx2 - sx*sx);
	a[0] = (sy - (a[1])*sx) / n;

	return a[1];
}

/************** 2nd order polynomial fit ************************************/
double fit_poly2(double *x, double *y, int n, double *a, double *chisq)
{
double sx, sy, sx2, sxy, sx3, sx4, syx2, temp, xx, mse;
int i;

sx=sy=sx2=sxy=sx3=sx4=syx2=0;
for(i=0;i<n;i++)
	{
	sx+=x[i];
	sy+=y[i];
	xx=x[i]*x[i];
	sx2+=xx;
	sxy+=x[i]*y[i];
	syx2+=y[i]*xx;
	sx3+=xx*x[i];
	sx4+=xx*xx;
	}

temp=(n*sx3*sx3-n*sx2*sx4+sx2*sx2*sx2+sx*sx*sx4-2.0*sx3*sx*sx2);
if(temp==0.0)
	return 0;
a[1] = (-sx*sx2*syx2+sx*sy*sx4+sx2*sx2*sxy+n*sx3*syx2-sx2*sx3*sy-n*sxy*sx4)/temp;
a[0] = -(-sy*sx3*sx3+sy*sx2*sx4+sx*sx3*syx2-sx*sxy*sx4-sx2*sx2*syx2+sx2*sx3*sxy)/temp;
a[2] = -(n*sx2*syx2+sx*sx2*sxy-sx2*sx2*sy-n*sx3*sxy-sx*sx*syx2+sx3*sx*sy)/temp;

*chisq=0;
for(i=0;i<n;i++)
	{
	temp=a[2]*x[i]*x[i]+a[1]*x[i]+a[0]-y[i];
	*chisq+=temp*temp;
	}
mse = *chisq;
*chisq /= n;
return mse;
}

/************** 3rd order polynomial fit ************************************/
double fit_poly3(double *x, double *y, int n, double *a, double *chisq)
{
double sx=0, sy=0, sx2=0, syx=0, sx3=0, sx4=0, syx2=0, temp=0, xx=0, mse;
double sx5=0,sx6=0,syx3=0, xxx=0, xxxx=0;
int i;

for(i=0;i<n;i++)
	{
	sx+=x[i];
	sy+=y[i];
	xx=x[i]*x[i];
	xxx=xx*x[i];
	xxxx=xx*xx;
	sx2+=xx;
	sx3+=xxx;
	sx4+=xxxx;
	sx5+=xxx*xx;
	sx6+=xxxx*xx;
	syx+=x[i]*y[i];
	syx2+=y[i]*xx;
	syx3+=y[i]*xxx;
	}

temp=(sx4*n*sx2*sx6+sx*sx*sx5*sx5-n*sx3*
sx3*sx6+2.0*sx5*sx2*sx2*sx3+2.0*sx*sx3*sx2*sx6+2.0*sx5*n*sx3*sx4-sx2*sx2*sx2*sx6-sx*sx*sx4*sx6+sx2*sx2
*sx4*sx4-2.0*sx*sx5*sx3*sx3-sx4*sx4*sx4*n-3.0*sx3*sx3*sx2*sx4-sx5*sx5*n*sx2-2.0*sx*sx5*sx2*sx4+2.0*
sx*sx4*sx4*sx3+sx3*sx3*sx3*sx3);

a[0] = (-sy*sx4*sx4*sx4-sx2*sx4*sx5*syx+sx2*sx3*syx*sx6-sx3*sx3*sx4*syx2-sx5*sx3*sx3*syx+sx3*sx3*sx3*syx3+2.0*sy*sx5*sx3*sx4-sx*syx*sx4*sx6+sx*sx3*syx2*sx6-sx*sx4*sx5*syx2+sx*syx*sx5*sx5+sx2*sx4*sx4*
syx2+sx3*sx4*sx4*syx-sy*sx3*sx3*sx6+sx*sx4*sx4*syx3-sx*sx3*sx5*syx3-sy*sx5*sx5*sx2-sx2*sx2*syx2*sx6-2.0*
sx3*sx4*sx2*syx3+sx2*sx2*sx5*syx3+sy*sx4*sx2*sx6+sx3*sx5*sx2*syx2)/temp;

a[1] = (sx*sy*sx5*sx5+sx*sx2*syx2*sx6-sx*sx5*sx2*syx3+sx*sx4*sx3*syx3-sx*sx3*sx5*syx2-sx*sy*sx4*sx6+sx3*sx4*sx4*sy-n*sx4*sx4*syx3-sx2*sx2*syx*sx6-sx2*sx3*sx3*syx3+n*sx3*sx5*syx3-sx2*sx4*sx5*sy+n*syx*sx4
*sx6+2.0*sx2*sx5*sx3*syx-n*sx3*syx2*sx6+n*sx4*sx5*syx2+sx2*sx2*sx4*syx3-sx3*sx3*sx5*sy-n*syx*sx5*sx5-
sx3*sx3*sx4*syx+sx3*sx3*sx3*syx2-sx2*sx3*sx4*syx2+sx2*sx3*sy*sx6)/temp;

a[2] = (-sx*sx*syx2*sx6+sx*sx*sx5*syx3-sx*sx4*sx2*syx3-sx*sx4*sx5*sy+2.0*sx*sx3*sx4*syx2-sx*sx3*sx3*syx3+sx*sx2*syx*sx6-sx*sx5*sx3*syx+sx*sx3*sy*sx6-sx4*sx4*n*syx2+sx4*sx4*sx2*sy-sx4*sx2*sx3*syx-sx4*sx3
*sx3*sy+sx4*n*sx3*syx3+sx4*sx5*n*syx+n*sx2*syx2*sx6-sx2*sx2*sy*sx6+sx5*sx3*sx2*sy-sx3*sx3*sx2*syx2+sx2*
sx2*sx3*syx3-n*sx3*syx*sx6+sx3*sx3*sx3*syx-sx5*n*sx2*syx3)/temp;

a[3] = -1/(sx4*n*sx2*sx6+sx*sx*sx5*sx5-n*sx3*sx3*sx6+2.0*sx5*sx2*sx2*sx3+2.0*sx*sx3*sx2*sx6+2.0*sx5*n*sx3*sx4-sx2*sx2*sx2*sx6-sx*sx*sx4*sx6+sx2*sx2*sx4*sx4-2.0*sx*sx5*sx3*sx3-sx4*sx4*sx4*n-3.0*
sx3*sx3*sx2*sx4-sx5*sx5*n*sx2-2.0*sx*sx5*sx2*sx4+2.0*sx*sx4*sx4*sx3+sx3*sx3*sx3*sx3)*(-sx*sx*sx5*syx2+
n*sx3*sx3*syx3-n*sx3*sx4*syx2+sx4*sx4*n*syx-sx*sx4*sx4*sy-sx4*n*sx2*syx3+sx*sx3*sx3*syx2+sx2*sx2*sx2*
syx3-sx3*sx3*sx3*sy-sx5*sx2*sx2*sy+sx5*n*sx2*syx2+sx3*sx3*sx2*syx-sx2*sx2*sx3*syx2-sx2*sx2*sx4*syx-sx*sx4
*sx3*syx+sx*sx*sx4*syx3+sx*sx5*sx2*syx+2.0*sx3*sx2*sx4*sy-2.0*sx*sx3*sx2*syx3-sx5*n*sx3*syx+sx*sx5*sx3*
sy+sx*sx2*sx4*syx2);


*chisq=0;
for(i=0;i<n;i++)
	{
	temp=a[3]*x[i]*x[i]*x[i]+a[2]*x[i]*x[i]+a[1]*x[i]+a[0]-y[i];
	*chisq+=temp*temp;
	}
mse = *chisq;
*chisq /= n;
return mse;
}

/******************************************************************************************/

double mfdfa(DFA_CONFIG* cfg, int nx, double* seq, 
	double qmin, double qmax, double dq, 
	double* H, double* tau, double* alpha, double* f)
{
	long n, boxsize, inc, i, j;
	double stat, chisq, a[2], f2snu, alphaDFA;
	double *temp, q, H2[MAXQ];
	int nq;

	memset(cfg->mse, 0, sizeof(cfg->mse));
	memset(cfg->dmse, 0, sizeof(cfg->dmse));
	memset(cfg->dmse2, 0, sizeof(cfg->dmse2));

	if (integrate){
		for (i = 1; i < nx; i++)
			seq[i] += seq[i - 1];
		}

	for (n = 0; n < cfg->nr; n++) {
//		printf("%d/%d\n", n + 1, cfg->nr);
		boxsize = cfg->rs[n];
		temp = (double *)malloc(boxsize*boxsize * sizeof(double));
		inc = cfg->sw ? 1 : boxsize; stat = 0;
		cfg->mse[n] = 0.0;
		for (i = 0; i < nx - boxsize; i += inc) {
			f2snu = polyfit(cfg->x, seq + i, boxsize, cfg->nfit);

			if (!f2snu)
				continue;
			nq = 0;
			for (q = qmin; q < qmax; q += dq) {
				cfg->dmse[nq][n] += pow(f2snu / boxsize, q / 2);
				cfg->dmse2[nq][n] += pow(f2snu / boxsize, (q+eps) / 2);
				nq++;
			}
				cfg->mse[n] += f2snu/boxsize;
				stat++;
		}

		if (!cfg->sw && cfg->goback) {
			// from the back
			for (i = nx - boxsize; i >= 0; i -= inc) {
				f2snu = polyfit(cfg->x, seq + i, boxsize, cfg->nfit);
				if (!f2snu)
					continue;
				nq = 0;
				for (q = qmin; q < qmax; q += dq) {
					cfg->dmse[nq][n] += pow(f2snu / boxsize, q / 2);
					cfg->dmse2[nq][n] += pow(f2snu / boxsize, (q + eps) / 2);
					nq++;
				}
				cfg->mse[n] += f2snu;
				stat++;
			}
		}

		nq = 0;
		for (q = qmin; q < qmax; q += dq) {
			cfg->dmse[nq][n]= pow(cfg->dmse[nq][n] / stat, 2 / q);
			cfg->dmse2[nq][n] = pow(cfg->dmse2[nq][n] / stat, 2 / (q+eps));
			nq++;
		}
		cfg->mse[n] /= stat*boxsize;

		free(temp);
	}
	nq = 0;
	for (q = qmin; q < qmax; q += dq) {
		H[nq] = fit_log(cfg->rs, cfg->dmse[nq], cfg->nr, a, &chisq);
		H2[nq] = fit_log(cfg->rs, cfg->dmse2[nq], cfg->nr, a, &chisq);
		nq++;
	}

	nq = 0;
	for (q = qmin; q < qmax; q += dq) {
		tau[nq] = q * H[nq] - 1.0;
		alpha[nq] = H[nq] + q * (H2[nq] - H[nq]) / eps;
		f[nq] = q * alpha[nq] - tau[nq];
		nq++;
	}

	alphaDFA = fit_log(cfg->rs, cfg->mse, cfg->nr, a, &chisq);
	return alphaDFA;
}

