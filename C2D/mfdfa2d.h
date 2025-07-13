// function declarations
double mse2d_grid(double* seq, int nx, int ny, int i0, int j0, int boxsize, double* xdispl);
double mfdfa_2d(DFA_CONFIG* cfg, int nx, int ny, double* seq,
	double qmin, double qmax, double dq, double eps, double* H, double* H2);

/******************************************************************************************/

double mfdfa_2d(DFA_CONFIG* cfg, int nx, int ny, double* seq, 
	double qmin, double qmax, double dq, double eps, double* H, double* H2)
{
	long n, boxsize, inc, i, j;
	double stat, chisq, a[2], f2snu, alpha;
	double* temp, q;
	int nq;

	memset(cfg->mse, 0, sizeof(cfg->mse));
	memset(cfg->dmse, 0, sizeof(cfg->dmse));
	memset(cfg->dmse2, 0, sizeof(cfg->dmse2));

	for (n = 0; n < cfg->nr; n++) {
		boxsize = cfg->rs[n];
		temp = (double*)malloc(boxsize * boxsize * sizeof(double));
		inc = cfg->sw ? 1 : boxsize; stat = 0;
		cfg->mse[n] = 0.0;
		for (i = 0; i < nx - boxsize; i += inc) {
			for (j = 0; j < ny - boxsize; j += inc) {
				f2snu = mse2d_grid(seq, nx, ny, i, j, boxsize, temp);	// fit to a plane, returns square deviation
				nq = 0;
				for (q = qmin; q < qmax; q += dq) {
					cfg->dmse[nq][n] += pow(f2snu / (boxsize * boxsize), q / 2);
					cfg->dmse2[nq][n] += pow(f2snu / (boxsize * boxsize), (q + eps) / 2);
					nq++;
				}
				cfg->mse[n] += f2snu;
				stat++;
			}
		}

		if (!cfg->sw && cfg->goback) {						// go through all four corners
			// second corner
			for (i = nx - boxsize; i >= 0; i -= inc) {
				for (j = 0; j < ny - boxsize; j += inc) {
					f2snu = mse2d_grid(seq, nx, ny, i, j, boxsize, temp);
					nq = 0;
					for (q = qmin; q < qmax; q += dq) {
						cfg->dmse[nq][n] += pow(f2snu / (boxsize * boxsize), q / 2);
						cfg->dmse2[nq][n] += pow(f2snu / (boxsize * boxsize), (q + eps) / 2);
						nq++;
					}
					cfg->mse[n] += f2snu;
					stat++;
				}
			}
			// third corner
			for (i = 0; i < nx - boxsize; i += inc) {
				for (j = nx - boxsize; j >= 0; j -= inc) {
					f2snu = mse2d_grid(seq, nx, ny, i, j, boxsize, temp);
					nq = 0;
					for (q = qmin; q < qmax; q += dq) {
						cfg->dmse[nq][n] += pow(f2snu / (boxsize * boxsize), q / 2);
						cfg->dmse2[nq][n] += pow(f2snu / (boxsize * boxsize), (q + eps) / 2);
						nq++;
					}
					cfg->mse[n] += f2snu;
					stat++;
				}
			}
			// fourth corner
			for (i = nx - boxsize; i >= 0; i -= inc) {
				for (j = ny - boxsize; j >= 0; j -= inc) {
					f2snu = mse2d_grid(seq, nx, ny, i, j, boxsize, temp);
					nq = 0;
					for (q = qmin; q < qmax; q += dq) {
						cfg->dmse[nq][n] += pow(f2snu / (boxsize * boxsize), q / 2);
						cfg->dmse2[nq][n] += pow(f2snu / (boxsize * boxsize), (q + eps) / 2);
						nq++;
					}
					cfg->mse[n] += f2snu;
					stat++;
				}
			}
		}

		nq = 0;
		for (q = qmin; q < qmax; q += dq) {
			cfg->dmse[nq][n] = pow(cfg->dmse[nq][n] / stat, 1 / q);
			cfg->dmse2[nq][n] = pow(cfg->dmse2[nq][n] / stat, 1 / (q + eps));
			nq++;
		}
		cfg->mse[n] = sqrt(cfg->mse[n] / (stat * (boxsize * boxsize)));

		free(temp);
	}

	nq = 0;
	for (q = qmin; q < qmax; q += dq) {
		H[nq] = 2*fit_log(cfg->rs, cfg->dmse[nq], cfg->nr, a, &chisq);
		H2[nq] = 2*fit_log(cfg->rs, cfg->dmse2[nq], cfg->nr, a, &chisq);
		nq++;
	}
	alpha = 2*fit_log(cfg->rs, cfg->mse, cfg->nr, a, &chisq);
	return alpha;
}

/********* fit data to a plane in (x,y) **************************************************************
the following Maple script produces formulas for coefficients of the fit: a + b * x + c * y:

f: = proc(x, y) options operator, arrow; a + b * x + c * y end proc;
s: = sum((v[i] - f(x[i], y[i])) ^ 2, i = 1 ..n);
eq1: = diff(s, a) = 0;
eq2: = diff(s, b) = 0;
eq3: = diff(s, c) = 0;
sol: = simplify(solve({ eq1, eq2, eq3 }, { a, b, c }));
	with(CodeGeneration);
ssol: = subs(sum(y[i] ^ 2, i = 1 ..n) = Y2, subs(sum(x[i] ^ 2, i = 1 ..n) = X2, subs(sum(v[i] * y[i], i = 1 ..n) = YV, subs(sum(v[i] * x[i], i = 1 ..n) = XV, subs(sum(x[i] * y[i], i = 1 ..n) = XY, subs(sum(v[i], i = 1 ..n) = V, subs(sum(x[i], i = 1 ..n) = X, subs(sum(y[i], i = 1 ..n) = Y, sol))))))));
	C(ssol[1]);
	C(ssol[2]);
	C(ssol[3]);
*****************************************************************************************************/

double mse2d_grid(double* seq, int nx, int ny, int i0, int j0, int boxsize, double* xdispl) {
	int i, j, n, m, ofs, curofs;
	double V = 0, VI = 0, VJ = 0, temp, val;
	double a, b, c, f, var, mse;
	double coef[3], avgdist, dist;

	for (i = 0; i < boxsize; i++) {
		temp = 0;
		ofs = (i + i0) * nx + j0;
		for (j = 0; j < boxsize; j++) {
			val = seq[ofs];
			V += val;
			VJ += val * (j + 1);
			temp += val;
			ofs++;
		}
		VI += temp * (i + 1);
	}

	n = m = boxsize;
	a = -(6 * m * VI - 7 * n * V * m + V * m + V * n - 6 * VI - 6 * VJ + 6 * n * VJ + 5 * V) / (m - 1) / (n - 1) / n / m;
	b = 6 * (2 * VI - V - V * n) / n / m / (-1 + n * n);
	c = -6 * (-2 * VJ + V * m + V) / n / m / (m * m - 1);

	mse = 0;	
	for (i = 0; i < boxsize; i++) {
		ofs = (i + i0) * nx + j0;
		for (j = 0; j < boxsize; j++) {
			val = seq[ofs];
			f = a + b * (i + 1) + c * (j + 1);
			dist = val - f;
			mse += dist * dist;
			ofs++;
		}
	}
	var = mse / (boxsize * boxsize);	// variance, not used, but might be useful for something...
	
	return mse;		// square deviation
}
