// function declarations
double mse3d_grid(double* seq, int nx, int ny, int nz, int i0, int j0, int k0, int boxsize);
double mfdfa_3d(DFA_CONFIG* cfg, int nx, int ny, int nz, double* seq,
	double qmin, double qmax, double dq, double eps, double* H, double* H2);

/******************************************************************************************/
double mfdfa_3d(DFA_CONFIG* cfg, int nx, int ny, int nz, double* seq, 
	double qmin, double qmax, double dq, double eps, double* H, double* H2)
{
	long n, boxsize, inc, i, j, k, nq;
	double stat, chisq, a[2], f2snu, alpha;
	double* temp, q, bscale;

	for (n = 0; n < cfg->nr; n++) {
		boxsize = cfg->rs[n];
		bscale = 1.0 / (boxsize * boxsize * boxsize);
		inc = cfg->sw ? 1 : boxsize; stat = 0;
		cfg->mse[n] = 0.0;

		printf("%d\t%d out of %d\n", boxsize, n, cfg->nr);

		// first corner
		for (i = 0; i < nx - boxsize; i += inc) {				//front
			for (j = 0; j < ny - boxsize; j += inc) {			//front
				for (k = 0; k < nz - boxsize; k += inc) {		//front
					f2snu = mse3d_grid(seq, nx, ny, nz, i, j, k, boxsize, integrate);
					nq = 0;
					for (q = qmin; q < qmax; q += dq) {
						cfg->dmse[nq][n] += pow(f2snu * bscale, q / 2);
						cfg->dmse2[nq][n] += pow(f2snu * bscale, (q + eps) / 2);
						nq++;
					}
					cfg->mse[n] += f2snu;
					stat++;
				}
			}
		}

		if (!cfg->sw && cfg->goback) {
			// second corner
			for (i = nx - boxsize; i >= 0; i -= inc) {			//back
				for (j = 0; j < ny - boxsize; j += inc) {		//front
					for (k = 0; k < nz - boxsize; k += inc) {	//front
						f2snu = mse3d_grid(seq, nx, ny, nz, i, j, k, boxsize, integrate);
						nq = 0;
						for (q = qmin; q < qmax; q += dq) {
							cfg->dmse[nq][n] += pow(f2snu * bscale, q / 2);
							cfg->dmse2[nq][n] += pow(f2snu * bscale, (q + eps) / 2);
							nq++;
						}
						cfg->mse[n] += f2snu;
						stat++;
					}
				}
			}

			// third corner
			for (i = 0; i < nx - boxsize; i += inc) {			//front
				for (j = ny - boxsize; j >= 0; j -= inc) {		//back
					for (k = 0; k < nz - boxsize; k += inc) {	//front
						f2snu = mse3d_grid(seq, nx, ny, nz, i, j, k, boxsize, integrate);
						nq = 0;
						for (q = qmin; q < qmax; q += dq) {
							cfg->dmse[nq][n] += pow(f2snu * bscale, q / 2);
							cfg->dmse2[nq][n] += pow(f2snu * bscale, (q + eps) / 2);
							nq++;
						}
						cfg->mse[n] += f2snu;
						stat++;
					}
				}
			}

			// fourth corner
			for (i = 0; i < nx - boxsize; i += inc) {			//front
				for (j = 0; j < ny - boxsize; j += inc) {		//front
					for (k = nz - boxsize; k >= 0; k -= inc) {	//back
						f2snu = mse3d_grid(seq, nx, ny, nz, i, j, k, boxsize, integrate);
						nq = 0;
						for (q = qmin; q < qmax; q += dq) {
							cfg->dmse[nq][n] += pow(f2snu * bscale, q / 2);
							cfg->dmse2[nq][n] += pow(f2snu * bscale, (q + eps) / 2);
							nq++;
						}
						cfg->mse[n] += f2snu;
						stat++;
					}
				}
			}

			// fifth corner
			for (i = nx - boxsize; i >= 0; i -= inc) {			//back
				for (j = ny - boxsize; j >= 0; j -= inc) {		//back
					for (k = 0; k < nz - boxsize; k += inc) {	//front
						f2snu = mse3d_grid(seq, nx, ny, nz, i, j, k, boxsize, integrate);
						nq = 0;
						for (q = qmin; q < qmax; q += dq) {
							cfg->dmse[nq][n] += pow(f2snu * bscale, q / 2);
							cfg->dmse2[nq][n] += pow(f2snu * bscale, (q + eps) / 2);
							nq++;
						}
						cfg->mse[n] += f2snu;
						stat++;
					}
				}
			}

			// sixth corner
			for (i = nx - boxsize; i >= 0; i -= inc) {			//back
				for (j = 0; j < ny - boxsize; j += inc) {		//front
					for (k = nz - boxsize; k >= 0; k -= inc) {	//back
						f2snu = mse3d_grid(seq, nx, ny, nz, i, j, k, boxsize, integrate);
						nq = 0;
						for (q = qmin; q < qmax; q += dq) {
							cfg->dmse[nq][n] += pow(f2snu * bscale, q / 2);
							cfg->dmse2[nq][n] += pow(f2snu * bscale, (q + eps) / 2);
							nq++;
						}
						cfg->mse[n] += f2snu;
						stat++;
					}
				}
			}

			// seventh corner
			for (i = 0; i < nx - boxsize; i += inc) {			//front
				for (j = ny - boxsize; j >= 0; j -= inc) {		//back
					for (k = nz - boxsize; k >= 0; k -= inc) {	//back
						f2snu = mse3d_grid(seq, nx, ny, nz, i, j, k, boxsize, integrate);
						nq = 0;
						for (q = qmin; q < qmax; q += dq) {
							cfg->dmse[nq][n] += pow(f2snu * bscale, q / 2);
							cfg->dmse2[nq][n] += pow(f2snu * bscale, (q + eps) / 2);
							nq++;
						}
						cfg->mse[n] += f2snu;
						stat++;
					}
				}
			}

			// eight corner
			for (i = nx - boxsize; i >= 0; i -= inc) {			//back
				for (j = ny - boxsize; j >= 0; j -= inc) {		//back
					for (k = nz - boxsize; k >= 0; k -= inc) {	//back
						f2snu = mse3d_grid(seq, nx, ny, nz, i, j, k, boxsize, integrate);
						nq = 0;
						for (q = qmin; q < qmax; q += dq) {
							cfg->dmse[nq][n] += pow(f2snu * bscale, q / 2);
							cfg->dmse2[nq][n] += pow(f2snu * bscale, (q + eps) / 2);
							nq++;
						}
						cfg->mse[n] += f2snu;
						stat++;
					}
				}
			}
		}

		nq = 0;
		for (q = qmin; q < qmax; q += dq) {
			cfg->dmse[nq][n] = pow(cfg->dmse[nq][n] / stat, 1 / q);
			cfg->dmse2[nq][n] = pow(cfg->dmse2[nq][n] / stat, 1 / (q + eps));
			nq++;
		}

		cfg->mse[n] = sqrt(cfg->mse[n] / (stat * boxsize * boxsize * boxsize));
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

/*****************************************************************************************************
the following Maple script produces formulas for coefficients of the fit: a+b*x+c*y+d*z :

f := proc (x, y, z) options operator, arrow; a+b*x+c*y+d*z end proc;
s := sum((v[i]-f(x[i], y[i], z[i]))^2, i = 1 .. n);
eq1 := diff(s, a) = 0;
eq2 := diff(s, b) = 0;
eq3 := diff(s, c) = 0;
eq4 := diff(s, d) = 0;
sol := simplify(solve({eq1, eq2, eq3, eq4}, {a, b, c, d}));
with(CodeGeneration);
sol := subs(sum(y[i]^2, i = 1 .. n) = Y2, subs(sum(x[i]^2, i = 1 .. n) = X2, subs(sum(y[i]*v[i], i = 1 .. n) = VY, subs(sum(x[i]*v[i], i = 1 .. n) = VX, subs(sum(x[i]*y[i], i = 1 .. n) = XY, subs(sum(v[i], i = 1 .. n) = V, subs(sum(x[i], i = 1 .. n) = X, subs(sum(y[i], i = 1 .. n) = Y, subs(sum(z[i], i = 1 .. n) = Z, subs(sum(z[i]^2, i = 1 .. n) = Z2, subs(sum(x[i]*z[i], i = 1 .. n) = XZ, subs(sum(y[i]*z[i], i = 1 .. n) = YZ, subs(sum(v[i]*z[i], i = 1 .. n) = VZ, sol)))))))))))));
C(sol[1]);
C(sol[2]);
C(sol[3]);
C(sol[4]);
*****************************************************************************************************/
double mse3d_grid(double* seq, int nx, int ny, int nz, int i0, int j0, int k0, int boxsize) {
	int i, j, k, n, m, l, ofs, curofs;
	double V = 0, VI = 0, VJ = 0, VK = 0, val;
	double a, b, c, d, f, var, mse;
	double coef[3], dist;

	for (i = 0; i < boxsize; i++) {
		for (j = 0; j < boxsize; j++) {
			ofs = (i + i0) * ny * nz + j0 * nz + k0;
			for (k = 0; k < boxsize; k++) {
				val = seq[ofs];
				V += val;
				VI += val * (i + 1);
				VJ += val * (j + 1);
				VK += val * (k + 1);
				ofs++;
			}
		}
	}

	n = m = l = boxsize;

	a = -2 * (3 * l * m * VI - 3 * l * VI + V * l - 3 * l * VJ - 5 * m * n * V * l + 3 * l * n * VJ + 2 * n * V * l + 2 * m * V * l + 3 * m * n * VK - 3 * m * VK - 3 * m * VI - 3 * n * VK - 3 * n * VJ + 3 * VJ - 4 * V + 3 * VI + 3 * VK + V * m + V * n + 2 * n * V * m) / (l - 1) / (m - 1) / (n - 1) / n / m / l;
	b = 6 * (2 * VI - V * n - V) / n / m / l / (-1 + n * n);
	c = -6 * (-2 * VJ + V * m + V) / n / m / l / (-1 + m * m);
	d = -6 * (-2 * VK + V * l + V) / n / m / l / (-1 + l * l);

	mse = 0;
	for (i = 0; i < boxsize; i++) {
		ofs = i * boxsize;
		for (j = 0; j < boxsize; j++) {
			ofs = (i + i0) * ny * nz + j0 * nz + k0;
			for (k = 0; k < boxsize; k++) {
				val = seq[ofs];
				f = a + b * (i + 1) + c * (j + 1) + d * (k + 1);
				dist = val - f;
				mse += dist * dist;
				ofs++;
			}
		}
	}

	var = mse / ((double)boxsize * boxsize * boxsize);	// variance, not used, but might be useful for something...	
	return mse;		// square deviation
}