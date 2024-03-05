#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <math.h>
#include "../mfdfa.h"

namespace py = pybind11;

int calc(py::array_t<double> data, int total,
	py::array_t<double> H, py::array_t<double> tau,  
	py::array_t<double> f,py::array_t<double> alpha,
	py::array_t<double> dmse,
	py::array_t<double> rs, int nrs,
	double qmin, double qmax, double dq,
	int minbox, int maxbox, double boxratio, int nfit, int sw) {

    DFA_CONFIG cfg = { 0 };		// configuration structure
    double* data_ptr = static_cast<double*>(data.request().ptr);
    double* H_ptr = static_cast<double*>(H.request().ptr);
    double* tau_ptr = static_cast<double*>(tau.request().ptr);
    double* f_ptr = static_cast<double*>(f.request().ptr);
    double* alpha_ptr = static_cast<double*>(alpha.request().ptr);
    double* rs_ptr = static_cast<double*>(rs.request().ptr);
	double* dmse_ptr = static_cast<double*>(dmse.request().ptr);
    int nq, nr;
    double q;

	cfg.minbox = minbox ? minbox : 4;
	cfg.maxbox = maxbox ? maxbox : total / 4;
	cfg.boxratio = boxratio ? boxratio : pow(2.0, 1.0 / 8.0);
	cfg.nfit = nfit ? nfit : 1;		// nfit==1 linear fit, nfit==2 2nd degree poly, etc.
	cfg.sw = sw ? 1 : 0;
	cfg.goback = 1;	// if non-overlapping windows, go backwards as well


	if (nrs) {				// user supplied scale
		int tot, n, i;
		cfg.nr = n = nrs;
		for (nr = 0; nr < cfg.nr; nr++) {
			cfg.rs[nr] = rs_ptr[nr];			// copy to cfg structure
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

    //call the MFDFA algorithm
    mfdfa(&cfg, total, data_ptr, qmin, qmax, dq, 
		H_ptr, tau_ptr, alpha_ptr, f_ptr);
	
	for (nr = 0; nr<cfg.nr; nr++) {
		rs_ptr[nr] = cfg.rs[nr];
	}
    
	nq = 0;
    for (q = qmin; q < qmax; q += dq) {
	  for(nr=0; nr<cfg.nr; nr++){
		dmse_ptr[nq*MAX_BOX+nr]=sqrt(cfg.dmse[nq][nr]);
	  }
	  nq++;
    }

    return cfg.nr;

}

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m)
{
    m.def("calc", &calc, "calc");
}