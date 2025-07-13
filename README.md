# 1D, 2D & 3D Multifractal detrended fluctuation analysis - MFDFA

## Introduction
This code repository implements standard 1D MFDFA algorithm (https://doi.org/10.1016/S0378-4371(02)01383-3) in C, R and Python, representing a unified platform for these software environments.

The repository was also generalized for 2D and 3D MFDFA (https://doi.org/10.1103/PhysRevE.74.061104), with code ONLY in C (no R and Python wrappers).

## Setup
The 1D library is written in C, self-contained in the ```mfdfa.h``` header with examples of integration in ```C```, ```R```, and ```Python```, found in the respective directories. Below are instructions for running examples in each environment.
* ```C```: Run the Microsoft Visual Studio project.
* ```Python```: Compile the dynamic library  ```python setup.py build```, copy it over ```cp build/lib/mfdfa*```, run the example ```python mfdfa_benchmark.py``` (precompiled Windows dynamic link library ```mfdfa.cp311-win_amd64.pyd``` is also provided)
* ```R```: Run the example ```mfdfa_benchmark.R``` using the precompiled Windows libraries (Visual Studio project and code for compiling the dlls is also provided for portability).

Besides using the common code in ```mfdfa.h```, the 2D and 3D library code is contained in ```mfdfa2d.h``` and ```mfdfa3d.h``` headers, the MSVC projects with code and data are in the C2D and C3D subdirectories. The data file for the example for 3D MFDFA, A1001_790_790_790.bin is too large for uploading here on GitHub (~4GB, 790X790X790 X-ray CT soil sample), it can be downloaded from the OneDrive link: https://1drv.ms/f/c/0f14c99298fa863e/ErkVkJrChB5Ntzkro4qHaaUBSGa78xlquAPMh3GfOD5YSA?e=htcXQm

Therefore, if you want to test the 3D MFDFA example, after downloading this repository and compiling the MSVC project, download the A1001_790_790_790.bin into the C3D/data subdirectory (or, download the whole repository from the above OneDrive link).

## Library
The 1D library exposes a single function that computes the multifractal spectrum. The API of that function is defined as follows and has similar arguments for the ```R``` and ```Python``` wrappers.
```
double mfdfa(DFA_CONFIG* cfg, int n, double* seq, 
	double qmin, double qmax, double dq, 
	double* H, double* tau, double* alpha, double* f);
```
and the 2D and 3D functions are declared as (see the corresponding example *.c files):
```
double mfdfa_2d(DFA_CONFIG* cfg, int nx, int ny, double* seq, 
	double qmin, double qmax, double dq, double eps, double* H, double* H2);
```
```
double mfdfa_3d(DFA_CONFIG* cfg, int nx, int ny, int nz, double* seq,
	double qmin, double qmax, double dq, double eps, double* H, double* H2);
```
* ```cfg```: Configuration structure (input)
* ```n```: Number of elements in the series (input)
* ```seq```: Sequence size (input)
* ```qmin, qmax, dq```: Scaling parameter range (input)
* ```H```: Generalized Hurst exponent (output)
* ```tau```: Renyi exponent (output)
* ```alpha```: alpha (output)
* ```f```: Spectrum f (output)

The configuration structure is defined in ```mfdfa.h``` header as:
```

#define MAX_BOX	200			// maximum 200 points on logarithmic scale...
#define MAXQ	201			// max q resolution -10,...,10 dq=0.1

typedef struct {			// configurateion structure, holds all parameters
	int minbox;			// minimum box size
	int maxbox;			// maximum box size
	double boxratio;		// multiplicative factor for box size
	int rs[MAX_BOX];		// box size array 
	double *x;			// absicssa for fitting
	double mse[MAX_BOX];		// fluctuation array
	double dmse[MAXQ][MAX_BOX];	//detailed fluctuation info
	double dmse2[MAXQ][MAX_BOX];	//detailed fluctuation info
	int nfit;			// order of the regression fit
	int nr;				// number of box sizes 
	int sw;				// sliding window flag
	int goback;			// go backwards if no sliding window
}DFA_CONFIG;
```
## Results for 1D MFDFA
Examples are provided on synthetic series generated for the Binomial multifractal model with a=0.75 and 2^16=65536 data points (https://doi.org/10.1016/S0378-4371(02)01383-3), included under ```data/```. Below are shown the results obtained with ```R``` and ```Python``` (the C example yields only numeric data in ```C/ser16_mfdfa.txt```), the red curves represent exact analytical results for infinite series.

| R| Python|
|:-------------------------:|:-------------------------:|
| <img width="" alt="" src="./R/ser16.gif">|<img width="" alt="" src="./PYTHON/ser16.png"> 

## Results for 2D ad 3D MFDFA
Examples of a 3D X-ray CT scan soil sample (Stosic et al., "Comparative Analysis of 1D, 2D, and 3D MFDFA for 3D Soil Structure Characterization via X-ray CT", submitted to Chaos, July 2025) are used as examples. For 2D MFDFA only the first 2D 790x790 slice A1001_layer1_790_790.bin can be found in the C2D/data subdirectory, while (as already mentioned) the full 3D image A1001_790_790_790.bin (double precision, 790x790x790 binary data file), can be downloaded from the above OneDrive link. 
| 2D| 3D|
|:-------------------------:|:-------------------------:|
| <img width="" alt="" src="./C2D/data/A1001_layer1_790_790_mfdfa.png">|<img width="" alt="" src="./C3D/data/A1001_790_790_790_mfdfa.png"> 

## Citation
If you use this work in academic research, citating the following reference would be appreciated:

```
@software{borkostosic2024MFDFA,
  author = {Stosic, Borko},
  title = {Multifractal detrended fluctuation analysis software},
  url = {https://github.com/borkostosic/mfdfa},
  version = {1.0.0},
  year = {2024},
}
```

## Contact
Borko Stosic (borkostosic@gmail.com)
