import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import time
import sys
def printf(format, *args):
    sys.stdout.write(format % args)

import numpy as np
import mfdfa

MAX_BOX=200		# maximum 200 points on logarithmic scale...
MAXQ=201		# max q resolution -10,...,10 dq=0.1

# Load data.
series = 'ser16'	# binomial multifractal series 2^16=65536 data points
data = np.loadtxt('../data/'+series+'.txt')
total = data.shape[0]
printf("Loaded %d data points from %s\n", total, '../data/'+series+'.txt')

# Parameters.
dq = 0.1
qmin = -10.000001	# small offset to avoid 0
qmax = 10
n = int(np.abs(qmax - qmin)/dq)+1

# arrays for results
H = np.empty(n)
tau = np.empty(n)
f = np.empty(n)
alpha = np.empty(n)
dmse = np.ndarray(shape = [MAXQ, MAX_BOX])	# maximum 200 points on logarithmic scale, and 201 maximum q values
nfit=1		# polynomial degree 1 (linear), 2 or 3
sw=0		# sw=1 for sliding segments, takes longer, use with care...

rs=np.zeros(MAX_BOX)			# segment size list
################ prepare multiplicative scale (equidistant on log-log plot)
nrs=0
minseg=4		# minimum segment size
maxseg=int(total/4)	# maximum segment size
boxratio=np.power(2,1.0/8)	# segment ratio (multiplicative scale)
################...OR prepare your own segment scale
#rs=np.arange(10, 101, 1) 	
#nrs=len(rs)

#run the MFDFA C dynamic link library code
t0 = time.time()
nr=mfdfa.calc(data,total,H,tau,f,alpha,dmse,rs,nrs,qmin,qmax,dq,minseg,maxseg,boxratio,nfit,sw)
t1 = time.time()
runtime = t1-t0
printf("MFDFA for q from %.2f to %.2f with step %.2f runtime=%.4f sec\n", qmin, qmax, dq, runtime)

# Do the plots
fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(2, 2)

qlist = []
for q in np.arange(qmin,qmax,dq):
    qlist.append(q)
q = np.asarray(qlist)

######################################## tau(q) plot ################################
ax = fig.add_subplot(gs[0, 0])
ax.plot(q,tau, marker='o')
ax.set_ylabel('tau(q)')
ax.set_xlabel('q')

######################################## F(s) plot ################################
ax = fig.add_subplot(gs[0, 1])
for i in range(0, n):
    ax.loglog(rs[0:nr], dmse[i,0:nr], 'o', label = q[i])
ax.set_ylabel('F(s)')
ax.set_xlabel('s')

############################## calculate theoretical curve ########################
tot=100
f_t = np.empty(tot)
alpha_t = np.empty(tot)
for i in range(0, tot):
    iq=qmin+i*(qmax-qmin)/tot
    f_t[i]= 0.2000000000e-9 * (0.2075187497e10 * iq * np.exp(-0.2876820725e0 * iq) + 0.1000000000e11 * iq * np.exp(-0.1386294361e1 * iq) + 0.7213475205e10 * np.log(np.exp(-0.2876820725e0 * iq) + np.exp(-0.1386294361e1 * iq)) * np.exp(-0.2876820725e0 * iq) + 0.7213475205e10 * np.log(np.exp(-0.2876820725e0 * iq) + np.exp(-0.1386294361e1 * iq)) * np.exp(-0.1386294361e1 * iq)) / (np.exp(-0.2876820725e0 * iq) + np.exp(-0.1386294361e1 * iq))
    alpha_t[i]= 0.1e1 / iq - np.log(np.power(0.75e0,iq) + np.power(0.25e0,iq)) / iq / np.log(0.2e1) + iq * (-np.power(iq,-0.2e1) - (-0.2876820725e0 * np.power(0.75e0,iq) - 0.1386294361e1 * np.power(0.25e0,iq)) / (np.power(0.75e0,iq) + np.power(0.25e0,iq)) / iq / np.log(0.2e1) + np.log(np.power(0.75e0,iq) + np.power(0.25e0,iq)) * np.power(iq,-0.2e1) / np.log(0.2e1))

#print(f_t)
#print(alpha_t)

######################################## f(alpha) plot ################################
ax = fig.add_subplot(gs[1, 0])
ax.scatter(alpha, f, s=40, facecolors='none', edgecolors='b')	# MFDFA result in blue
ax.plot(alpha_t, f_t, 'r')					# theoretical curve red
ax.set_ylabel('f(alpha)')
ax.set_xlabel('alpha')

######################################## H(q) plot ################################
ax = fig.add_subplot(gs[1, 1])
ax.plot(q,H, marker='o')
ax.set_ylabel('H(q)')
ax.set_xlabel('q')

plt.show()
fig.savefig(series + '.png')

