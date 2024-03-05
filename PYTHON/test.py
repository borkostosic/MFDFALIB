import numpy as np

qmin = -10.000001	# small offset to avoid 0
qmax = 10

tot=100
f_t = np.empty(tot)
alpha_t = np.empty(tot)
for i in range(1, tot):
    q=qmin+i*(qmax-qmin)/tot
    f_t[i]<- 0.2000000000e-9 * (0.2075187497e10 * q * np.exp(-0.2876820725e0 * q) + 0.1000000000e11 * q * np.exp(-0.1386294361e1 * q) + 0.7213475205e10 * np.log(np.exp(-0.2876820725e0 * q) + np.exp(-0.1386294361e1 * q)) * np.exp(-0.2876820725e0 * q) + 0.7213475205e10 * np.log(np.exp(-0.2876820725e0 * q) + np.exp(-0.1386294361e1 * q)) * np.exp(-0.1386294361e1 * q)) / (np.exp(-0.2876820725e0 * q) + np.exp(-0.1386294361e1 * q))

    alpha_t[i]<- 0.1e1 / q - np.log(np.power(0.75e0,q) + np.power(0.25e0,q)) / q / np.log(0.2e1) + q * (-np.power(q,-0.2e1) - (-0.2876820725e0 * np.power(0.75e0,q) - 0.1386294361e1 * np.power(0.25e0,q)) / (np.power(0.75e0,q) + np.power(0.25e0,q)) / q / np.log(0.2e1) + np.log(np.power(0.75e0,q) + np.power(0.25e0,q)) * np.power(q,-0.2e1) / np.log(0.2e1))
