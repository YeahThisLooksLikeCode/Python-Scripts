from cmath import exp,pi
import numpy as np
from matplotlib import pyplot as plt

def fourier(x):
	N = np.size(x)
	T = np.empty(N)
	for i in range(0,N-1):
		for k in range(0,N-1):
			T[i] += exp(-2j*k*i/N) * x[i]
	return T

sample_rate = 1000
dt = 1.0/sample_rate
t = np.arange(sample_rate)*dt  # 1 second of samples
freq = 8
freq2 = 15
amp = 4.0
amp2 = 0.5
sine1 = amp*np.sin(2*np.pi*freq*t)
sine2 = amp2*np.sin(2*np.pi*freq2*t)
sinsum = sine1 + sine2 #combined sine waves to get signal
arr = fourier(sinsum)
plt.plot(arr)
plt.show()