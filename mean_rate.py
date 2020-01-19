#  mean_rate.py
#  
#  Calculation of mean detection rate of a SPAD exhibiting recovery
#  time, afterpulsing, and twilight pulsing.

import numpy as np
from scipy.fftpack import fft,ifft

# time unit - 1 microsecond

binWidth = 81e-6	# AP histogram bin width

# precision setting - sufficient values depends on other parameters
no_iter = 30		# number of iterations
ilace_n = 3			# number of interlacing zeroes to enhance sampling

# need to recalculate binWidth
binWidth /= (1+ilace_n)

meanRate = 1.		# incident rate
AP_mean = 0.0171	# mean number of afterpulses
rec_time = 0.024	# recovery time
alpha = 0.0032		# twilight constant

# import afterpulsing distribution (normalized)
AP_dist = np.loadtxt("AP_distributions/apDist_3605H.txt")


# pad n zeroes to the right
def ZeroPad(ar, n):
	return np.pad(ar, (0, n), 'constant', constant_values=0)

# the convolution integral
def OverheadConvolve(f, mu, d, APmean, a, apDist):
	# interlace with zeroes to achieve finer sampling
	apM = apDist.reshape((1,-1))
	apM = np.concatenate((apM,np.zeros((ilace_n,len(apDist)))))
	apDist = apM.T.flatten()
	# AP intensity
	nuArray = apDist*APmean
	# time values
	t = np.arange(len(f))*binWidth
	# convolution integrands
	convol_elem_1 = ZeroPad(f[d:], len(f)+d)
	convol_elem_2 = ZeroPad((mu*binWidth + f) * \
					np.exp(-mu*t - np.cumsum(f)), len(f))
	# convolution using FFT
	convol = ifft(fft(convol_elem_1)*np.conjugate(fft(convol_elem_2)))
	# remove numerical artifacts
	convol = np.real_if_close(convol[:len(f)])
	res = (1-a*mu)*convol + a*mu*ZeroPad(f[d:], d) + nuArray
	return res

# calculate mean rate from f
def MeanRateFromOverhead(f, mu, dt, a):
	# time values
	t = np.arange(len(f))*binWidth
	# stationary interarrival PDF
	p_bar = (mu*binWidth + f)*np.exp(-mu*t - np.cumsum(f))
	# calculate an auxiliary term x
	x = np.dot(t,p_bar) + ((1./mu + (t[-1]+binWidth)) * \
		np.exp(-(t[-1]+binWidth)*mu - sum(f)))
	return (1./(x*(1.-a*mu) + dt))

def MeanRate(mu,dt,APmean,a,apDist):
	# initial f
	f = np.zeros(len(apDist)*(1+ilace_n), dtype=np.double)
	# recovery time expressed in bins
	dead_index = int(round(rec_time/binWidth))
	for i in range(no_iter):
		f = OverheadConvolve(f, mu,dead_index,APmean,a,apDist)
	return MeanRateFromOverhead(f, mu, dt, a)

print MeanRate(meanRate,rec_time,AP_mean,alpha,AP_dist)
