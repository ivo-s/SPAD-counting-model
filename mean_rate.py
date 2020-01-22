#  mean_rate.py
#  
#  Calculation of mean detection rate of a SPAD exhibiting recovery
#  time, afterpulsing, and twilight pulsing.

import numpy as np
from scipy.signal import correlate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# time unit - 1 microsecond

### ~~ CONSTANTS AND DATA ~~ ###

binWidth = 81e-6	# AP histogram bin width

# precision setting - sufficient values depends on other parameters
no_iter = 30		# number of iterations
ilace_n = 3			# number of interlacing zeroes to enhance sampling

AP_mean = 0.0171	# mean number of afterpulses
rec_time = 0.024	# recovery time
alpha = 0.0032		# twilight constant

# import afterpulsing distribution (normalized)
AP_dist = np.loadtxt("AP_distributions/apDist_3605H.txt")

# implement interlacing
binWidth /= (1+ilace_n)
apM = AP_dist.reshape((1,-1))
apM = np.concatenate((apM,np.zeros((ilace_n,len(AP_dist)))))
AP_dist = apM.T.flatten()

### ~~ FUNCTIONS ~~ ###

# pad n zeroes to the right
def ZeroPad(ar, n):
	return np.pad(ar, (0, n), 'constant', constant_values=0)

# the integral transform
def OverheadIterate(f, mu, d, APmean, a, apDist):
	# AP intensity
	nuArray = apDist*APmean
	# time values
	t = np.arange(len(f))*binWidth
	# integrands
	integrand_1 = ZeroPad(f[d:], d)
	integrand_2 = (mu*binWidth + f) * np.exp(-mu*t - np.cumsum(f))
	# correlation using FFT
	corr = correlate(integrand_1,integrand_2, method='fft')[(len(f)-1):]
	# remove numerical artifacts
	res = (1-a*mu)*corr + a*mu*ZeroPad(f[d:], d) + nuArray
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

# mean rate of the full model
def MeanRate(mu,dt,APmean,a,apDist):
	# initial f
	f = np.zeros(len(apDist), dtype=np.double)
	# recovery time expressed in bins
	dead_index = int(round(rec_time/binWidth))
	for i in range(no_iter):
		f = OverheadIterate(f, mu,dead_index,APmean,a,apDist)
	return MeanRateFromOverhead(f, mu, dt, a)

# mean rate of the inter-arrival model
def MeanRate1(mu,dt,APmean,a,apDist):
	t = np.arange(len(apDist))*binWidth
	pa = 1. - np.exp(-APmean)
	return 1./((1./mu - a)*(1. - pa*np.dot(apDist,np.exp(-mu*t))) + dt)

# mean rate of the simple model of instant afterpulses
def MeanRate2(mu,dt,APmean,a,apDist):
	pa = 1. - np.exp(-APmean)*(1-mu*a)
	return 1./((1 - pa)/mu + dt)


### ~~ GRAPH PLOTTING ~~ ###

# sample incident rates from 0.001 to 100.
# (1k to 100M)
incLogRates = np.arange(-3.,2.01,0.2)
incRates = np.power(10., incLogRates)

# mean rate formula combining the commonly used correction factors
def MeanRateNaive(mu,dt,APmean,a,apDist):
	pa = 1. - np.exp(-APmean)*(1-mu*a)
	return (1. + pa)/(1/mu + dt)

# calculating the detection rates for all formulas
results = [ [func(r,rec_time,AP_mean,alpha,AP_dist) for r in incRates]\
			for func in (MeanRate,MeanRate1,MeanRate2,MeanRateNaive) ]
results = np.asarray(results)

# interpolate in log-scale for a nicer curve
ip = [	interp1d(incLogRates,res/results[0], 'quadratic')\
		for res in results]

# finer sampling for the interpolation
xlogvals = np.arange(-3.,2.01,0.05)
plotVals = [i(xlogvals) for i in ip]
labels = [	r"$\mu_{\mathrm{det}}$", r"$\mu_{\mathrm{det}}^{(1)}$",\
			r"$\mu_{\mathrm{det}}^{(2)}$","naive correction"]
colors = ['#003049', '#d62828', '#f77f00', '#fcbf49']

fig, ax = plt.subplots()

# plot each dataset
for (y, pts, l, clr) in zip(plotVals, results/results[0], labels, colors):
	ax.plot(xlogvals, y, color=clr, label=l)
	ax.plot(incLogRates, pts, color=clr, marker='o', linestyle='')

# set labels and ticks
ax.set(xlabel='incident rate', ylabel='detection rate (rel.)')
ax.set_ylim(0.99875,1.0015)
ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.4f}"))
ax.xaxis.set_major_locator(ticker.FixedLocator(range(-3,3)))
ax.set_xticklabels([r"$10^{%d}$" % (x) for x in range(-3,3)])

ax.set(title="Mean detection rate")

ax.legend(loc='upper left')

plt.tight_layout()

plt.show()
