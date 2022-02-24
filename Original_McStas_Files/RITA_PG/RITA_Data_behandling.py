import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as scio


data_290K = np.loadtxt('~/Documents/Bachelor_Projekt/RITA_PG/RITA-2-PG_as_analyzer_20210312_102231/mccode.dat')

Ei_290K = []
I_290K = []
I_err_1 = []

for i in data_290K:
    Ei_290K.append(i[0]-5)
 
for i in data_290K:
    if i[11] <= 0:
        I_290K.append(1) 
    elif i[11] > 0:
        I_290K.append(i[11]/i[1])

for i in data_290K:
    if i[12] == 0:
        I_err_1.append(i[12]/i[2]) 
    elif i[12] > 0:
        I_err_1.append(i[12]/i[2])    

I_err_290K = []
for n in I_err_1:
    if n < 0.005:
        I_err_290K.append(0.1)
    else:
        I_err_290K.append(n)

""" 10 K data """
data_10K = np.loadtxt('~/Documents/Bachelor_Projekt/RITA_PG/RITA-2-PG_as_analyzer_20210312_102231/mccode.dat')

Ei_10K = []
I_10K = []
I_err_2 = []

for i in data_10K:
    Ei_10K.append(i[0]-5)
 
for i in data_10K:
    if i[11] <= 1:
        I_10K.append(1) 
    elif i[11] > 0:
        I_10K.append(i[11]/i[1])

for i in data_10K:
    if i[12] == 0:
        I_err_2.append(i[12]/i[2]) 
    elif i[12] > 0:
        I_err_2.append(i[12]/i[2])

I_err_10K = []
for n in I_err_2:
    if n < 0.00005:
        I_err_10K.append(0.1)
    else:
        I_err_10K.append(n)

""" Fitting og plotting"""

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

#popt_290K, pcov_290K = scio.curve_fit(gauss_function, Ei_290K, I_290K, sigma=I_err_290K, absolute_sigma=True, p0=[1400000, 5, 0.25])
#popt_10K, pcov_10K = scio.curve_fit(gauss_function, Ei_10K, I_10K, sigma=I_err_10K, absolute_sigma=True, p0=[1400000, 5, 0.25])


#a_290K = popt_290K[0]
#x0_290K = popt_290K[1]
#sigma_290K = popt_290K[2]

#a_10K = popt_10K[0]
#x0_10K = popt_10K[1]
#sigma_10K = popt_10K[2]

energy_range = np.arange(4.5,5.5,0.001)


plt.figure(1)   
plt.errorbar(Ei_290K, I_290K, fmt='b.', capsize=2, label='Simulated data 290K')
plt.errorbar(Ei_10K, I_10K, fmt='r.', capsize=2, label='Simulated data 10K')
#plt.plot(energy_range, gauss_function(energy_range,a_290K,x0_290K,sigma_290K), 'b-', label='fit 290K')
#plt.plot(energy_range, gauss_function(energy_range,a_10K,x0_10K,sigma_10K), 'r-', label='fit 10K')
plt.legend(loc='upper left')
plt.ylabel('I [au]')
plt.xlabel('$\hbar \omega$ [meV]')
plt.yscale('log', nonposy='clip')
plt.show()



