import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
from iminuit import Minuit
from ExternalFunctions import Chi2Regression
from scipy import stats


import matplotlib as mpl
mpl.rcParams['figure.figsize']   = (18,10)
mpl.rcParams['font.size']        = 25 # standard er 45
mpl.rcParams['lines.color']      = 'r'
mpl.rcParams['lines.markersize'] = 15
plt.rcParams['figure.constrained_layout.use'] = True


data_290K = np.loadtxt('290K_Escan_frimar19y21_2/mccode.dat')

col_name = ['Ef', 'psd_monitor_mono_I', 'psd_monitor_mono_ERR', 'psd_monitor_after_mono_I', 'psd_monitor_after_mono_ERR', 'psd_monitor_sample_I', 'psd_monitor_sample_ERR', 'e_monitor_sample_I', 'e_monitor_sample_ERR', 'psd_monitor_bef_analyzer_I', 'psd_monitor_bef_analyzer_ERR', 'psd_monitor2_I', 'psd_monitor2_ERR']

frame_290 = pd.DataFrame(data_290K, columns=col_name)

plt.plot(frame_290['Ef'], frame_290['psd_monitor2_I']/frame_290['psd_monitor_sample_I'],'.')
plt.show()