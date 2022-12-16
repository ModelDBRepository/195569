"""

fig2a_show.py

Make Fig. 2A in Huang S, Hong S, and De Schutter E. 2015. Non-linear leak currents affect mammalian neuron physiology. Front Cell Neurosci 9: 4602-10, based on simulated data created by fig2a_run.py.

Written by Shiwei Huang, Computational Neuroscience Unit, OIST (shiwei.huang@hotmail.com.au)
Updated 2015

"""


#import matplotlib
#matplotlib.use('SVG')
import matplotlib.pyplot as plt
import os
import pickle
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def get_index_of_prev_element(np_array, curr_element):
    for i, val in enumerate(np_array):
        if(val >= curr_element):
            return i-1
            break

allfiles  = os.listdir('data')
pklfiles  = [os.path.join('data', x) for x in allfiles if( 'Pr_1.341_ecl_-67.9_CONC_10_EL_-85_IHOLD_0.3_ISTEP_' in x and ('0.015' in x) )]
pklfiles1 = [os.path.join('data', x) for x in allfiles if( ('Pr_1.341_ecl_-67.9_CONC_10_EL_-85_IHOLD_0_ISTEP_' in x) and ('0.015' in x) )]
pklfiles2 = [os.path.join('data', x) for x in allfiles if( 'Pr_1.341_ecl_-67.9_CONC_10_EL_-85_IHOLD_-0.3_ISTEP_' in x and ('0.015' in x) )]

fig = plt.figure()
ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)
ax6 = fig.add_subplot(326)

def load_and_plottraces(dataset, ax1, ax2):
    for n, fname in enumerate(dataset):
        f = open(fname, 'rb')
        gdat = pickle.load(f)
        odat = pickle.load(f)

        v0 = gdat['v']
        t  = gdat['t']
        t0 = get_index_of_prev_element(t, 3999)
        v  = gdat['v'] - gdat['v'][t0]

        vo0 = odat['v']
        to  = odat['t']
        to0 = get_index_of_prev_element(to, 3999)
        vo  = odat['v'] - odat['v'][to0]

        if (gdat['ISTEP'] < 0):
            v = abs(v)
            ax2.plot(t, v, color='orange', linestyle='dotted', linewidth=3)
            vo = abs(vo)
            ax1.plot(to, vo, color='blue', linestyle='dotted', linewidth=3)

        else:
            ax2.plot(t, v, 'orange')
            ax1.plot(to, vo, 'blue')

        for l in ax1, ax2:
            l.set_ylim([0, 3])
            l.set_xlim([3500, 7000])

        f.close()

load_and_plottraces(pklfiles, ax1, ax2)
load_and_plottraces(pklfiles1, ax3, ax4)
load_and_plottraces(pklfiles2, ax5, ax6)

for ax in (ax5, ax6): ax.set_xlabel('Time (ms)')

plt.show()
