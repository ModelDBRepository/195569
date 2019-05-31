"""

fig2a_run.py

Run simulation for Fig. 2A in Huang S, Hong S, and De Schutter E. 2015. Non-linear leak currents affect mammalian neuron physiology. Front Cell Neurosci 9: 4602-10.

Written by Shiwei Huang, Computational Neuroscience Unit, OIST, 2015

"""


import neuron
from neuron import h
from numpy import *
import pickle
import time
import os

os.makedirs('data')

cvode = h.CVode() ; cvode_is_active = cvode.active(1)
MEMBRANERESIST = 120000. #float
AXIALRESIST    = 140.
MEMBRANECAP    = 0.8
h.celsius      = 34
tstop          = 1

soma = h.Section();
for sec in [soma]:
    sec.diam = 262.6117
    sec.L    = 262.6117
    sec.cm   = MEMBRANECAP
    sec.Ra   = AXIALRESIST


# MEMBRANE PROPERTIES
section_count = 0
for i, sec in enumerate(h.allsec()):
    sec.cm   = MEMBRANECAP
    sec.Ra   = AXIALRESIST

    sec.insert("kcl")
    sec.gsum_kcl = 1/MEMBRANERESIST
    sec.veq_kcl  = -85
    sec.cli = 10
    sec.clo = 130
    sec.ki  = 150
    sec.ko  = 2.5
    sec.nai = 10
    sec.nao = 150

    sec.insert("pas")
    sec.e_pas = -85
    sec.g_pas = 0

    section_count = i

#for multicompartmental model only
spine_dens = 13
spine_area = 1.33
def spinecomp(sec): #thin dendrite corrections
    a = 0
    maxdiam = 0
    for seg in sec.allseg():
	    if (seg.diam > maxdiam):
	        maxdiam = seg.diam
    if (maxdiam <= 3.17):
        a = 0
        for seg in sec.allseg():
            zeta = h.area(seg.x)
            a = a + zeta
        F = (sec.L * spine_area * spine_dens + a)/a
    else:
        F = 1

    return F
#

def get_index(np_array, curr_element):
    for i, val in enumerate(np_array):
        if(val >= curr_element):
            return i-1
            break


def get_params(section):
    """ returns a dictionary containing useful section parameters."""
    vector = {}
    vector['diam'] = section.diam
    vector['L']    = section.L
    vector['Ra']   = section.Ra
    vector['cm']   = section.cm
    vector['ek']   = round(section.ek,2)
    vector['ecl']  = round(section.ecl,2)
    vector['ena']  = round(section.ena,2)
    vector['ko']   = round(section.ko,1)
    vector['ki']   = round(section.ki,1)
    vector['clo']  = round(section.clo,1)
    vector['cli']  = round(section.cli,1)
    vector['nao']  = round(section.nao,1)
    vector['nai']  = round(section.nai,1)

    return vector

def get_iclamp(section, init_mV, holding_i, stepping_i):
    ichold       = h.IClamp(0.5, sec=section)
    ichold.delay = 2000
    ichold.dur   = 100000
    ichold.amp   = holding_i
    icstep       = h.IClamp(0.5, sec=section)
    icstep.delay = 4000
    icstep.dur   = 2000
    icstep.amp   = stepping_i

    hocvec = {}
    for var in 't', 'v', 'icl', 'ik':
        hocvec[var] = h.Vector()
    hocvec['t'].record(h._ref_t)
    hocvec['v'].record(section(0.5)._ref_v)
    hocvec['icl'].record(section(0.5)._ref_icl_kcl)
    hocvec['ik'].record(section(0.5)._ref_ik_kcl)

    h.finitialize(init_mV)
    TSTOP = icstep.delay + icstep.dur + 1000
    neuron.run(TSTOP)

    return hocvec

#generate parameter set
EKCL   = array([-85])
CONC   = linspace(10, 10, 1) #intracellular Cl concentration
IHOLD = linspace(-0.5, 2, 51)
ISTEP = linspace(-0.06, 0.06, 9)
param_set = []
for a in EKCL:
    for b in CONC:
        for c in IHOLD:
            for d in ISTEP:
                line = array([a,b,c,d])
                param_set.append(line)


# PARALLEL RUN
pc = h.ParallelContext();
myid = int(pc.id());
np = int(pc.nhost());
n = len(param_set)

for i in range(myid,n,np):
    print 'myid: ' + str(myid)+':'+str(i); pctime = str(time.time())[5:-3]

    for sec in h.allsec():
        sec.veq_kcl = param_set[i][0]
        sec.e_pas   = param_set[i][0]
        sec.cli     = param_set[i][1]
    h.finitialize(soma.veq_kcl);

    for sec in h.allsec():
        sec.g_pas = 0
        sec.gsum_kcl = 1/MEMBRANERESIST

    #Save useful GHK model parameters
    h.finitialize(soma.veq_kcl);
    vector = get_params(soma)
    vector['gsum_kcl']   = soma.gsum_kcl
    vector['Pratio_kcl'] = round(soma.Pratio_kcl, 3)
    vector['veq_kcl']    = soma.veq_kcl
    vector['e']          = soma.veq_kcl
    vector['IHOLD'] = param_set[i][2]
    vector['ISTEP'] = param_set[i][3]

    hvec_d = get_iclamp(soma, soma.veq_kcl, vector['IHOLD'], vector['ISTEP'])
    for key in hvec_d:
	    vector[key] = array( hvec_d[key] );

    #repeat simulation with ohmic leak current
    for sec in h.allsec():
        sec.gsum_kcl = 0
        sec.g_pas = 1/MEMBRANERESIST

    # save useful ohmic model paramters
    h.finitialize(soma.e_pas)
    vector2  = get_params(soma)
    vector2['g_pas']  = soma.g_pas
    vector2['e_pas']  = round(soma.e_pas, 2)
    vector2['e']      = round(soma.e_pas, 2)
    vector2['IHOLD'] = vector['IHOLD']
    vector2['ISTEP'] = vector['ISTEP']

    hvec_d2 = get_iclamp(soma, soma.e_pas, vector2['IHOLD'], vector2['ISTEP'])
    for key in hvec_d2:
        vector2[key] = array( hvec_d[key] )

    fname='Pr_%g_ecl_%g_CONC_%g_EL_%g_IHOLD_%g_ISTEP_%g.pkl' \
	  % (vector['Pratio_kcl'], vector['ecl'], vector['cli'], vector['veq_kcl'], \
	     vector['IHOLD'], vector['ISTEP'])
    fname = os.path.join('data', fname)
    pklfile  = open(fname, 'wb')
    pickle.dump(vector, pklfile)
    pickle.dump(vector2, pklfile)
    pklfile.close()

    vector  = None
    vector2 = None
pc.done();
