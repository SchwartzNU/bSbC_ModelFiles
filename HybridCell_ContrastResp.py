# imports
import os
os.chdir('C:/Users/sophi/Documents/Schwartz/NEURON/bSbC Paper Model/IndiraModel2/')
import csv
import numpy
import scipy.io

from neuron import h
from neuron.units import ms, mV
from neuron import gui
h.load_file('stdrun.hoc')

import matplotlib.pyplot as plt

from HybridCell import HybridCell

temp = 32
h.celsius = temp

vLeak = -50
v_init = -61.5
NaVRatio = [0.4, 0]
# Na Density has the largest effect on block
NaDensity =  [0.003671929, 0.002279877] #0.002975903
KDensity = 0.004841545 #[0.003726023, 0.005957068] #0.004841545
SomaDiam = 17.5
DendLen = 513.0843283
HillLen = 20.75
AISLen = [20.5, 14.5]
Multiplier = 30

h.dt = 1/10 #1/ms


OFFsACell = HybridCell(vLeak, NaVRatio[0], Multiplier, NaDensity[0], KDensity, SomaDiam, DendLen, HillLen, AISLen[0])
bSbCCell = HybridCell(vLeak, NaVRatio[1], Multiplier, NaDensity[1], KDensity, SomaDiam, DendLen, HillLen, AISLen[1])

# measure IR
ImpA = h.Impedance()
ImpA.loc(0.5, sec=OFFsACell.soma)
ImpA.compute(0)
IR_val_A = ImpA.input(0.5, sec=OFFsACell.soma)
print(IR_val_A)
ImpB = h.Impedance()
ImpB.loc(0.5, sec=bSbCCell.soma)
ImpB.compute(0)
IR_val_B = ImpB.input(0.5, sec=bSbCCell.soma)
print(IR_val_B)


# load mat files of conductances
filepath = 'C:/Users/sophi/Documents/Schwartz/analysis/DynamicClampConductances/'
Alpha_n100_exc = scipy.io.loadmat(filepath+'Sophia_Alpha_cm100_Exc'+'.mat')
Alpha_n100_inh = scipy.io.loadmat(filepath+'Sophia_Alpha_cm100_Inh'+'.mat')
bSbC_n100_exc = scipy.io.loadmat(filepath+'Sophia_Bursty_cm100_Exc'+'.mat')
bSbC_n100_inh = scipy.io.loadmat(filepath+'Sophia_Bursty_cm100_Inh'+'.mat')

Factor =0.6 #  0.48 # average multiplier for DC recordings 0.4786
A_Exc = h.Vector(1000/Alpha_n100_exc['conductances'][0]) * Factor # needs to be in MOhms
A_Inh = h.Vector(1000/Alpha_n100_inh['conductances'][0]) * Factor
B_Exc = h.Vector(1000/bSbC_n100_exc['conductances'][0]) * Factor
B_Inh = h.Vector(1000/bSbC_n100_inh['conductances'][0]) * Factor

# add in the conductances measured for dynamic clamp and input them
OFFsA_E = h.SEClamp(OFFsACell.soma(0.4))
OFFsA_E.dur1 = 1e9
OFFsA_E.amp1 = 0 # set to reversal potential
OFFsA_I = h.SEClamp(OFFsACell.soma(0.6))
OFFsA_I.dur1 = 1e9
OFFsA_I.amp1 = -70

bSbC_E = h.SEClamp(bSbCCell.soma(0.4))
bSbC_E.dur1 = 1e9
bSbC_E.amp1 = 0 # set to reversal potential
bSbC_I = h.SEClamp(bSbCCell.soma(0.6))
bSbC_I.dur1 = 1e9
bSbC_I.amp1 = -70

# run conductances
A_Exc.play(OFFsA_E._ref_rs, h.dt)
A_Inh.play(OFFsA_I._ref_rs, h.dt)
B_Exc.play(bSbC_E._ref_rs, h.dt)
B_Inh.play(bSbC_I._ref_rs, h.dt)

TimeVec = numpy.linspace(0, 2.5, num=25001, endpoint=True) # sample rate 10000
SomaV_A = h.Vector().record(OFFsACell.soma(0.5)._ref_v)
SomaV_B = h.Vector().record(bSbCCell.soma(0.5)._ref_v)
h.finitialize(-60 * mV)
h.continuerun(2500 * ms)

plt.figure(1)
plt.plot(TimeVec, list(SomaV_A))
plt.xlabel('t (s)')
plt.ylabel('v (mV)')
plt.title('OFFsA, bSbC Conductances')
plt.figure(2)
plt.plot(TimeVec, list(SomaV_B))
plt.xlabel('t (s)')
plt.ylabel('v (mV)')
plt.title('bSbC, OFFsA Conductances')


scipy.io.savemat('Voltages_ContrastResponse_matched.mat', dict(OFFsA = numpy.array(SomaV_A), bSbC = numpy.array(SomaV_B), TimeVec = numpy.array(TimeVec)))
print('done')


### try a different version
#IClamp_A = h.IClamp(OFFsACell.soma(0.5))
#IClamp_A.delay = 0
#IClamp_A.dur = 1e9
#IClamp_B = h.IClamp(bSbCCell.soma(0.5))
#IClamp_B.delay = 0
#IClamp_B.dur = 1e9
#
#NSteps = 25000
#TimeVec = h.Vector().record(h._ref_t)
#SomaV_A = h.Vector().record(OFFsACell.soma(0.5)._ref_v)
#SomaV_B = h.Vector().record(bSbCCell.soma(0.5)._ref_v)
#
#ExcR = 0 # reversal potential of excitatory synapses
#InhR = -60 # reversal potential of inhibitory synapses
#Factor = 0.3
#
#h.finitialize(v_init * mV)
#for i in range(0, NSteps):
#    # calculate the amount of current to inject based on the conductances
#    # I = g(v-e)
#    v_A = SomaV_A[i]
#    excG_A = A_Exc[i] * Factor
#    inhG_A = A_Inh[i] * Factor
#    I_A = (excG_A * (v_A - ExcR) + inhG_A * (v_A - InhR)) / 1000 # needs to be in nA
#    IClamp_A.amp = -I_A
#    v_B = SomaV_B[i]
#    excG_B = B_Exc[i] * Factor
#    inhG_B = B_Inh[i] * Factor
#    I_B = (excG_B * (v_B - ExcR) + inhG_B * (v_B - InhR)) / 1000 # needs to be in nA
#    IClamp_B.amp = -I_B
#    h.fadvance()
#
#plt.figure(1)
#plt.plot(list(TimeVec), list(SomaV_A))
#plt.xlabel('t (ms)')
#plt.ylabel('v (mV)')
#
#plt.figure(2)
#plt.plot(list(TimeVec), list(SomaV_B))
#plt.xlabel('t (ms)')
#plt.ylabel('v (mV)')
