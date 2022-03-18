# imports
import os
os.chdir('C:/Users/Documents/Schwartz/NEURON/bSbC Paper Model/')
import csv
import numpy
import scipy.io

from neuron import h
from neuron.units import ms, mV
from neuron import gui

h.load_file('stdrun.hoc')
temp = 32
h.celsius = temp
dt = 0.01
h.dt = dt

## Notes:
## using mod files from (Hu et al. 2009) written by Zach Mainen, Salk Institute, 1994, zach@salk.edu
## some values from Fohlmeister 1997 multiplied by 1E-9

class HybridCell:
    def __init__(self, v_init, NaVRatio, Multiplier, NaDensity, KDensity, SomaDiam, DendLen, HillLen, AISLen):
        self._setup_morphology(SomaDiam, DendLen, HillLen, AISLen)
        self._setup_biophysics(v_init, NaVRatio, Multiplier, NaDensity, KDensity)
    
    def _setup_morphology(self, SomaDiam, DendLen, HillLen, AISLen):
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.hill = h.Section(name='hill', cell=self) # hillock
        self.AIS = h.Section(name='AIS', cell=self)
        self.axon = h.Section(name='axon', cell=self)
        self.all = [self.soma, self.dend, self.hill, self.AIS, self.axon]
        self.dend.connect(self.soma(0.0))
        self.hill.connect(self.soma(1.0))
        self.AIS.connect(self.hill(1.0))
        self.axon.connect(self.AIS(1.0))
        # give all fifteen segments
        self.soma.nseg = 15
        self.dend.nseg = 15
        self.hill.nseg = 15
        self.AIS.nseg = 15
        self.axon.nseg = 15
        # define physical params (microns)
        self.soma.L = self.soma.diam = SomaDiam # cylinder w/ SA with the sphere SA Not CORRECT
        self.dend.L = DendLen
        self.dend.diam = 1
        self.hill.L = HillLen
        self.hill.diam = 1
        self.AIS.L = AISLen
        self.AIS.diam = 1
        self.axon.L = 1000 # arbitrarily set
        self.axon.diam = 1
        
    def _setup_biophysics(self, v_init, NaVRatio, Multiplier, NaDensity, KDensity):
        for sec in self.all:
            sec.cm = 1 # Membrane capacitance in micro Farads / cm^2
            sec.Ra = 200 # 100 # ohm-cm Schachter 2010, Abbas 2013
        print(v_init)
        
        Calc_Rm = 2000 # 9500 - 70000 Ohm cm2, Freed Sterling 1992
        self.soma.insert('pas')
        for seg in self.soma:
            seg.pas.g = 1/Calc_Rm
            seg.pas.e = v_init
        self.dend.insert('pas')
        for seg in self.dend:
            seg.pas.g = 1/Calc_Rm
            seg.pas.e = v_init
        self.hill.insert('pas')
        for seg in self.hill:
            seg.pas.g = 1/Calc_Rm
            seg.pas.e = v_init
        self.AIS.insert('pas')
        for seg in self.AIS:
            seg.pas.g = 1/Calc_Rm
            seg.pas.e = v_init
        self.axon.insert('pas')
        for seg in self.axon:
            seg.pas.g = 1/Calc_Rm
            seg.pas.e = v_init
        
        # values from Fried Frontiers, supplemental table 1
        na16Ratio = NaVRatio
        na12Ratio = 1-na16Ratio
        AISMulti = Multiplier
        # input both at nS/cm^2
        NaDens = NaDensity/1000*100000000 # 547.9397037 # (pS/um2)
        KDens = KDensity/1000*100000000 # (pS/um2)

        ## soma
        self.soma.insert('kv')
        for seg in self.soma:
            seg.kv.gbar = KDens
        self.soma.insert('na12')
        self.soma.insert('na16')
        for seg in self.soma:
            seg.na12.gbar = NaDens * na12Ratio
            seg.na16.gbar = NaDens * na16Ratio
            
        ## dendrites
        self.dend.insert('kv')
        for seg in self.dend:
            seg.kv.gbar = KDens
        self.dend.insert('na12')
        self.dend.insert('na16')
        for seg in self.dend:
            seg.na12.gbar = NaDens * na12Ratio
            seg.na16.gbar = NaDens * na16Ratio
        
        ## hillock
        self.hill.insert('kv')
        for seg in self.hill:
            seg.kv.gbar = KDens
        self.hill.insert('na12')
        self.hill.insert('na16')
        for seg in self.hill:
            seg.na12.gbar = NaDens * na12Ratio
            seg.na16.gbar = NaDens * na16Ratio
        
        ## AIS
        self.AIS.insert('kv')
        for seg in self.AIS:
            seg.kv.gbar = AISMulti * KDens
        # add specific na channels
        self.AIS.insert('na12')
        self.AIS.insert('na16')
        for seg in self.AIS:
            seg.na12.gbar = AISMulti * NaDens * na12Ratio
            seg.na16.gbar = AISMulti * NaDens * na16Ratio
        
        ## axon
        self.axon.insert('kv')
        for seg in self.axon:
            seg.kv.gbar = KDens
        self.axon.insert('na12')
        self.axon.insert('na16')
        for seg in self.axon:
            seg.na12.gbar = NaDens * na12Ratio
            seg.na16.gbar = NaDens * na16Ratio
        
        Ena = 30
        Ek = -90
        
        self.dend.ena = Ena
        self.soma.ena = Ena
        self.hill.ena = Ena
        self.AIS.ena = Ena
        self.axon.ena = Ena
        
        self.dend.ek = Ek
        self.soma.ek = Ek
        self.hill.ek = Ek
        self.AIS.ek = Ek
        self.axon.ek = Ek
    
    def __repr__(self):
        return 'AISCell[{}]'
        

