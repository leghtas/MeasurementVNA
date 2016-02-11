# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 14:09:29 2015

@author: Tennessee
"""

from measurementVna import MeasurementVna
from mclient import instruments
import numpy as np
import matplotlib.pyplot as plt
#plt.close('all')

#vna = instruments.create('vna1', 'Agilent_E5071C', address='GPIB0::29')
vna = instruments.create('vna1', 'Agilent_E5071C', address='GPIB::29')
ag1 = instruments.create('ag1', 'Agilent_N5183A', address='GPIB::5')

#ag1.set_rf_on(True)

vna_exp= MeasurementVna(vna, 'Readout_Spec', gen=ag1,
#                        vnafreqstartstop = (7.05e9, 7.1457e9), 
#                        vnapowers=np.linspace(-45, -20, 1), 
                        vnaavg=50, vnapts=201, ifbw=1000,
                        genpowers=np.asarray([15]),
                        genfreqs=np.linspace(9.5e9, 9.6e9, 21),
                        fitdata=False, f_0_guess=7.175e9, kc_guess=300e3, ki_guess=10e3)
vna_exp.measure()

#ag1.set_rf_on(False)