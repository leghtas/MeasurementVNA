# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 11:31:08 2015
@author: Tennessee
"""
import time
import h5py
import config
import os
import analysis_functions as af
import numpy as np
import matplotlib.pyplot as plt
import datetime

class MeasurementVna(object):
    '''
    class for measurement using VNA and generators,
    CW spec,
    sweeping powers and freqs
    '''
    def __init__(self, vna, exptname, vnaavg=None, ifbw=None, vnapts=None,
                 vnapowers=None, vnafreqstartstop=(None, None), gen=None, genpowers=[None],
                 genfreqs=[None], fitdata=False,
                 f_0_guess=None, kc_guess=None, ki_guess=None,
                 a_in_guess=None, T_guess=0e-9):
        self.vna = vna
        self.exptname = exptname
        self.vnaavg = vnaavg
        self.ifbw = ifbw
        self.vnapts = vnapts
        self.vnapowers = vnapowers
        self.vnafreqstartstop = vnafreqstartstop
        self.gen = gen
        self.genpowers=genpowers
        self.genfreqs=genfreqs
        self.fitdata=fitdata
        self.f_0_guess=f_0_guess
        self.kc_guess=kc_guess
        self.ki_guess=ki_guess
        self.a_in_guess=a_in_guess
        self.T_guess=T_guess
        self.data_filename = None
        self.image_filename = None
        self.sweepparam = 'None'
        self.sweepbounds = ['None', 'None']

        self.set_vna_params()


    #def get_vna_params():

    def set_vna_params(self):
        vna=self.vna
        gen=self.gen
        if self.vnaavg is None:
            self.vnaavg=vna.get_average_factor()
        else:
            vna.set_average_factor(self.vnaavg)
        if self.ifbw is None:
            self.ifbw=vna.get_if_bandwidth()
        else:
            vna.set_if_bandwidth(self.ifbw)
        if self.vnapts is None:
            self.vnapts=vna.get_points()
        else:
            vna.set_points(self.vnapts)
        if self.vnapowers is None:
            self.vnapowers=[vna.get_power()]
        else: vna.set_power(self.vnapowers[0])
        if self.vnafreqstartstop[0] is None or self.vnafreqstartstop[1] is None:
            self.vnafreqstartstop=vna.get_start_freq(), vna.get_stop_freq()
        else:
            vna.set_start_freq(self.vnafreqstartstop[0])
            vna.set_stop_freq(self.vnafreqstartstop[1])
        if gen is not None:
           if self.genpowers[0] is None:
               self.genpowers=[gen.get_power()]
           else:
               gen.set_power(self.genpowers[0])
           if self.genfreqs[0] is None:
               self.genfreqs=[gen.get_frequency()]
           else:
               gen.set_frequency(self.genfreqs[0])
               
        self.sweepTime = vna.get_sweep_time()

    def set_filenames(self):
        date = time.strftime('%Y%m%d', time.localtime())
        rootDir = config.datadir+'/VNA/'
        ts = time.localtime()
        tstr = time.strftime('%Y%m%d\%H%M%S', ts)
        imageFile = os.path.join(rootDir, 'images\%s_%s'%(tstr, self.exptname))
        imageDir = os.path.split(imageFile)[0]
        if not os.path.isdir(imageDir): os.makedirs(imageDir)
        if not os.path.isdir(rootDir + '/' + self.exptname): os.makedirs(rootDir + '/' + self.exptname)
        self.data_filename = rootDir + '/' + self.exptname + '/' + date + '.h5'
        tree = h5py.File(self.data_filename)
        if tree.keys():
            self.expNum = str(max(map(int, tree.keys())) + 1)
        else:
            self.expNum = '0'

        imageFile += '_vnaPow_' + str(self.vna.get_power())
        if self.gen is not None: imageFile += '_genPow_' + str(self.gen.get_power()) + '_genFreq_' + str(self.gen.get_frequency())
                
        self.image_filename = imageFile + '.png'
        self.image_filename2D = imageFile + '_2D.png'

    def save_data(self, data):
        '''
        exptm paramaters after setting None to default values
        '''
        h5Data = {
        'vnaaddress': self.vna.get_address(),
        'exptname': self.exptname,
        'vnaavg': self.vnaavg,
        'ifbw': self.ifbw,
        'vnapts': self.vnapts,
        'vnapower': self.vna.get_power(),
        'vnafreqs': self.vnafreqs,
        'data': data,
        'timestamp': time.strftime('%Y%m%d-%H%M%S', time.localtime()),
        'sweepparam': self.sweepparam,
        'sweepbounds': self.sweepbounds,
        }
        if self.fitdata:
            h5Data['fitdata']=self.fitdata
            h5Data['fitparams']=self.fitparams
        if self.gen is not None:
            h5Data['gen']=self.gen.get_name()
            h5Data['genpower']=self.gen.get_power()
            h5Data['genfreq']=self.gen.get_frequency()

        tree = h5py.File(self.data_filename)
        expTree = tree.create_group(self.expNum)

        for name, data in h5Data.items():
            expTree[name] = data
        tree.close()

    def saveimage(self, fig, image_filename):
        fig.savefig(image_filename)
            
    def plot1D(self, data, fig):
        vnafreqs = self.vnafreqs
        
        plt.figure(fig.number)
        plt.clf()
        fig.suptitle(self.image_filename+'\n'+self.data_filename+'/'+self.expNum)

        ax1 = fig.add_subplot(221)
        ax1.plot(np.real(data), np.imag(data))
        ax1.set_xlabel('Re(a)')
        ax1.set_ylabel('Im(a)')
        ax1.ticklabel_format(useOffset=False)
        ax2 = fig.add_subplot(222)
        ax2.plot(vnafreqs*1e-9, np.abs(data))
        ax2.set_xlabel('VNA freq (GHz)')
        ax2.set_ylabel('|a|')
        ax2.ticklabel_format(useOffset=False)

        ax3 = fig.add_subplot(223)
        ax3.plot(vnafreqs*1e-9, np.angle(data, deg=True))
        ax3.set_xlabel('VNA freq (GHz)')
        ax3.set_ylabel('angle(a)')
        ax3.ticklabel_format(useOffset=False)
        ax4 = fig.add_subplot(224)
        ax4.plot(vnafreqs*1e-9, 20 * np.log(np.abs(data)) / np.log(10))
        ax4.set_xlabel('VNA freq (GHz)')
        ax4.set_ylabel('20 log|a|')
        ax4.ticklabel_format(useOffset=False)

        if self.fitdata and len(self.fitparams)>0:
            popt = self.fitparams
            data_fit = af.complex_a_out(vnafreqs, *popt)
            ax1.plot(np.real(data_fit),np.imag(data_fit),'--')
            ax2.plot(vnafreqs*1e-9,np.abs(data_fit),'--')
            ax3.plot(vnafreqs*1e-9,np.angle(data_fit, deg=True),'--')
            ax4.plot(vnafreqs*1e-9, 20 * np.log(np.abs(data_fit)) / np.log(10))

        plt.show()
        self.saveimage(fig, self.image_filename)

        plt.pause(0.1)

    def plot2D(self, dataArray, fig):
        #vnafreqs=self.vnafreqs
        plt.figure(fig.number)
        plt.clf()

        fig.suptitle(self.image_filename2D+'\n'+self.data_filename+'/'+self.expNum)
        ax1 = fig.add_subplot(211)
        ax1.imshow(np.abs(dataArray), interpolation='nearest', origin='lower', aspect='auto',
                   extent=(self.vnafreqstartstop[0]*1e-9, self.vnafreqstartstop[1]*1e-9, self.sweepbounds[0], self.sweepbounds[1]))
        ax1.set_title('|a|')
        ax1.set_ylabel(self.sweepparam)
        ax1.set_xlabel('VNA frequency (GHz)')
        ax1.ticklabel_format(useOffset=False)

        ax2 = fig.add_subplot(212)
        ax2.imshow(np.angle(dataArray, deg=True), interpolation='nearest', origin='lower', aspect='auto',
                   extent=(self.vnafreqstartstop[0]*1e-9, self.vnafreqstartstop[1]*1e-9, self.sweepbounds[0], self.sweepbounds[1]))
        ax2.set_title('Phase')
        ax2.set_ylabel(self.sweepparam)
        ax2.set_xlabel('VNA frequency (GHz)')
        ax2.ticklabel_format(useOffset=False)
        plt.show()
        plt.pause(0.1)

    def plotfits(self):
        f0Vec = [self.fitparamsVec[k][0] for k in range(np.size(self.sweepaxis))]
        kcVec = [self.fitparamsVec[k][1] for k in range(np.size(self.sweepaxis))]
        kiVec = [self.fitparamsVec[k][2] for k in range(np.size(self.sweepaxis))]
        
        fig = plt.figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax1.plot(self.sweepaxis, np.asarray(f0Vec)*1e-9)
        ax1.set_ylabel('f0 (GHz)')
        ax1.set_xlabel(self.sweepparam)
        ax2.plot(self.sweepaxis, np.asarray(kcVec)*1e-6)
        ax2.set_ylabel('kappaC/2pi (MHz)')
        ax2.set_xlabel(self.sweepparam)
        ax3.plot(self.sweepaxis, np.asarray(kiVec)*1e-6)
        ax3.set_ylabel('kappaI/2pi (MHz)')
        ax3.set_xlabel(self.sweepparam)
        ax1.ticklabel_format(useOffset=False)
        ax2.ticklabel_format(useOffset=False)
        ax3.ticklabel_format(useOffset=False)        
        
        plt.show()
        
    def get_vna_trace(self):
        vna=self.vna
        vna.do_enable_averaging()
        vna.set_averaging_trigger(False)
        time.sleep(self.sweepTime*self.vnaavg*1.1)

        ## Data acquisition
        self.vnafreqs = vna.do_get_xaxis()
        data_re, data_im = vna.do_get_yaxes()
        data = data_re+1j*data_im
        if self.fitdata:
            try:
                self.fitparams = af.fit_complex_a_out(self.vnafreqs, data, f_0=self.f_0_guess, kc=self.kc_guess, ki=self.ki_guess, a_in=self.a_in_guess, T=self.T_guess)
                if self.plotsweep: self.fitparamsVec.append(self.fitparams)
            except SyntaxError:
                self.fitparams = []
        return data


    def measure(self):
        self.vna.set_format('POL')
        vnapowers = self.vnapowers
        genfreqs = self.genfreqs
        genpowers = self.genpowers

        print 'This is going to take (h:m:s) '  + str(datetime.timedelta(seconds=max(int(self.vnaavg*np.size(vnapowers)*np.size(genfreqs)*np.size(genpowers)*self.sweepTime),1)))

        if self.gen is not None: 
            self.gen.set_rf_on(True)

        fig1 = plt.figure(figsize=(11.8 , 11.85))

        self.plotsweep = np.size(vnapowers)>1 or np.size(genfreqs)>1 or np.size(genpowers)>1
        plot_vnapows = np.size(vnapowers)>1
        plot_genfreqs = np.size(vnapowers)==1 and np.size(genfreqs)>1
        plot_genpows =  np.size(vnapowers)==1 and np.size(genfreqs)==1 and np.size(genpowers)>1

        if self.plotsweep:
            fig2 = plt.figure(figsize=(11.8 , 11.85))
            if self.fitdata:
                self.fitparamsVec = []

        if plot_genpows:
            dataArray = []
            self.sweepparam = 'genpower (dBm)'
            self.sweepaxis = genpowers
        for genpower in genpowers:
            if plot_genfreqs :
                dataArray = []
                self.sweepparam = 'genfreq (Hz)'
                self.sweepaxis = genfreqs
            if self.gen is not None and genpower is not None: self.gen.set_power(genpower)

            for genfreq in genfreqs:
                if plot_vnapows:
                    dataArray = []
                    self.sweepparam = 'vnapower (dBm)'
                    self.sweepaxis = vnapowers
                if self.gen is not None and genfreq is not None : self.gen.set_frequency(genfreq)

                for vnapower in vnapowers:
                    self.vna.set_power(vnapower)
                    data = self.get_vna_trace() # fits data
                    self.set_filenames()
                    self.save_data(data)
                    self.plot1D(data, fig1)

                    if self.plotsweep: dataArray.append(data)
                    if plot_vnapows:
                        self.sweepbounds = [self.vnapowers[0], vnapower]
                        self.plot2D(dataArray, fig2)
                if plot_vnapows: self.saveimage(fig2, self.image_filename2D)
                if plot_genfreqs:
                    self.sweepbounds = [self.genfreqs[0], genfreq]
                    self.plot2D(dataArray, fig2)
            if plot_genfreqs: self.saveimage(fig2, self.image_filename2D)
            if plot_genpows:
                self.sweepbounds = [self.genpowers[0], genpower]
                self.plot2D(dataArray, fig2)
        if plot_genpows: self.saveimage(fig2, self.image_filename2D)
        if self.fitdata and self.plotsweep:
            self.plotfits()
