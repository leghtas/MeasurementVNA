# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 16:38:20 2013

@author: Zaki
"""

from pylab import *

import pylab
from scipy.optimize import curve_fit, leastsq
from math import *
import operator
import scipy.special
import matplotlib.pyplot as plt
import scipy.stats
import time
import numpy as np
from scipy.interpolate import splrep, sproot, splev


font = {'family' : 'arial',
        'weight' : 'normal',
        'size'   : 22}


# Constants
def constants():
    fs=24
    fsTicks=22
    Vg4     =0.0016498356460426227
    Ve4     = 0.01410316219620108
    Vg8ns   = -0.001479
    Ve8ns   = -0.000434
    Vg200   = -0.001492
    Ve200   = -0.000551
    Vg400   = -0.0003583683581309346
    Ve400   = 0.00030092166965822341
    VoltageCal= 4.0
    return (Vg4, Ve4, fs, fsTicks)

f_storage   = 7.578615
f_readout   = 7.15
f_qubit     = 4.900619
f_pump_1    = 7.26172
f_pump_2    = 4.999
f_pump_3    = 7.694

#
def lorentzian(x,x0,y0,A,B):
    numerator =  A
    denominator = ( x - x0 )**2 + B
    y = y0 + (numerator/denominator)
    return y

def reso_params(popt):
    (x0,y0,A,B) = popt

    return "f0: %e\nbw: %e" % (x0,2*sqrt(B))

def fit_lorentzian(x,y):
    index, ymax = max(enumerate(y), key=operator.itemgetter(1))
    x0 = x[index]
    y0 = min(y)

    spline = splrep(x, y-(y0+ymax)/2)
    roots = sproot(spline)
    whm = roots[-1]-roots[0]

   # y0=-0.00072
    #x0=7.68377*1e9
    #whm=0.1*1e6
    B = whm*whm/4.0
    A = B*(ymax-y0)
    #print x0,y0,B,A
    #x0=7.68374*1e9;y0=0.0019;B=((100*1e3)**2)/4;A=B*(max(y)-min(y))


    popt, pcov  = curve_fit(lorentzian,x,y,(x0,y0,A,B))
    #popt=[x0,y0,A,B]

    return popt



def fit_resonance(peak,bg,show_plot=False):
    freq = peak[0]
    gammadB = peak[1]-bg[1]
    gamma = pow(10.0,0.1*gammadB)
    print fit_lorentzian(freq,gamma,show_plot)

def estimate_phase(spec):
    n=len(spec)
    m=len(spec[0])
    x=real(spec)
    y=imag(spec)


def load_file_VNA(file):
    a = loadtxt(file)
    points_per_record = len(a[0])
    recs = int(len(a)/3)
    return reshape(a, (recs, 3, points_per_record))

def plot_all(data):
    for sweep in data:
        fit_lorentzian(data[0],sweep)

def load_file(filename):
    filename_re=filename+"_re_end.dat"
    filename_im=filename+"_im_end.dat"
    return loadtxt(filename_re) + 1j* loadtxt(filename_im)

def loadspec(filename,fstart,fstop,ttranspose=False,rrotate=True):
    spec=load_file(filename)
    if ttranspose:
        spec=transpose(spec)
    if rrotate:
        spec=IQrotate(spec)
    freq = linspace(fstart,fstop,size(spec))
    return (freq,spec)

def loadautorabi(filename):
    data=loadtxt(filename)
    n=floor(size(data)/2)
    data=data[0:n]+1j*data[n:2*n]
    theta=IQrotate_hist(data)
    data=np.exp(1j*theta)*data
    return data

def loadT1T2E(filename,dt,numsteps):
    y=loadautorabi(filename)
    x=np.linspace(0,numsteps*dt,numsteps)
    return (x,y)

def loadautorabinoIQrotate(filename):
    data=loadtxt(filename,dtype='float32')
    n=floor(size(data)/2)
    data=data[0:n]+1j*data[n:2*n]
    #data=IQrotate(data)
    return data

def cosine(x,A,omega,y0,x0):
    return y0+A*cos(omega*(x-x0))

def expcosine(t,A,Tdet,y0,t0,T2):
    y=y0+A*exp(-t/T2)*cos(2*pi*(t-t0)/Tdet)
    return y

def arctan_fit(x,y0,x0,A,B,C):
    y=y0+A*np.arctan((x-x0)/B)+C*x
    return y

def gausscosine(t,A,Tdet,y0,t0,T2):
    y=y0+A*exp(-t**2/(2*T2**2))*cos(2*pi*(t-t0)/Tdet)
    return y

def Previval(t,A,chi,y0,t0,T2,nbar):
    y=exp(nbar*(cos(chi*(t-t0))-1))
    y=y*cos(nbar*sin(chi*(t-t0)))
    y=(1-y)
    y=A*exp(-t/T2)*y
    y=y+y0
    return y

def Previval_detune(t,t0,A,chi,y0,T2,nbar, detune):
    y=exp(nbar*(cos(chi*(t-t0))-1))
    y=y*cos(nbar*sin(chi*(t-t0)) + detune*(t - t0))
    y=(1-y)
    y=A*exp(-t/T2)*y
    y=y+y0
    return y

# def fit_revival_detune(t,y,chi,T2,nbar, detune):
#     A=max(y)-min(y)
#     t0=0
#     y0=mean(y)
#     #popt, pcov = curve_fit(Previval,t,y,(A,chi,y0,t0,T2,nbar))
#     popt, pcov = curve_fit(lambda t,A,chi,y0,nbar:Previval(t,A,chi,y0,t0,T2,nbar, detune), t, y ,(A,chi,y0,nbar))
#     return popt

def fit_revival_detune(t,y,chi,T2,nbar, detune):
    A=max(y)-min(y)
    t0=0
    y0=mean(y)
    #popt, pcov = curve_fit(Previval,t,y,(A,chi,y0,t0,T2,nbar))
    popt, pcov = curve_fit(lambda t,t0,A,chi,y0,nbar:Previval_detune(t,t0,A,chi,y0,T2,nbar, detune), t, y ,(t0,A,chi,y0,nbar))
    return popt

def fit_revival(t,y,chi,T2,nbar):
    A=max(y)-min(y)
    t0=0
    y0=mean(y)
    #popt, pcov = curve_fit(Previval,t,y,(A,chi,y0,t0,T2,nbar))
    popt, pcov = curve_fit(lambda t,A,chi,y0,nbar:Previval(t,A,chi,y0,t0,T2,nbar), t, y ,(A,chi,y0,nbar))
    return popt

def fit_arctan(x,y,x0,B):
    A=(-max(y)+min(y))/np.pi
    y0=mean(y)
    C=-2*np.pi*20/(3*1e8)
    popt, pcov = curve_fit(arctan_fit,x,y,(y0,x0,A,B,C),maxfev=100000)
    #popt=[y0,x0,A,B]
    return popt

def fit_expcosine(t,y,Tdet,T2):
    y0=y[size(t)-1]
    A=max(y)-min(y)
    t0=0
    popt , pcov = curve_fit(expcosine,t,y,(A,Tdet,y0,t0,T2),maxfev=100000)
    #popt=[A,Tdet,y0,t0,T2]
    return popt

def fit_gausscosine(t,y,Tdet,T2):
    y0=y[size(t)-1]
    A=max(y)-min(y)
    t0=0
    popt , pcov = curve_fit(gausscosine,t,y,(A,Tdet,y0,t0,T2))
    #popt=[A,Tdet,y0,t0,T2]
    return popt

def fit_cosine(x,y,omegaguess):
    index, ymax = max(enumerate(y), key=operator.itemgetter(1))
    x0 = x[index]
    y0 = mean(y)
    A = (ymax-y0)
    #spline = splrep(x,y-y0)
    #roots = sproot(spline)
    #T = (roots[1]-roots[0])*2
    #omega=2*pi/T
    omega=omegaguess
    popt , pcov = curve_fit(cosine,x,y,(A,omega,y0,x0))
    #popt=(A,omega,y0,x0)
    return popt


def load2dsweep(filename,fstart,fstop,pstart,pstop):
    spec=load_file(filename)
    spec=IQrotate(spec)
    freq = linspace(fstart,fstop,size(spec[0]))
    powers=linspace(pstart,pstop,size(spec[:,0]))
    return (spec,freq,powers)

def load2dsweepnoIQrotate(filename,xstart,xstop,ystart,ystop):
    data=load_file(filename)
    x = linspace(xstart,xstop,size(data[0,:]))
    y=linspace(ystart,ystop,size(data[:,0]))
    return (data,x,y)

def analyze2dsweep(spec,freq,powers,nfits,fitAmp=False):
    bw=[]
    frequency=[]
    ymax=[]
    rangeMax=len(powers);rangeMin=rangeMax-nfits;
    for i in range(rangeMin,rangeMax):
        spectofit=abs(spec[i])**2
        if fitAmp:
            spectofit=real(spec[i])

        popt=fit_lorentzian(freq,spectofit,show_plot=False)
        figure(3)
        plot(freq,spectofit/sum(spectofit)+i*(max(spectofit)-min(spectofit)))
        yfit = lorentzian(freq,popt[0],popt[1],popt[2],popt[3])
        plot(freq,yfit/(sum(yfit))+i*(max(spectofit)-min(spectofit)))
        xlabel('Frequency (Hz)')
        ylabel('S21')
        bw.append(2*sqrt(popt[3])/1e6)
        frequency.append(popt[0]/1e9)
        ymax.append(max(spectofit))

#    figure(2)
#    contourf(freq,powers,abs(spec)**2,200)
#    xlabel('Frequency (GHz)')
#    ylabel('Powers (dBm)')
#    title('Spectroscopy - power sweep')

    figure(4)
    plot(powers[len(powers)-nfits:len(powers)],bw,'o')
    xlabel('Power (dBm)')
    ylabel('Bandwidth (MHz)')

    figure(5)
    plot(powers[len(powers)-nfits:len(powers)],frequency,'o')
    xlabel('Power (dBm)')
    ylabel('Frequency (GHz)')

    figure(6)
    plot(powers[len(powers)-nfits:len(powers)],ymax,'o')
    xlabel('Power (dBm)')
    ylabel('Max amplitude (V)')
    return (frequency,bw)


def IQrotate(spec):
    numcol=size(spec[0]);numligns=size(spec)/numcol;
    specrot=spec
    if numcol>1:
        for i in range(numligns):
            speci=spec[i]
            x=real(speci);y=imag(speci);popt=numpy.polyfit(x, y, 1);
            yfit=popt[0]*x+popt[1]
            phase=arctan(popt[0])
            specrot[i]=speci*exp(-1j*phase)
            print phase
    else:
        x=real(spec);y=imag(spec);popt=numpy.polyfit(x, y, 1);
        yfit=popt[0]*x+popt[1]
        phase=arctan(popt[0])
        specrot=spec*exp(-1j*phase)
        print phase*180/pi
    return (specrot)

def IQrotate2D(data):
    dataVec=data.flatten()
    x=real(dataVec);y=imag(dataVec);popt=numpy.polyfit(x, y, 1);
    yfit=popt[0]*x+popt[1]
    phase=arctan(popt[0])
    datarot=data*exp(-1j*phase)
    print phase*180/pi
    return datarot


def loadT2R(filename,dt,numsteps):
    y=loadtxt(filename)
    y=y[0:numsteps]+1j*y[numsteps:2*numsteps]
    y=IQrotate(y)
    x=linspace(0,numsteps*dt,numsteps)
    return(x,y)

def exponential(t,t0,y0,A,tau):
    y=y0+A*np.exp(-(t-t0)/tau)
    return (y)

def fit_exponential(x,y,guess,T2E=False,plotfig=False, t0=0):
    A=max(y)-min(y)
    tau=guess
    y0=y[-1]
    popt , pcov = curve_fit(exponential, x, y,(t0,y0,A,tau), maxfev = 10000)
    yfit=exponential(x,popt[0],popt[1],popt[2],popt[3])
    tau=popt[3]
    data_string = r"T1 ($\mu$s) : %.5f" % (round(tau*1e6,5))
    if T2E:
        data_string = r"T2E ($\mu$s) : %.5f" % (round(tau*1e6,5))

    if plotfig:
        f, ax = subplots()
        ax.plot(x/1e-6,y,'.')
        ax.plot(x/1e-6,yfit,'--')
        text(0.01, 0.5, data_string, bbox=dict(facecolor='red', alpha=0.5),transform=ax.transAxes)
        xlabel('Time (mus)')
        ylabel('Signal')
        print 'T1 = '
        print popt[3]
    return popt

def gaussian(t, A, sigma,t0):
    y=A*np.exp(-((t-t0)/sigma)**2/2)
    return y

def fit_gaussian(x, y, sigma):
    x0=0
    A=max(y)-min(y)
    popt, pcov = curve_fit(gaussian, x, y, (A, sigma,x0))
    return popt

def gauss_exp(t, A, sigma, tau, t0, y0):
    y = exponential(t, t0, 0, 1, tau) * gaussian(t, 1, sigma, t0)
    y *= A
    y += y0
    return y

def fit_gauss_exp(x, y, sigma, tau, t0=0, plotfig=True):
    A = (max(y) - min(y)) / 2
    y0 = y[-1]
    popt, pcov = curve_fit(gauss_exp, x, y, (A, sigma, tau, t0, y0))

    if plotfig:
        yfit = gauss_exp(x, *popt)
        f, ax = subplots()
        ax.plot(x/1e-6,y,'.')
        ax.plot(x/1e-6,yfit,'--')
        data_string = r"$\tau$ ($\mu$s): %.5f" % round(popt[2]*1e6,5) + "\n" + r"$\sigma$ ($\mu$s): %.5f" % round(popt[1]*1e6, 5)
        text(0.01, 0.5, data_string, bbox=dict(facecolor='red', alpha=0.5),transform=ax.transAxes)
        xlabel('Time (mus)')
        ylabel('Signal')
        print 'Tau = ', popt[2], '\nSigma = ', popt[1]

    return popt

def paritydecay(t,kappa,A,y0,nbar0):
    P=y0+A*exp(-2*nbar0*exp(-kappa*t))
    return P

def paritydecay2(t,nbar0,kappa,A):
    P=A*exp(-nbar0*exp(-kappa*t))*cos(nbar0*exp(-kappa*t))
    return P

def Qdecay(t,kappa,A,y0):
    P=y0+A*exp(-nbar0*exp(-kappa*t))
    return P

def fit_paritydecay(t,P,nbar0,kappa):
    y0=P[0]
    A=max(P)-min(P)
    popt , pcov = curve_fit(lambda t,kappa,A,y0:paritydecay(t,kappa,A,y0,nbar0), t, P ,(kappa,A,y0))
    return popt

def fit_paritydecay2(t,P,nbar0,kappa):
    A=max(P)-min(P)
    popt , pcov = curve_fit(paritydecay2, t, P ,(nbar0,kappa,A))
    return popt

def fit_Qdecay(t,P,nbar0,kappa):
    y0=P[0]
    nbar0=nbar0
    A=max(P)-min(P)
    popt , pcov = curve_fit(Qdecay, t, P ,(kappa,A,y0))
    return popt

def getphase(spec):
    numcol=size(spec[0]);numligns=size(spec)/numcol;
    specrot=spec
    if numcol>1:
        for i in range(numligns):
            speci=spec[i]
            x=real(speci);y=imag(speci)
            for j in range(numcol):
                specrot[i][j]=arctan(y[j]/x[j])
    else:
        x=real(spec);y=imag(spec)
        for i in range(numligns):
            specrot[i]=arctan(y[i]/x[i])
    return (specrot)

def getWorstT1(filename,dt,numsteps,fmin,fmax):
    T=numsteps*dt
    (data,time,freq)=load2dsweep(filename,0,T,fstart,fstop)
    data=real(data)
    T1=[]
    for i in range(len(freq)):
        popt=fit_exponential(time,data[i,:],1000*1e-9)
        T1.append(popt[3]*1e6)
    worstT1index, min_value = min(enumerate([abs(t) for t in T1]), key=operator.itemgetter(1))
    return (data[worstT1index,:],freq[worstT1index],T1,data,time,freq)

def nbarfitformula(time,A,B,kappa1ph,kappa2ph,nbar0):
    y=A+B*(exp(-kappa1ph*time)/((1./nbar0)+2*(kappa2ph/kappa1ph)*(1-exp(-kappa1ph*time))))
    return y

def fit_nbar(time,data,kappa1ph,kappa2ph,nbar0):
    A=data[-1]
    B=(data[0]-A)/nbar0
    popt , pcov = curve_fit(nbarfitformula, time, data ,(A,B,kappa1ph,kappa2ph,nbar0))
    return popt

def nbarCircle(time,A,B,C,kappa):
    y=[A+B/(1+C*exp(-kappa*t)) for t in time]
    return y

def fit_nbarCircle(time,data,kappa,C):
    A=data[0]
    B=data[-1]-data[0]
    #popt , pcov = curve_fit(nbarCircle, time, data ,(A,B,C,kappa))
    popt=(A,B,C,kappa)
    return popt


def calc_R(xc, yc):
    """ calculate the distance of each data points from the center (xc, yc) """
    return sqrt((x-xc)**2 + (y-yc)**2)

def f_2b(c):
    """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

def Df_2b(c):
    """ Jacobian of f_2b
    The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
    xc, yc     = c
    df2b_dc    = empty((size(c), x.size))

    Ri = calc_R(xc, yc)
    df2b_dc[0] = (xc - x)/Ri                   # dR/dxc
    df2b_dc[1] = (yc - y)/Ri                   # dR/dyc
    df2b_dc    = df2b_dc - df2b_dc.mean(axis=1)[:, newaxis]

    return df2b_dc

def getnbarfromdata(t,data):
    x=real(data)
    y=imag(data)

    method_2b  = "leastsq with jacobian"


    center_estimate = x_m, y_m
    center_2b, ier = optimize.leastsq(f_2b, center_estimate, Dfun=Df_2b, col_deriv=True)

    xc_2b, yc_2b = center_2b
    Ri_2b        = calc_R(*center_2b)
    R_2b         = Ri_2b.mean()
    residu_2b    = sum((Ri_2b - R_2b)**2)

    circle_fit=plt.Circle((xc_2b,yc_2b),R_2b,color='r')

    fig, ax = subplots()
    ax.plot(real(data),imag(data),'o',0,0)
    fig.gca().add_artist(circle_fit)
    ax.set_aspect('equal')
    ax.set_xlim([xc_2b-R_2b, xc_2b+R_2b])
    ax.set_ylim([yc_2b-R_2b, yc_2b+R_2b])

    fig2, ax2 = subplots()
    ax2.plot(x,y, 'o', 0,0)


    circle_fit=plt.Circle((R_2b,0),R_2b,color='r')

    fig, ax = subplots()
    ax.plot(real(data)-xc_2b+R_2b,imag(data)-yc_2b,'o',0,0)
    fig.gca().add_artist(circle_fit)
    ax.set_aspect('equal')
    ax.set_xlim([-R_2b/10, 2*R_2b])
    ax.set_ylim([-R_2b,R_2b])

    data=data-xc_2b+R_2b-1j*(yc_2b)
    data=1/data
    data=IQrotate(data)
    data=real(data)
    nbar=data

    return nbar

def gaussian_sum(t, A1, sigma1, t01, A2, sigma2, t02):
    y = np.zeros(np.size(t))
    for ii in range(size(y)):
        y[ii] = A1*exp(-((t[ii]-t01)/sigma1)**2/2)+A2*exp(-((t[ii]-t02)/sigma2)**2/2)
    return y

def fit_gaussian_sum(x, y):
    popt, pcov = curve_fit(gaussian_sum, x, y, (A1, sigma1,t01, A2, sigma2,t02),maxfev=5000)
    return popt

#def IQrotate_hist_old(data):
#    re=real(data);
#    im=imag(data);
#    rot_data = mdp.pca(np.array(zip(re,im))-np.array([np.mean(re),np.mean(im)]))
#    re_rot = rot_data[:,0]
#    im_rot = rot_data[:,1]
#    return re_rot+ 1j*im_rot

def fit_hist(data_postselect_raw, numbins=100, plot_hist=True, logplot=True):
    numbins=numbins
    xmin=min(min(np.real(data_postselect_raw)),min(np.imag(data_postselect_raw)))
    xmax =max(max(np.real(data_postselect_raw)),max(np.imag(data_postselect_raw)))

    histi_re, bin_edges_re = np.histogram(np.real(data_postselect_raw), bins=numbins, density=True)
    histi_im, bin_edges_im = np.histogram(np.imag(data_postselect_raw), bins=numbins, density=True)
    hist2D, xedges, yedges = np.histogram2d(np.real(data_postselect_raw), np.imag(data_postselect_raw), [numbins-1, numbins-1])#, range=[[xmin, xmax], [xmin, xmax]])

    x=xedges
    y=histi_re

    gaussroots=[]
    Aroots=y.max()
    #dAroots=np.multiply(y.max(),0.2)
    while np.size(gaussroots)<4:
        Aroots=Aroots/2
        yroots=y-Aroots
        spline = scipy.interpolate.splrep(x, yroots)
        gaussroots = scipy.interpolate.sproot(spline)

    gaussroots=np.sort(gaussroots)
    print gaussroots
    t01=(gaussroots[1]+gaussroots[0])/2
    t02=(gaussroots[3]+gaussroots[2])/2
    listforindex01=abs(xedges-t01)
    index01=listforindex01.argmin()
    listforindex02=abs(xedges-t02)
    index02=listforindex02.argmin()
    A1=max(y[0:index01])
    A2=max(y[index02:-1])
    sigma1=max((gaussroots[3]-gaussroots[2])/2,(gaussroots[1]-gaussroots[0])/2)
    sigma2=sigma1
    popt, pcov = scipy.optimize.curve_fit(gaussian_sum, x, y, (A1, sigma1, t01, A2, sigma2, t02), maxfev=1000000)
    thresholdVec = [np.abs(+math.erf((t-popt[2])/(np.sqrt(2)*np.abs(popt[1])))+math.erf((t-popt[5])/(np.sqrt(2)*np.abs(popt[4])))) for t in xedges]
    thresholdIdx = np.argmin(thresholdVec)
    threshold= xedges[thresholdIdx]
    #threshold=(popt[2]+popt[5])/2

    if plot_hist:
        fs=14
        fsTicks=14
        gaussfit=gaussian_sum(xedges,*popt)
        Pgth=0
        for i in range(np.size(xedges)-1):
            if xedges[i]<threshold:
                Pgth=Pgth+histi_re[i]*(xedges[1]-xedges[0])

        Scurve=np.zeros((numbins))
        for i in range(numbins-1):
            Scurve[i+1]=Scurve[i]+histi_re[i]
        Scurve=Scurve/Scurve[-1]

        Peg=0.5*(1-math.erf((threshold-popt[2])/np.sqrt(2)/np.abs(popt[1])))
        Pge=0.5*(1+math.erf((threshold-popt[5])/np.sqrt(2)/np.abs(popt[4])))

        data_string= r'${\rm P_{g|e}}=$' + str(round(Pge*100,3))+r'%'
        data_string_0= r'${\rm P_{e|g}}=$' + str(round(Peg*100,3))+r'%'
        data_string2= r'P$_{gth}$='+np.str(np.round(Pgth,3)*100)+'%'
        data_string3= r'1-P$_{gth}$='+np.str(100-np.round(Pgth,3)*100)+'%'
        data_string_01= r'${\rm Threshold}=$' + str(round(threshold,3))
        fig, ax = plt.subplots(2, 2, figsize=(10,10))
        if logplot:
            ax[0][0].pcolor(xedges, yedges, np.transpose(np.log(1+hist2D)), cmap='afmhot')
        else:
            ax[0][0].pcolor(xedges, yedges, np.transpose(hist2D), cmap='afmhot')
        ax[0][0].axis('equal')
        ax[0][0].axis('tight')
        ax[1][0].plot(xedges, histi_re,'.',xedges,gaussfit,'--',linewidth=2.0)
        ax[1][0].text(0.5, 0.4, data_string, bbox=dict(facecolor='red', alpha=0.5),transform=ax[1][0].transAxes)
        ax[1][0].text(0.2, 0.4, data_string_0, bbox=dict(facecolor='red', alpha=0.5),transform=ax[1][0].transAxes)
        ax[1][0].text(0.5, 0.2, data_string_01, bbox=dict(facecolor='red', alpha=0.5),transform=ax[1][0].transAxes)
        ax[0][0].text(0.5, 0.7, data_string2, bbox=dict(facecolor='red', alpha=0.5),transform=ax[0][0].transAxes)
        ax[0][0].text(0.5, 0.3, data_string3, bbox=dict(facecolor='red', alpha=0.5),transform=ax[0][0].transAxes)
        ax[1][0].axis('tight')
        ax[1][0].set_ylim([histi_re.max(), 0])
        ax[0][1].plot(histi_im, yedges)
        ax[0][1].axis('tight')
        #ax[1][1].plot(xedges, Scurve)
        #ax[1][1].set_ylim([0,1])
        #ax[1][1].axis('tight')
        plt.show()
    return popt, threshold

def IQrotate_hist(data):
    I=np.real(data)
    Q=np.imag(data)
    Cov=np.cov(I,Q)
    A=scipy.linalg.eig(Cov)
    eigvecs=A[1]
    if A[0][1]>A[0][0]:
        eigvec1=eigvecs[:,0]
    else:
        eigvec1=eigvecs[:,1]
    theta=np.arctan(eigvec1[0]/eigvec1[1])
    return theta

def reS11(x,x0,y0,A,B):
    numerator =  A - ( x - x0 )**2
    denominator =  B + ( x - x0 )**2
    y = y0 + (numerator/denominator)
    return y

def fit_reS11(x,y, B0=1e6**2):
    index, ymax = max(enumerate(y), key=operator.itemgetter(1))
    x0 = x[index]
    y0 = min(y)+1

    #spline = splrep(x, y-(y0+ymax)/2)
    #roots = sproot(spline)
    #whm = roots[-1]-roots[0]

   # y0=-0.00072
    #x0=7.68377*1e9
    #whm=0.1*1e6
#    B = whm*whm/2.0
    A0 = B0*(ymax-y0)
    #print x0,y0,B,A
    #x0=7.68374*1e9;y0=0.0019;B=((100*1e3)**2)/4;A=B*(max(y)-min(y))


    popt, pcov  = curve_fit(reS11,x,y,(x0,y0,A0,B0))
    #popt=[x0,y0,A,B]

    return popt

def reS21_dump(x, x0, y0, A, B, u):
    return reS11(x,x0,y0,A,B) + u*(x - x0)

def fit_reS21_dump(x, y, B0=1e3**2, u0=4e-6):
    index, ymax = max(enumerate(y), key=operator.itemgetter(1))
    x0 = x[index]
    y0 = min(y)+1
    A0 = B0*(ymax-y0)
    popt, pcov  = curve_fit(reS21_dump,x,y,(x0,y0,A0,B0,u0))
    #popt=[x0,y0,A,B,u]

    return popt

def analyze_reflection(freq, ropows, data, f_0=7.226e9, kc=380e3, ki=380e3, a_in=None, T=-400e-9):
    if np.size(ropows)==1:
        data=np.array([data])

    f0Vec=[]
    kappaCover2piVec=[]
    kappaIover2piVec=[]
    ainVec=[]
    TVec = []

    fig1 = plt.figure()
    ax12 = fig1.add_subplot(2,2,2)
    ax11 = fig1.add_subplot(1,2,1)
    ax13 = fig1.add_subplot(2,2,4)
    fig2 = plt.figure()
    ax21 = fig2.add_subplot(1,2,1)
    ax22 = fig2.add_subplot(2,2,2)
    ax23 = fig2.add_subplot(2,2,4)

    for ii in range(np.size(ropows)):
        popt = fit_complex_a_out(freq, data[ii], f_0=f_0, kc=kc, ki=ki, a_in=a_in, T=T)
        freqfit = np.linspace(freq[0], freq[-1], 10 * len(freq))
        fitfunc = complex_a_out(freqfit, *popt)
        ax11.plot(np.real(fitfunc), np.imag(fitfunc))
        ax11.plot(np.real(data[ii]), np.imag(data[ii]),'o')
        ax11.axis('equal')
        ax12.plot(freqfit, np.angle(fitfunc, deg=True))
        ax12.plot(freq, np.angle(data[ii], deg=True), 'o')
        ax13.plot(freqfit, np.abs(fitfunc))
        ax13.plot(freq, np.abs(data[ii]), 'o')
        ax11.set_xlabel('Real')
        ax11.set_ylabel('Imag')
        ax12.set_ylabel('Phase (deg)')
        ax13.set_xlabel('Frequency (GHz)')
        ax13.set_ylabel('Abs')
        ax12.ticklabel_format(useOffset=False)
        ax13.ticklabel_format(useOffset=False)

        fitfunc_mod = fitfunc / popt[3] * np.exp(-1j*popt[4]*(np.asarray(freqfit) - popt[0]))
        data_mod = data[ii] / popt[3] * np.exp(-1j*popt[4]*(np.asarray(freq) - popt[0]))
        ax21.plot(np.real(fitfunc_mod), np.imag(fitfunc_mod))
        ax21.plot(np.real(data_mod), np.imag(data_mod),'o')
        ax21.axis('equal')
        ax22.plot(freqfit, np.angle(fitfunc_mod, deg=True))
        ax22.plot(freq, np.angle(data_mod, deg=True), 'o')
        ax23.plot(freqfit, np.abs(fitfunc_mod))
        ax23.plot(freq, np.abs(data_mod), 'o')
        ax21.set_xlabel('Real')
        ax21.set_ylabel('Imag')
        ax22.set_ylabel('Phase (deg)')
        ax23.set_xlabel('Frequency (GHz)')
        ax23.set_ylabel('Abs')
        ax22.ticklabel_format(useOffset=False)
        ax23.ticklabel_format(useOffset=False)

        print 'Frequency_' + str(ii) + ' = ' + str(popt[0]*1e-9) + ' GHz'
        print '\kappa/2\pi_' + str(ii) + ' = ' + str((popt[1]+popt[2])*1e-6) + ' MHz'
        print 'Electrical delay_' + str(ii) + ' = ' + str(popt[4]*1e9) + ' ns'
        print '\kappa_c/2\pi' + str(ii) + ' = ' + str((popt[1])*1e-6) + ' MHz'
        print '\kappa_i/2\pi' + str(ii) + ' = ' + str((popt[2])*1e-6) + ' MHz'
        print 'T1internal' + str(ii) + ' = ' + str(1./2/np.pi/popt[2]*1e6) + 'us'
        f0Vec.append(popt[0])
        kappaCover2piVec.append(popt[1])
        kappaIover2piVec.append(popt[2])
        ainVec.append(popt[3])
        TVec.append(popt[4])

    if np.size(ropows)>1:
        fig, ax = plt.subplots(5)
        ax[0].plot(ropows, np.asarray(f0Vec)*1e-9)
        ax[0].set_ylabel('f0 (GHz)')
        ax[1].plot(ropows, np.asarray(kappaCover2piVec)*1e-6)
        ax[1].set_ylabel('kappaC/2pi (MHz)')
        ax[2].plot(ropows, np.asarray(kappaIover2piVec)*1e-6)
        ax[2].set_ylabel('kappaI/2pi (MHz)')
        ax[3].plot(ropows, 10*np.log(np.abs(ainVec)**2))
        ax[3].set_ylabel('20*log(|ain|)')
        ax[4].plot(ropows, np.asarray(TVec)*1e9)
        ax[4].set_ylabel('Electrical delay (ns)')
        ax[4].set_xlabel('Power (dBm)')

    return f0Vec, kappaCover2piVec, kappaIover2piVec, ainVec, TVec

def complex_a_out(f, f_0, kc, ki, a_in, T): #kc and ki are kappas/2pi
    D = f - f_0
    #D=f_0-f
    num = - 1j*D + (kc - ki)/2
    den = 1j*D + (kc+ki)/2
    if kc>0 and ki>0 and f_0>0:
        return num/den*a_in*np.exp(1j*D*T)
    else:
        return np.Inf

def fit_complex_a_out(f, a_out, f_0=7.226e9, kc=500e3, ki=500e3, a_in=None, T=0e-9):
    def aux(f, f_0, kc, ki, re_a_in, im_a_in, T):
        return complex_a_out(f, f_0, kc, ki, re_a_in + 1j*im_a_in, T)
    if a_in is None:
        a_in = a_out[0]
    popt = complex_fit(aux, f, a_out, (f_0, kc, ki, np.real(a_in), np.imag(a_in), T))
    return [popt[0], popt[1], popt[2], popt[3] + 1j*popt[4], popt[5]]

def complex_fit(f, xData, yData, p0, weights=None, bounds=()):
    if np.isscalar(p0):
        p0 = np.array([p0])

    def residuals(params, x, y):
        if weights is not None:
            diff = weights * f(x, *params) - y
        else:
            diff = f(x, *params) - y
        flatDiff = np.zeros(diff.size * 2, dtype=np.float64)
        flatDiff[0:flatDiff.size:2] = diff.real
        flatDiff[1:flatDiff.size:2] = diff.imag
        return flatDiff

    popt, bar = leastsq(residuals, p0, args=(xData, yData))
    return popt

def langevinEnvelope(t, t0, t1, a_in, k_c, D, k):
    k_c = np.abs(k_c)
    k = k_c + np.abs(k - k_c)
    z = k/2 + 1j * D
    if  np.isscalar(t):
        t = np.array(t)
    
    i0=np.argmin(np.abs(t-t0))
    i1=np.argmin(np.abs(t-t1))
    a_out = np.zeros(len(t), dtype=np.complex)
    a_out[:i0]=0
    a_out[i0:i1]=a_in * (1 + k_c / z * (np.exp(-z * (t[i0:i1] - t0)) - 1))
    a_out[i1:]=a_in * k_c / z * (np.exp(-z * (t[i1:] - t0)) - np.exp(-z * (t[i1:] - t1)))

    if len(a_out) == 1:
        a_out = a_out[0]

    return a_out

def fitLangevinEnvelope(x, y, t0, t1, k_c=200e3 * 2*np.pi, D=1e6 * 2*np.pi, k=500e3*2*np.pi, cutEdges=False, plotfig=True):
    a_in = y[np.argmin(np.abs(x-t0)) + 2]
    def auxFunc(x, Rea_in, Ima_in, k_c, D, k):
        return langevinEnvelope(x, t0, t1, Rea_in + 1j*Ima_in, k_c, D, k)
    if cutEdges:
        weights = [0 if t >= t1 or t <= t0 else 1 for t in x]
        print weights
        popt = complex_fit(auxFunc, x, y, (a_in.real, a_in.imag, k_c, D, k), weights=weights)
    else:
        popt = complex_fit(auxFunc, x, y, (a_in.real, a_in.imag, k_c, D, k))

    if plotfig:
        Rea_in, Ima_in, k_c, D, k = popt
        a_in = Rea_in + 1j*Ima_in
        a_out = langevinEnvelope(x, t0, t1, a_in, k_c, D, k)

        fig = plt.figure(figsize=(23.5, 8), tight_layout=True)
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 4)

        ax1.plot(np.real(y),np.imag(y))
        ax1.set_xlabel('real')
        ax1.set_ylabel('imag')
        ax1.set_title('aout')
        ax2.plot(x*1e9,np.real(y))
        ax2.set_xlabel('time (ns)')
        ax2.set_ylabel('real')
        ax2.set_title('aout')
        ax3.plot(x*1e9,np.imag(y))
        ax3.set_xlabel('time (ns)')
        ax3.set_ylabel('imag')
        ax3.set_title('aout')

        ax1.plot(np.real(a_out),np.imag(a_out))
        ax1.set_xlabel('real')
        ax1.set_ylabel('imag')
        ax1.set_title('aout')
        ax2.plot(x*1e9,np.real(a_out))
        ax2.set_xlabel('time (ns)')
        ax2.set_ylabel('real')
        ax2.set_title('aout')
        ax3.plot(x*1e9,np.imag(a_out))
        ax3.set_xlabel('time (ns)')
        ax3.set_ylabel('imag')
        ax3.set_title('aout')

    print "Fit results:", "\n", "Detuning: ", popt[3] / 2 / np.pi / 1e3, " kHz", "\n", "ka_c / 2pi: ", (popt[2]) / 2 / np.pi / 1e3, " kHz", "\n", "ka_i / 2pi: ", (popt[4] - popt[2]) / 2 / np.pi / 1e3, " kHz"
    return popt

def jitterLangevin(t, t0, t1, a_in, k_c, D, k, s):
    '''
    Warning: the sign convention here isn't consistent with the sign convention in the Langevin function defined above.
    '''

    z = k/2 - 1j * D
    if  np.isscalar(t):
        t = np.array(t)

    a_out = np.zeros(len(t), dtype=np.complex)
    for i in range(len(t)):
        if t[i] < t0:
            a_out[i] = 0
        elif t[i] < t1:
            a = s / sqrt(2)
            b = z / s / sqrt(2)
            a_out[i] = -a_in * (1 + k_c * np.exp(- b**2) * (scipy.special.erf(-b) - scipy.special.erf(-a * t[i] - b)))
        else:
            a = s / sqrt(2)
            b = z / s / sqrt(2)
            a_t1 = -sqrt(k_c) * np.exp(- b**2) * (scipy.special.erf(-b) - scipy.special.erf(-a * t1 - b))
            a_out[i] = sqrt(k_c) * np.exp(- z * t[i]) * np.exp(- t[i]**2 * s**2 / 2) * a_t1

    if len(a_out) == 1:
        a_out = a_out[0]

    return a_out

def fitJitterLangevin(x, y, t0, t1, k_c=200e3 * 2*np.pi, D=5e5 * 2*np.pi, k=500e3*2*np.pi, s=10e3 * 2 * np.pi, cutEdges=False, plotfig=True):
    a_in = y[np.argmin(np.abs(x-t0)) + 2]
    def auxFunc(x, Rea_in, Ima_in, k_c, D, k, s):
        return jitterLangevin(x, t0, t1, Rea_in + 1j*Ima_in, k_c, D, k, s)
    if cutEdges:
        weights = [0 if t >= t1 or t <= t0 else 1 for t in x]
        print weights
        popt = complex_fit(auxFunc, x, y, (a_in.real, a_in.imag, k_c, D, k, s), weights=weights)
    else:
        popt = complex_fit(auxFunc, x, y, (a_in.real, a_in.imag, k_c, D, k, s))

    if plotfig:
        Rea_in, Ima_in, k_c, D, k, s = popt
        a_in = Rea_in + 1j*Ima_in
        a_out = jitterLangevin(x, t0, t1, a_in, k_c, D, k, s)

        fig = plt.figure(figsize=(23.5, 8), tight_layout=True)
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 4)

        ax1.plot(np.real(y),np.imag(y))
        ax1.set_xlabel('real')
        ax1.set_ylabel('imag')
        ax1.set_title('aout')
        ax2.plot(x*1e9,np.real(y))
        ax2.set_xlabel('time (ns)')
        ax2.set_ylabel('real')
        ax2.set_title('aout')
        ax3.plot(x*1e9,np.imag(y))
        ax3.set_xlabel('time (ns)')
        ax3.set_ylabel('imag')
        ax3.set_title('aout')

        ax1.plot(np.real(a_out),np.imag(a_out))
        ax1.set_xlabel('real')
        ax1.set_ylabel('imag')
        ax1.set_title('aout')
        ax2.plot(x*1e9,np.real(a_out))
        ax2.set_xlabel('time (ns)')
        ax2.set_ylabel('real')
        ax2.set_title('aout')
        ax3.plot(x*1e9,np.imag(a_out))
        ax3.set_xlabel('time (ns)')
        ax3.set_ylabel('imag')
        ax3.set_title('aout')

    print "Fit results:", "\n", "Detuning: ", popt[3] / 2 / np.pi / 1e3, " kHz", "\n", "ka_c / 2pi: ", (popt[2]) / 2 / np.pi / 1e3, " kHz", "\n", "ka_i / 2pi: ", (popt[4] - popt[2]) / 2 / np.pi / 1e3, " kHz", "\n", "Sigma_Delta / 2pi: ", popt[5] / 2 / np.pi / 1e3
    return popt

def jitterLangevin_debug(t, t0, t1, a_in, k_c, D, k, s):
    '''
    Warning: the sign convention here isn't consistent with the sign convention in the Langevin function defined above.
    '''

    dtau=1e-9
    z = k/2 - 1j * D
    if  np.isscalar(t):
        t = np.array(t)
    i0=0
    a_out = np.zeros(len(t), dtype=np.complex)
    for i in range(len(t)):
        if t[i] < t0:
            a_out[i] = 0
            i0+=1
        elif t[i] < t1:
            for tau in np.linspace(0,t[i],int(t[i]/dtau)):
                a_out[i] += -dtau*np.exp(z*(tau-t[i]))*np.exp(-(tau-t[i])**2*s**2/2)*a_in*k_c
            a_out[i]-=a_in
            i0+=1
        else:
            a_out[i] = np.exp(- z * t[i]) * np.exp(- t[i]**2 * s**2 / 2) * (a_out[i0]+a_in)

    if len(a_out) == 1:
        a_out = a_out[0]

    return a_out
    
def autocorrelation(x):
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    result = r/(variance*(np.arange(n, 0, -1)))
    return result