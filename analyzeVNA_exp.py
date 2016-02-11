import analysis_functions as af
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
data_dir = "/Users/leghtas/Documents/Travail/Projects/github/MeasurementVNA/data/"
data_filename = "029.dat"
plotdata = False
debug_fits = False

data = np.loadtxt(data_dir + data_filename, skiprows=1)
VNApts_times_genpoints,n=np.shape(data)

VNApts = 101
genpoints = VNApts_times_genpoints/VNApts
genfreq_start = 1e9
genfreq_stop = 20e9

vna_freqs = data[:,0][0:VNApts]*1e9

genfreqs = np.linspace(genfreq_start, genfreq_stop, genpoints)
data_re = data[:,1]
data_im = data[:,2]

data_re = data_re.reshape(genpoints,VNApts)
data_im = data_im.reshape(genpoints,VNApts)
#

#
data_complex = data_re + 1j*data_im
#

freqs_fit = np.zeros(genpoints)
kc_fit = np.zeros(genpoints)
ki_fit = np.zeros(genpoints)
f_0_guess = 13052250602
kc_guess = 18736539
ki_guess = 5465667
for n in range(genpoints):
    data_n = data_complex[n]
    popt = af.fit_complex_a_out(vna_freqs, data_n, f_0=None, kc=kc_guess, ki=ki_guess)
    freqs_fit[n]=popt[0] 
    kc_fit[n]=popt[1]
    ki_fit[n]=popt[2]
    if debug_fits:
        data_n_fit = af.complex_a_out(vna_freqs, *popt)
        fig, ax = plt.subplots(2,2)
        ax[0,0].plot(vna_freqs, np.abs(data_n_fit))
        ax[0,0].plot(vna_freqs, np.abs(data_n))
        ax[1,0].plot(vna_freqs, np.angle(data_n_fit))
        ax[1,0].plot(vna_freqs, np.angle(data_n))
        ax[1,1].plot(np.real(data_n_fit),np.imag(data_n_fit))
        ax[1,1].plot(np.real(data_n_fit),np.imag(data_n))

fig, ax =plt.subplots(3,1)
ax[0].plot(genfreqs*1e-9, freqs_fit*1e-9)
ax[0].set_ylabel('Freq (GHz)')
ax[0].ticklabel_format(useOffset=False)
ax[1].plot(genfreqs*1e-9, kc_fit*1e-6)
ax[1].set_ylabel('kc (MHz)')
ax[1].ticklabel_format(useOffset=False)
ax[2].plot(genfreqs*1e-9, ki_fit*1e-6)
ax[2].set_ylabel('ki (MHz)')
ax[2].set_xlabel('Pump frequency (GHz)')
ax[2].ticklabel_format(useOffset=False)


if plotdata:
    fig, ax = plt.subplots(2)
    ax[0].imshow(np.angle(data_complex, deg=1),aspect='auto', extent=(vna_freqs[0]*1e-9, vna_freqs[-1]*1e-9, genfreqs[0]*1e-9, genfreqs[-1]*1e-9))
    ax[0].set_title('phase(a)')
    #ax[0].set_xlabel('VNA frequency (GHz)')
    ax[0].set_ylabel('Generator frequency (GHz)')
    ax[0].ticklabel_format(useOffset=False)
    
    ax[1].imshow(np.abs(data_complex),aspect='auto', extent=(vna_freqs[0]*1e-9, vna_freqs[-1]*1e-9, genfreqs[0]*1e-9, genfreqs[-1]*1e-9))
    ax[1].set_title('|a|')
    ax[1].set_xlabel('VNA frequency (GHz)')
    ax[1].set_ylabel('Generator frequency (GHz)')
    ax[1].ticklabel_format(useOffset=False)