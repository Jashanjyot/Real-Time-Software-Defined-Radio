#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

"""
The command-line instructions for recording RF data are only for those
who have the RF dongle (Nooelec NESDR Smart v4 Bundle) available at home.
After you have installed the drivers to work with the RF dongle,
the 8-bit unsigned values for the I/Q pairs can be recorded as follows:
rtl_sdr -f 99.9M -s 2.4M - > iq_samples.raw
The above assumes that we are tuned to the FM station at 99.9 MHz,
we use an RF sample rate of 2.4 Msamples/sec and our file is called
iq_samples.raw (change as you see fit).
For the above use case, the data acquisition runs indefinitely,
hence the recording needs to be stopped by pressing Ctrl+C.
If we wish to stop it after a pre-defined number of samples,
e.g., 12 million I/Q pairs (5 seconds at 2.4 Msamples/sec),
we can use an extra argument:
rtl_sdr -f 99.9M -s 2.4M -n 12000000 - > iq_samples.raw
To check if the raw I/Q data has been recorded properly, place the file
in the "data" sub-folder from your project repo and run this Python file
from the "model" sub-folder. It should produce both the .png image files
(of the PSD estimates) for a couple of blocks, as well as the .wav file.
In the source code below (check lines 90+) you can observe where the
raw_data is read and the normalization of the 8-bit unsigned I/Q samples
to 32-bit float samples (in the range -1 to +1) is done; while the
32-bit floats and the range -1 to +1 are optional choices (used by
many third-party SDR software implementations), it is at the discretion
of each project group to decide how to handle the 8-bit unsigned I/Q samples
in their Python model and C++ implementation.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
import math
import cmath

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD

# From constraints document, the settings for group 18
# 						Mode 0 	Mode 1	Mode 2	Mode 3
# RF Fs (Ksamples/s)	2400	1152	2400	1920	#radio-frequency (RF) sampling rate
# IF Fs (Ksamples/s)	240		288		240		320
# Audio Fs (Ksamples/s)	48		48		44.1	44.1

rf_Fc = 100e3	# the cutoff frequency to extract the FM channel from raw IQ data
if_Fc = 16e3
rf_taps = 151 # the number of taps for the low-pass filter to extract the FM channel

# Mode 0
rf_Fs0 = 2400e3
if_Fs0 = 240e3
rf_U0 = 1
rf_D0 = 10
if_U0 = 1
if_D0 = 5
audio_Fs0 = 48e3

# Mode 1
up1 = 1
down1 = 6
gain1 = 1

# Mode 2
up2 = 147
down2 = 800
gain2 = 147

def bpf(Fb,Fe, Fs,N_taps):
	Norm_c= ((Fe+Fb)/(2))/(Fs/2)
	Norm_p = ((Fe-Fb))/(Fs/2)
	h= np.zeros(N_taps, dtype = complex)
	for i in range(N_taps):
		if i == ((N_taps-1)/2):
			h[i] = Norm_p
		else:
			h[i] = Norm_p*(cmath.sin(cmath.pi* Norm_p*(i-(N_taps-1)/2))/(cmath.pi*Norm_p*(i-(N_taps-1)/2)))
		h[i] = h[i]*(cmath.cos(i*cmath.pi*Norm_c))
		h[i] = h[i]*(cmath.sin(i*cmath.pi/N_taps))**2

	return h


def my_own_coeff(Fc, Fs,N_taps):
	Normco= Fc/(Fs/2)
	h= np.zeros(N_taps, dtype = complex)
	for i in range(N_taps):
		if i == ((N_taps-1)/2):
			h[i] = Normco
		else:
			h[i] = Normco*(cmath.sin(cmath.pi* Normco*(i-(N_taps-1)/2))/(cmath.pi*Normco*(i-(N_taps-1)/2)))
		h[i] = h[i]*(cmath.sin(i*cmath.pi/N_taps))**2

	return h

def fastConv(h,x,U,D):
	y= np.zeros(math.floor((len(h)+len(x)-1)*U/D))

	for n in range(len(y)):
		k = (n*D)%U
		for s in range(math.ceil((len(h)-k)/U)):
			i = k + s*U
			if n*D-i >= 0 and n*D-i < len(x)*U:
				y[n] += h[i] * x[int((n*D-i)/U)]

	return y

if __name__ == "__main__":
	# read the raw IQ data from the recorded file
	in_fname = "../data/99-9fmFs2-4M.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8') # 8-bits unsigned (and interleaved)
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	iq_data = (np.float32(raw_data) - 128.0)/128.0 # IQ data is normalized between -1 and +1 in 32-bit float format
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs0/2), window=('hann'))

	# filter to extract the FM channel (I samples are even, Q samples are odd)
	i_ds = fastConv(rf_coeff, iq_data[0::2],rf_U0,rf_D0)
	q_ds = fastConv(rf_coeff, iq_data[1::2],rf_U0,rf_D0)

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# FM demodulator (check the library)
	fm_demod, dummy = fmDemodArctan(i_ds, q_ds)	# we use a dummy because there is no state for this single-pass model

	# PSD after FM demodulation
	fmPlotPSD(ax0, fm_demod, (rf_Fs0/rf_D0)/1e3, subfig_height[0],\
			'Demodulated FM (full recording)')

	audio_coeff = my_own_coeff(if_Fc,if_Fs0,rf_taps)
	audio_data = fastConv(audio_coeff, fm_demod, if_U0,if_D0)

	# PSD after decimating mono audio
	fmPlotPSD(ax2, audio_data, audio_Fs0/1e3, subfig_height[2], 'Downsampled Mono Audio')

	# save PSD plots
	fig.savefig("../data/fmMonoBasic.png")
	plt.show()

	# write audio data to file (assumes audio_data samples are -1 to +1)
	out_fname = "../data/fmMonoBasic.wav"
	wavfile.write(out_fname, int(audio_Fs0), np.int16((audio_data/2)*32767)) #scaling specific for mono audio
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")
