import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
import math
import cmath
from fmSupportLib import fmDemodArctan, fmPlotPSD

# From constraints document, the settings for group 18
# 						Mode 0 	Mode 1	Mode 2	Mode 3
# RF Fs (Ksamples/s)	2400	1152	2400	1920	#radio-frequency (RF) sampling rate
# IF Fs (Ksamples/s)	240		288		240		320
# Audio Fs (Ksamples/s)	48		48		44.1	44.1

rf_Fc = 100e3	# the cutoff frequency to extract the FM channel from raw IQ data
if_Fc = 16e3
rf_taps = 151 	# the number of taps for the low-pass filter to extract the FM channel

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
rf_Fs2 = 2400e3
if_Fs2 = 240e3
rf_U2 = 1
rf_D2 = 10
if_U2 = 147
if_D2 = 800
audio_Fs2 = 44.1e3

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

def fastConvWithState(h,x,U,D,state):
	y= np.zeros(math.floor((len(x))*U/D))

	# print(math.floor((len(x))*U/D))
	for n in range(len(y)):
		k = (n*D)%U
		for s in range(math.ceil((len(h)-k)/U)):
			i = k + s*U
			if n*D-i >= 0:
				y[n] += np.real(h[i] * x[int((n*D-i)/U)])
			else:
				# print(((len(state)-1)+(n*D-i))/U)
				y[n] += np.real(h[i]*state[int(((len(state)-1)+(n*D-i))/U)])
	state = x[(len(x)-len(h)+1):]

	return y,state

def fmDemodDeriv(I,Q,prev_I,prev_Q):
	fmDedom = np.empty(len(I))
	print(len(I))
	for k in range(len(I)):
		fmDedom[k]= ((I[k] * (Q[k]-prev_Q) )- (Q[k]*(I[k]-prev_I)))/(I[k]**2+Q[k]**2)
		prev_Q = Q[k]
		prev_I = I[k]
	return fmDedom, I[k], Q[k]


def fmDemodArctan(I, Q, prev_phase = 0.0):
#
# the default prev_phase phase is assumed to be zero, however
# take note in block processing it must be explicitly controlled

	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))

	# iterate through each of the I and Q pairs
	for k in range(len(I)):

		# use the atan2 function (four quadrant version) to detect angle between
		# the imaginary part (quadrature Q) and the real part (in-phase I)
		current_phase = math.atan2(Q[k], I[k])

		# we need to unwrap the angle obtained in radians through arctan2
		# to deal with the case when the change between consecutive angles
		# is greater than Pi radians (unwrap brings it back between -Pi to Pi)
		[prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

		# take the derivative of the phase
		fm_demod[k] = current_phase - prev_phase

		# save the state of the current phase
		# to compute the next derivative
		prev_phase = current_phase

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_phase

if __name__ == "__main__":
	in_fname = "../data/99-9fmFs2-4M.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8') # 8-bits unsigned (and interleaved)
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	iq_data = (np.float32(raw_data) - 128.0)/128.0 # IQ data is normalized between -1 and +1 in 32-bit float format
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# Create LPF Coeffs
	rf_coeff = signal.firwin(rf_taps * rf_U2, rf_Fc/(rf_Fs2/2), window=('hann'))
	audio_coeff = my_own_coeff(if_Fc, if_Fs2,rf_taps * rf_U2)

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# Initializing block processing variables
	# set block size
	block_size = 1024 * int(rf_U2/rf_D2) * int(rf_U2/if_D2) * 2 	#select a block_size that is a multiple of KB and a multiple of decimation factors
	block_count = 0
	# init state vars
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	prev_I = 0
	prev_Q = 0
	state_mono_lpf_16k = np.zeros(rf_taps-1)	# add state as needed for the mono channel filter
	# audio buffer that stores all the audio blocks
	audio_data = np.array([])

	while (block_count+1)*block_size < len(iq_data):	#if the number of samples in the last block is less than the block size, it is ignored
		# print('Processing block ' + str(block_count))
		# RF Processing
		i_ds, state_i_lpf_100k = fastConvWithState(rf_coeff, iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
													rf_U2,rf_D2,state_i_lpf_100k)
		q_ds, state_q_lpf_100k = fastConvWithState(rf_coeff, iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
													rf_U2,rf_D2,state_q_lpf_100k)

		# Demodulation and IF Processing
		fm_demod,prev_I,prev_Q = fmDemodDeriv(i_ds,q_ds,prev_I,prev_Q)
		audio_block,state_mono_lpf_16k = fastConvWithState(audio_coeff,fm_demod,if_U2,if_D2,state_mono_lpf_16k)

		audio_data = np.concatenate((audio_data, audio_block))

		if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
			ax0.clear()
			fmPlotPSD(ax0, fm_demod, (rf_Fs2/rf_D2)/1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after downsampling mono audio
			fmPlotPSD(ax2, audio_block, audio_Fs2/1e3, subfig_height[2], \
					'Downsampled Mono Audio (block ' + str(block_count) + ')')

			# save figure to file
			fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = "../data/fmMonoBlock.wav"
	wavfile.write(out_fname, int(audio_Fs2), np.int16((audio_data/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	plt.show()
