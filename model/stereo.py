
import numpy as np
from scipy.io import wavfile
from scipy import signal
import math
import cmath

rf_Fc = 100e3
if_Fc = 16e3
rf_taps = 51
stereo_taps = 51
#
# // Mode 0
rf_Fs0 = 2400e3
if_Fs0 = 240e3
rf_U0 = 1
rf_D0 = 10
if_U0 = 1
if_D0 = 5
audio_Fs0 = 48e3

def bpf(Fb,Fe, Fs,N_taps):
	Norm_c= ((Fe+Fb)/(2))/(Fs/2)
	Norm_p = ((Fe-Fb))/(Fs/2)
	h= np.zeros(N_taps, dtype = complex)
	for i in range(N_taps):
		if i == ((N_taps-1)/2):
			h[i] = Norm_p
		else:
			h[i] = Norm_p*(cmath.sin(cmath.pi* Norm_p/2*(i-(N_taps-1)/2))/(cmath.pi*Norm_p/2*(i-(N_taps-1)/2)))
		h[i] = h[i]*(cmath.cos(i*cmath.pi*Norm_c))
		h[i] = h[i]*(cmath.sin(i*cmath.pi/N_taps))**2

	return h


def fmPll(pllIn, freq, Fs, ncoScale = 2.0, phaseAdjust = 0.0, normBandwidth = 0.01):

	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp
	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	# output array for the NCO
	ncoOut = np.empty(len(pllIn)+1)

	# initialize internal state
	integrator = 0.0
	phaseEst = 0.0
	feedbackI = 1.0
	feedbackQ = 0.0
	ncoOut[0] = 1.0
	trigOffset = 0
	# note: state saving will be needed for block processing

	for k in range(len(pllIn)):

		# phase detector
		errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential

		# four-quadrant arctangent discriminator for phase error detection
		errorD = math.atan2(errorQ, errorI)

		# loop filter
		integrator = integrator + Ki*errorD

		# update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator

		# internal oscillator
		trigOffset += 1
		trigArg = 2*math.pi*(freq/Fs)*(trigOffset) + phaseEst
		feedbackI = math.cos(trigArg)
		feedbackQ = math.sin(trigArg)
		ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)

	# for stereo only the in-phase NCO component should be returned
	# for block processing you should also return the state
	return ncoOut
	# for RDS add also the quadrature NCO component to the output


def mixer(v1, v2):
	product = np.zeros(len(v2))
	for i in range(len(v2)):
		product[i] = v1[i] * v2[i] * 2
	return product



def combiner(mono_data,stereo_data):
	left_channel= np.zeros(len(stereo_data))
	right_channel= np.zeros(len(stereo_data))

	for k in range(min(len(stereo_data),len(mono_data))):
		left_channel[k] = stereo_data[k] + mono_data[k]
		right_channel[k] = stereo_data[k] - mono_data[k]

	return left_channel,right_channel



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
		for s in range((n*D)%U,len(h),U):
			x_idx = int(n*D-s/U)
			if x_idx >= 0 and x_idx < len(x)*U:
				y[n] += h[s] * np.real(x[x_idx]) * U

	return y

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

def fmDemodDeriv(I,Q,prev_I = 0,prev_Q = 0):
	fmDedom = np.empty(len(I))
	for k in range(len(I)):
		fmDedom[k]= ((I[k] * (Q[k]-prev_Q) )- (Q[k]*(I[k]-prev_I)))/(I[k]**2+Q[k]**2)
		if(math.isnan(fmDedom[k]) or math.isinf(fmDedom[k])): fmDedom[k] = 0
		prev_Q = Q[k]
		prev_I = I[k]
	return fmDedom, I[k], Q[k]



if __name__ == "__main__":
	# read the raw IQ data from the recorded file
	in_fname = "../data/stereo_l0_r9.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8') # 8-bits unsigned (and interleaved)
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	iq_data = (np.float32(raw_data) - 128.0)/128.0 # IQ data is normalized between -1 and +1 in 32-bit float format
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs0/2), window=('hann'))

	# filter to extract the FM channel (I samples are even, Q samples are odd)
	i_ds = fastConv(rf_coeff, iq_data[0::2],rf_U0,rf_D0)
	q_ds = fastConv(rf_coeff, iq_data[1::2],rf_U0,rf_D0)

	# FM demodulator (check the library)
	fm_demod, dummy,dummy = fmDemodDeriv(i_ds, q_ds)	# we use a dummy because there is no state for this single-pass model
	print("RF finished")
	audio_coeff = my_own_coeff(if_Fc,if_Fs0,rf_taps)
	audio_data = fastConv(audio_coeff, fm_demod, if_U0,if_D0)
	print("Mono finished")
	bpf_pilot = bpf(Fb = 18.5e3,Fe =19.5e3, Fs = if_Fs0, N_taps = stereo_taps)
	carrier = fastConv(bpf_pilot, fm_demod,1,1)
	ncoOut = fmPll(carrier, freq= 19e3, Fs=if_Fs0, ncoScale = 2.0, phaseAdjust = 0.0, normBandwidth = 0.01)

	bpf_LR = bpf(Fb =22e3,Fe = 54e3, Fs=if_Fs0, N_taps = stereo_taps)
	channel = fastConv(bpf_LR, fm_demod,1,1)

	mixerProduct  = mixer(ncoOut,channel)
	lpf_dig_filt = my_own_coeff(15000, if_Fs0, N_taps = stereo_taps)
	stereo_block = fastConv(lpf_dig_filt,mixerProduct, if_U0 , if_D0 )
	print("Stereo finished")
	np.pad(audio_data, (int((stereo_taps-1)/2),0), 'constant')
	left, right  = combiner(audio_data,stereo_block)
	
	out_fname = "../data/stereo_l0_r9.wav"
	# for i in range(len(left)):
	# 	left[i] = np.int16((left[i]/2)*32767)
	# 	right[i] = np.int16((right[i]/2)*32767)
	LR_data = np.vstack((left,right))
	LR_data = LR_data.transpose()
	wavfile.write(out_fname, int(audio_Fs0), np.int16((LR_data/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")
