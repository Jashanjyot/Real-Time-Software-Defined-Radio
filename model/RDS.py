
import numpy as np
from scipy.io import wavfile
from scipy import signal
import math
import cmath
from bitstring import BitStream

import matplotlib.pyplot as plt

from stereo import bpf, mixer, my_own_coeff, fastConv, fmDemodDeriv
from fmRRC import impulseResponseRootRaisedCosine

rf_Fc = 100e3
if_Fc = 16e3
rf_taps = 11
stereo_taps = 11

#  Mode 0
rf_Fs0 = 2400e3
if_Fs0 = 240e3
rf_U0 = 1
rf_D0 = 10
if_U0 = 1
if_D0 = 5
audio_Fs0 = 48e3

# RDS variables
RDS_extract_Fb = 54e3
RDS_extract_Fe = 60e3
RDS_carrier_Fb = 113.5e3
RDS_carrier_Fe = 114.5e3
RDS_PLL_freq = 114e3
RDS_ncoScale = 0.5
RDS_Fc = 3e3
RDS_taps = 11
RDS_gain = 1
SPS_0 = 15
SPS_2 = 38

def fmPll(pllIn, freq, Fs, ncoScale = 0.5, phaseAdjust = 0.0, normBandwidth = 0.01):

	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp
	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	# output arrays for the NCO
	ncoOut_i = np.empty(len(pllIn)+1)
	ncoOut_q = np.empty(len(pllIn)+1)

	# initialize internal state
	integrator = 0.0
	phaseEst = 0.0
	feedbackI = 1.0
	feedbackQ = 0.0
	ncoOut_i[0] = 1.0
	ncoOut_q[0] = 1.0
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
		ncoOut_i[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
		# the quadrature output is the in-phase tone delayed by pi/2
		ncoOut_q[k+1] = math.cos(trigArg*ncoScale + phaseAdjust + math.pi/2)

	# for block processing you should also return the state
	return ncoOut_i, ncoOut_q

# functionality: delay signal/ pad with zeros
def allpass(x, state):
	delayed_x = state + x[0:len(x)-len(state)]	#delay x
	state = x[len(x)-len(state):len(x)]			#save state

	return delayed_x,state

# int main:
# 	while(block)
# 	state = np.zeros((N_taps-1)/2*U/D)
# 	delayed_x, state = allpass(x, state, N_taps,U,D)


def gcd(a, b): 
	while (b):     
		temp = b
		b = (a % b)
		a = temp
	return a

def rationalResampler(input, Fs_in, SPS_mode):

	# if_Fs0 = 240e3 to Fs_out = 15*2375 = 35.625e3

	# desired output frequency Fs_out
	RDS_symbol_rate = SPS_mode*2375
	factor = gcd(Fs_in,RDS_symbol_rate)
	
	RDS_U = RDS_symbol_rate/factor
	RDS_D = Fs_in/factor

	# apply up and down sampling to input signal
	output_signal = input*RDS_U*RDS_D

	return RDS_U, RDS_D, output_signal, RDS_symbol_rate

# • Construct a Quadrature Path for debugging purpose
# • Plot the constellation diagram before Manchester Decoding
# • Identify constellation issues, and tune PLL / APF until close to 100% power is In-Phase
# • Then, you can move on to Manchester Decoding + Frame Synchronization
# #differential decoding function:
diff_last_bit=0 
def manchester_decode(encoded_bits, diff_last_bit):
    
    decoded_bits = [] # list to store the decoded bits 
    
    for i in range(len(encoded_bits)):
        if i == 0:#xor first bit with this parameter to recover the original bit
            decoded_bits.append(encoded_bits[i] ^ diff_last_bit) # append the decoded bit to the list
        else: #if it's not the first bit,
            decoded_bits.append(encoded_bits[i] ^ decoded_bits[i-1])# XOR current w the previously decoded bit to recover the original bit
	
    diff_last_bit = encoded_bits[-1]   # make it the last bit of the current encoded symbol
    decoded_data = BitStream(decoded_bits) # create a new BitStream object from the decoded_bits list.
   
    return decoded_data, diff_last_bit

def frame_synchronization(decoded_bits):

	return 0

if __name__ == "__main__":
	# read the raw IQ data from the recorded file
	in_fname = "C:/Users/clara/OneDrive/Documents/McMaster/Year 3/Winter 2023/COMPENG 3DY4/Labs/3dy4-project-group18-prj-wednesday/data/stereo_l0_r9.raw"
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

	print("Starting RDS")

	# *** RDS Channel Extraction ***
    # bandpass filter
	RDS_extract_coeff = bpf(RDS_extract_Fb, RDS_extract_Fe, if_Fs0, RDS_taps)
	RDS_channel = fastConv(RDS_extract_coeff, fm_demod, 1, 1)

	# *** RDS Carrier Recovery ***
    # squaring nonlinearity - pointwise multiplication of each sample with itself
    # basically the mixer function with both inputs as RDS_channel
	RDS_channel_sq = np.zeros(len(RDS_channel))

	for i in range (len(RDS_channel)):
		RDS_channel_sq[i] = RDS_channel[i]*RDS_channel[i]

    # bandpass filter
	RDS_carrier_coeff  = bpf(RDS_carrier_Fb, RDS_carrier_Fe, if_Fs0, RDS_taps)
	RDS_carrier = fastConv(RDS_carrier_coeff, RDS_channel_sq, 1, 1)

    # PLL and NCO
	RDS_ncoOut_i, RDS_ncoOut_q = fmPll(RDS_carrier, RDS_PLL_freq, if_Fs0, RDS_ncoScale)
	# **use quadrature only to tune the PLL phase**
	# **duplicate mixer, LPF, resampler, RRC for i and q samples

    # all-pass filter
	state = np.zeros(int((RDS_taps-1)/2*if_U0/if_D0))
	RDS_channel_delayed, state = allpass(RDS_channel, state)

    # *** RDS Demodulation ***
    # mixer (multiply pll output with filtered RDS channel)
	RDS_mixer_prod_i = mixer(RDS_ncoOut_i, RDS_channel_delayed)
	RDS_mixer_prod_q = mixer(RDS_ncoOut_q, RDS_channel_delayed)
    
    # low pass filter
	RDS_LPF_carrier_coeff = my_own_coeff(RDS_Fc, if_Fs0, RDS_taps)
	mixed_filtered_result_i = fastConv(RDS_mixer_prod_i, RDS_LPF_carrier_coeff, 1, 1)
	mixed_filtered_result_q = fastConv(RDS_mixer_prod_q, RDS_LPF_carrier_coeff, 1, 1)

    # apply rational resampler to signal after LPF
	RDS_U_i, RDS_D_i, resampled_signal_i, RDS_symbol_rate_i = rationalResampler(mixed_filtered_result_i, if_Fs0, SPS_0)
	RDS_U_q, RDS_D_q, resampled_signal_q, RDS_symbol_rate_q = rationalResampler(mixed_filtered_result_q, if_Fs0, SPS_0)

    # root-raised cosine filter
	RRC_coeff_i = impulseResponseRootRaisedCosine(RDS_symbol_rate_i, RDS_taps)
	RRC_filtered_result_i = fastConv(resampled_signal_i, RRC_coeff_i, 1, 1)

	RRC_coeff_q = impulseResponseRootRaisedCosine(RDS_symbol_rate_q, RDS_taps)
	RRC_filtered_result_q = fastConv(resampled_signal_q, RRC_coeff_q, 1, 1)

	# plot i and q samples pre clock data recovery -> constellation diagram
	your_figure = plt.figure()
	your_figure.scatter(RRC_filtered_result_i, RRC_filtered_result_q, s=10)
	your_figure.set_xlabel('Q')
	your_figure.set_ylabel('I')
	your_figure.set_title('IQ Constellation diagram')
	plt.show()

    # clock and data recovery
	# check the local max and min of the input wave - take derivative?
	# check if H-L -> 1, if L-H -> 0
	# output is a digital clock signal?
	clock_data_recovery =0


	# *** RDS Data Processing ***
    # Manchester and differenital decoding
	decoded_bits, diff_last_bit = manchester_decode(clock_data_recovery,diff_last_bit)

    # frame synchronization and error detection
    
	# RDS application layer




