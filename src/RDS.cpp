#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "dy4.h"
#include "fourier.h"
#include "filter.h"
#include "RDS.h"
#include "stereoProcessing.h"

int gcd(int a, int b){
    while(b){
        int temp = b;
        b= (a%b);
        a= temp; 
    }
    return a;
}

std::vector<double> rationalResampler(std::vector<double> input, int Fs_in, int SPS_mode, int length, std::vector<double> RDS_U,std::vector<double> output_signal, int RDS_symbol_rate){
    RDS_symbol_rate = SPS_mode*2375;
    int factor = gcd(Fs_in, RDS_symbol_rate);
    RDS_U = RDS_symbol_rate/factor;
    RDS_D = Fs_in/factor;
    // apply up and down sampling to input signal
    for (int i = 0; i < length; i++) {
        output_signal[i] = input[i] * RDS_U * RDS_D;
    }

    
}

int main()
{
	// *** RDS Channel Extraction ***
    // bandpass filter
    std::vector<float> RDS_extract_coeff;
    std::vector<float> RDS_channel;
    std::vector<float> RDS_channel_state;

    bpf(RDS_extract_Fb, RDS_extract_Fe, params.if_Fs, bpf_taps, RDS_extract_coeff);
    fastConvWithState(RDS_channel, fm_Demod, RDS_extract_coeff, RDS_channel_state, 1, 1, 1);


	// *** RDS Carrier Recovery ***
    // squaring nonlinearity - pointwise multiplication of each sample with itself
    // basically the mixer function with both inputs as RDS_channel?
    for(int i=0; i < RDS_channel.size(); i++){
        RDS_channel[i] = RDS_channel[i]*RDS_channel[i];
    }

    // bandpass filter
    std::vector<float> RDS_carrier_coeff;
    std::vector<float> RDS_carrier;
    std::vector<float> RDS_carrier_state;

    bpf(RDS_carrier_Fb, RDS_carrier_Fe, params.if_Fs, bpf_taps, RDS_carrier_coeff);
    fastConvWithState(RDS_carrier, RDS_channel, RDS_carrier_coeff, RDS_carrier_state, 1, 1, 1);

    // all-pass filter
    int apf_taps = 51;
    std::vector<float> RDS_carrier_delayed(RDS_carrier.size());
    std::vector<float> delayed_state((apf_taps-1)/2);
    std::copy(delayed_state.begin(),delayed_state.end(),RDS_carrier_delayed.begin());	                    //copy in state
    std::copy(RDS_carrier.begin(),RDS_carrier.end()-delayed_state.size(),RDS_carrier_delayed.begin()+delayed_state.size());	//copy relevant part of current block
    std::copy(RDS_carrier.end()-delayed_state.size(),RDS_carrier.end(),delayed_state.begin());	            //save state

    // PLL and NCO
    // *need state input*
    fmPllBlock(RDS_carrier_delayed, RDS_PLL_freq, params.if_Fs, pllState &state,  RDS_ncoScale);
    
    // *** RDS Demodulation ***
    // mixer
    std::vector<float> RDS_mixer_prod;
    mixer(RDS_carrier, RDS_channel, RDS_mixer_prod);
    
    // low pass filter
    std::vector<float> RDS_LPF_carrier_coeff;
    std::vector<float> RDS_LPF_state;
    std::vector<float> mixed_filtered_res;

    impulseResponseLPF(RDS_Fc, params.if_Fs, bpf_taps, RDS_LPF_carrier_coeff);
    fastConvWithState(mixed_filtered_res, RDS_mixer_prod, RDS_LPF_carrier_coeff, RDS_LPF_state, 1, 1, 1);

    // rational resampler
    rationalResampler(mixed_filtered_res, params.if_Fs, SPS_mode, length, RDS_U, output_signal, RDS_symbol_rate);

    // root-raised cosine filter

    // clock and data recovery


	// *** RDS Data Processing ***
    // Manchester and differenital decoding
    // frame synchronization and error detection
    // RDS application layer

	return 0;
}
