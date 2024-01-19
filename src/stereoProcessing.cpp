#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "dy4.h"
#include "fourier.h"
#include "filter.h"
#include "stereoProcessing.h"
#include "mono.h"
#include "rfFrontend.h"

// From constraints document, the settings for group 18
// 						Mode 0 	Mode 1	Mode 2	Mode 3
// RF Fs (Ksamples/s)	2400	1152	2400	1920	#radio-frequency (RF) sampling rate
// IF Fs (Ksamples/s)	240		288		240		320
// Audio Fs (Ksamples/s)	48		48		44.1	44.1

float rf_Fc = 100e3;	// the cutoff frequency to extract the FM channel from raw IQ data
float if_Fc = 16e3;
int rf_taps = 151; 	// the number of taps for the low-pass filter to extract the FM channel

// Mode 0
float rf_Fs0 = 2400e3;
float if_Fs0 = 240e3;
int rf_U0 = 1;
int rf_D0 = 10;
int if_U0 = 1;
int if_D0 = 5;
float audio_Fs0 = 48e3;

// Mode 1
int up1 = 1;
int down1 = 6;
int gain1 = 1;

// Mode 2
int up2 = 147;
int down2 = 800;
int gain2 = 147;

Stereo::Stereo(Settings &params){
    Fb_LR = 22e3;
    Fe_LR = 54e3;
    Fb_pilot = 18.5e3;
    Fe_pilot = 19.5e3;
    Fm = 15e3;
    lpf_gain = params.if_U;
    bpf_taps = 151;
    lpf_taps = 151*lpf_gain;
    state_LR_bpf.clear(); state_LR_bpf.resize(bpf_taps-1);
    state_pilot_bpf.clear(); state_pilot_bpf.resize(bpf_taps-1);
    pllFreq = 19e3;
    state_dig_filt_lpf.clear(); state_dig_filt_lpf.resize(lpf_taps);
    ncoScale = 2.0;
}

void Stereo::fmPllBlock(std::vector<float>& pllIn, float freq, float Fs, pllState &state, float ncoScale = 2.0,
                        float phaseAdjust = 0.0, float normBandwidth = 0.01) 
{
    float Cp = 2.666;
    float Ci = 3.555;

    float Kp = normBandwidth * Cp;
    float Ki = (normBandwidth * normBandwidth) * Ci;

    state.ncoOut.resize(pllIn.size() + 1);
    state.ncoOut[0] = state.ncoOut0;

    for (unsigned int k = 0; k < pllIn.size(); k++) {

        // phase detector
        float errorI = pllIn[k] * (+state.feedbackI); // complex conjugate of the
        float errorQ = pllIn[k] * (-state.feedbackQ); // feedback complex exponential

        // four-quadrant arctangent discriminator for phase error detection
        float errorD = std::atan2(errorQ, errorI);

        // loop filter
        state.integrator = state.integrator + Ki * errorD;

        // update phase estimate
        state.phaseEst = state.phaseEst + Kp * errorD + state.integrator;

        // internal oscillator
        ++state.trigOffset;
        float trigArg = 2 * PI * (freq / Fs) * (state.trigOffset) + state.phaseEst;
        state.feedbackI = std::cos(trigArg);
        state.feedbackQ = std::sin(trigArg);
        state.ncoOut[k + 1] = std::cos(trigArg * ncoScale + phaseAdjust);
    }
    state.ncoOut0 = state.ncoOut[state.ncoOut.size()-1];
}

void Stereo::initPll(){
    //init state
    state.integrator = 0.0;
    state.phaseEst = 0.0;
    state.feedbackI = 1.0;
    state.feedbackQ = 0.0;
    state.ncoOut.resize(1);
    state.ncoOut[0] = 1.0;
    state.trigOffset = 0;
    state.ncoOut0 = 1.0;
}

void Stereo::mixer(std::vector<float> &v1, std::vector<float> &v2, std::vector<float> &product)
{
    product.clear(); product.resize(v2.size());

    for(unsigned int i=0;i<v1.size();i++){
        product[i] = v1[i] * v2[i] * 2; //since cos(α)cos(𝛽) = (cos(α + 𝛽) + cos(α − 𝛽) ) / 2, need to multiply by 2
    }
}

void Stereo::combiner(std::vector<float> &mono_data,std::vector<float> &stereo_data,
                        std::vector<float> &left_channel,std::vector<float> &right_channel){
    left_channel.clear(); left_channel.resize(stereo_data.size());
    right_channel.clear(); right_channel.resize(stereo_data.size());

    for(unsigned int k=0; k < std::min(stereo_data.size(), mono_data.size()); k++){
        left_channel[k] = stereo_data[k] + mono_data[k];
        right_channel[k] = stereo_data[k] - mono_data[k];
    }
}

void Stereo::processing(Settings &params, std::vector<float>& fm_Demod, Mono mono){
    // * Stereo Carrier Recovery
    // extract the 19kHz pilot tone through BPF
    bpf(Fb_pilot, Fe_pilot, params.if_Fs, bpf_taps, bpf_pilot);
    fastConvWithState(carrier, fm_Demod, bpf_pilot, state_pilot_bpf,1,1,1); //U and D = 1 1 since no resampling needed
    
    // then synchronize to it using a PLL and multiply the output of PLL by ncoScale=2
    initPll();
    fmPllBlock(carrier, pllFreq, params.if_Fs, state, ncoScale,0,0.01);
    // carrier goes through pll and outputs state.ncoOut

    // * Stereo Channel Extraction
    // apply BPF
    bpf(Fb_LR, Fe_LR, params.if_Fs, bpf_taps, bpf_LR);
    fastConvWithState(channel, fm_Demod, bpf_LR, state_LR_bpf,1,1,1);

    // * Stereo Processing
    // Mixer: multiply the output signal from stereo carrier recovery with bpf signal?
    mixer(state.ncoOut, channel, mixerProduct);
    // Digital filtering (LPF, sample rate conversion)
    impulseResponseLPF(Fm, params.if_Fs*lpf_gain, lpf_taps, lpf_dig_filt);
    fastConvWithState(stereo_block, mixerProduct, lpf_dig_filt,state_dig_filt_lpf, params.if_U, params.if_D, lpf_gain);
    
    // combine stereo data with mono audio to get left and right audio channels
    // pad with zeros to account for delay
    mono.state.resize((bpf_taps-1)/2*params.if_U/params.if_D);
    mono.delayed.resize(mono.processed_block.size());
    std::copy(mono.state.begin(),mono.state.end(),mono.delayed.begin());	//copy in state of old mono block
    std::copy(mono.processed_block.begin(),mono.processed_block.end()-mono.state.size(),mono.delayed.begin()+mono.state.size());	//copy relevant part of mono block
    std::copy(mono.processed_block.end()-mono.state.size(),mono.processed_block.end(),mono.state.begin());	//save mono state

    // combine stereo data with mono audio to get left and right audio channels
    combiner(mono.delayed, stereo_block, left_channel, right_channel);

}


