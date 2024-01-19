#include "filter.h"
#include "settings.h"
#include "rfFrontend.h"
#include <cmath>


rfFrontend::rfFrontend(Settings &params){
    // Init rf variables
    blockSize = params.blockSize;
    lpf_gain = params.rf_U;
    lpf_taps = rf_lpf_taps_const*lpf_gain;
    lpf_Fc = rf_Fc_const;
    raw_data8.clear(); raw_data8.resize(blockSize);
	raw_data32.clear(); raw_data32.resize(blockSize);
    state_i_lpf.clear(); state_i_lpf.resize(lpf_taps-1);
	state_q_lpf.clear(); state_q_lpf.resize(lpf_taps-1);
    prev_I = 0;
	prev_Q = 0;
    i_samples.clear(); i_samples.resize(blockSize/2);
	q_samples.clear(); q_samples.resize(blockSize/2);
    
    //create LPF
	impulseResponseLPF(lpf_Fc, params.rf_Fs*lpf_gain, lpf_taps, lpf_coeff);
}

void rfFrontend::getIQSamples(){
    // retrieving i and q samples, which are even and odd indices, respectively of audio data
    for (int k=0; k<blockSize;k+=2){
			i_samples[k/2] = raw_data32[k];
			q_samples[k/2] = raw_data32[k+1];
		}
}

void rfFrontend::processing(Settings &params){
    fastConvWithState(i_processed, i_samples, lpf_coeff, state_i_lpf, params.rf_U, params.rf_D,lpf_gain);
    fastConvWithState(q_processed, q_samples, lpf_coeff, state_q_lpf, params.rf_U, params.rf_D,lpf_gain);
    fmDemodDeriv(i_processed, q_processed, prev_I, prev_Q, fm_Demod);
}

void rfFrontend::fmDemodDeriv(std::vector<float> &I,std::vector<float> &Q,float &prev_I,float &prev_Q, std::vector<float> &fmDemod)
{
	fmDemod.clear();
	fmDemod.resize(I.size(), 0.0);

	for (unsigned int k = 0; k < I.size(); k++){
		fmDemod[k] = ((I[k] * (Q[k] - prev_Q)) - (Q[k] * (I[k] - prev_I))) / (std::pow(I[k], 2) + std::pow(Q[k], 2));
        if(std::isnan(fmDemod[k]) || std::isinf(fmDemod[k])) fmDemod[k] = 0;
		prev_Q = Q[k];
		prev_I = I[k];
	}
}