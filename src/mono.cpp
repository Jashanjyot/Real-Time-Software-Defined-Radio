#include "filter.h"
#include "settings.h"
#include "mono.h"
#include "rfFrontend.h"

Mono::Mono(Settings &params){
    lpf_gain = params.if_U;
    lpf_taps = mono_lpf_taps_const*lpf_gain;
    state_mono_lpf.clear(); state_mono_lpf.resize(lpf_taps-1);
    lpf_Fc = mono_lpf_Fc_const;

	impulseResponseLPF(lpf_Fc, params.if_Fs*lpf_gain, lpf_taps, lpf_coeff);
}

void Mono::processing(Settings &params, std::vector<float>& fm_Demod){
    fastConvWithState(processed_block, fm_Demod, lpf_coeff, state_mono_lpf, params.if_U, params.if_D,lpf_gain);
}