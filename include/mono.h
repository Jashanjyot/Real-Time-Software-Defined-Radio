#ifndef DY4_MONO_H
#define DY4_MONO_H

#include <vector>
#include "settings.h"
#include "rfFrontend.h"

const int mono_lpf_Fc_const = 16e3;      // the cutoff frequency to extract the IF from RF
const int mono_lpf_taps_const = 151; 	// the number of taps for the low-pass filter to extract the audio channel

class Mono {
    public:
        int lpf_gain;
        std::vector<float> state_mono_lpf;
        std::vector<float> processed_block;
        std::vector<float> lpf_coeff;
        std::vector<float> state;
        std::vector<float> delayed;    //mono_block phase delayed
        int lpf_Fc;
        int lpf_taps;

        Mono(Settings &params);
        void processing(Settings &params, std::vector<float>& fm_Demod);
};

#endif