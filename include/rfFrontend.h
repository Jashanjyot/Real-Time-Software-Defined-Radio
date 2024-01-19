#ifndef DY4_RFFRONTEND_H
#define DY4_RFFRONTEND_H

#include <vector>

const int rf_Fc_const = 100e3;	    // the cutoff frequency to extract the FM channel from raw IQ data
const int rf_lpf_taps_const = 151; 	// the number of taps for the low-pass filter to extract the FM channel

class rfFrontend{
    public:
        float blockSize;
        int lpf_taps;
        int lpf_Fc;
        int lpf_gain;
        std::vector<char> raw_data8;		//8-bit unsigned integer IQ data
        std::vector<float> raw_data32;		//32-bit float IQ data
        std::vector<float> lpf_coeff;
        std::vector<float> state_i_lpf;		// add state for LPF of RF Processing
        std::vector<float> state_q_lpf;
        float prev_I;
        float prev_Q;
        std::vector<float> i_samples; 
        std::vector<float> q_samples;
        std::vector<float> i_processed;
        std::vector<float> q_processed;
        std::vector<float> fm_Demod;

        rfFrontend(Settings &params);
        void getIQSamples();
        void processing(Settings &params);
        void fmDemodDeriv(std::vector<float> &I,std::vector<float> &Q,
                                    float &prev_I,float &prev_Q, std::vector<float> &fmDemod);
};

#endif