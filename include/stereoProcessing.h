#include <vector>
#include "settings.h"
#include "rfFrontend.h"
#include "mono.h"

struct pllState {
            std::vector<float> ncoOut;
            float integrator;
            float phaseEst;
            float feedbackI;
            float feedbackQ;
            int trigOffset;
            float ncoOut0;
};

class Stereo{
    public:
        pllState state;
        float Fb_LR;
        float Fe_LR;
        float Fb_pilot;
        float Fe_pilot;
        int Fm;
        unsigned short int bpf_taps;
        int lpf_taps;
        int lpf_gain;
        std::vector<float> bpf_LR;
        std::vector<float> bpf_pilot;
        std::vector<float> state_LR_bpf;
        std::vector<float> state_pilot_bpf;
        float pllFreq;
        std::vector<float> lpf_dig_filt;
        std::vector<float> state_dig_filt_lpf;
        std::vector<float> carrier;
        std::vector<float> channel;
        std::vector<float> mixerProduct;
        std::vector<float> stereo_block;
        std::vector<float> left_channel;
        std::vector<float> right_channel;
        float ncoScale;

        Stereo(Settings &params);
        void processing(Settings &params, std::vector<float>& fm_Demod, Mono mono);
        void initPll();
        void mixer(std::vector<float> &v1, std::vector<float> &v2, std::vector<float> &product);
        void fmPllBlock(std::vector<float>& pllIn, float freq, float Fs, pllState &state, float ncoScale, 
                        float phaseAdjust, float normBandwidth);
        void combiner(std::vector<float> &mono_data,std::vector<float> &stereo_data,
                        std::vector<float> &left_channel,std::vector<float> &right_channel);
};