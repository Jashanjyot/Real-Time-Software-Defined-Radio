#include <vector>
#include "settings.h"

const float RDS_extract_Fb = 54e3;
const float RDS_extract_Fe = 60e3;
const float bpf_taps = 151;
const float RDS_carrier_Fb = 113.5e3;
const float RDS_carrier_Fe = 114.5e3;
const float RDS_PLL_freq = 114e3;
const float RDS_ncoScale = 0.5;
const float RDS_Fc = 3e3;

class RDS{
    public:

        RDS(Settings &params);

};