/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>
#include <math.h>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void bpf(float Fb,float Fe, float Fs, unsigned short int N_taps,  std::vector<float> &h);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void fastConvWithState(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, 
                        std::vector<float> &state, int U, int D, int gain);

#endif // DY4_FILTER_H
