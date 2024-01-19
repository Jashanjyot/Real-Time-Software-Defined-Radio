/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/
#include <iostream>
#include <vector>
#include <math.h>

#include "dy4.h"
#include "filter.h"

// function to compute the impulse response "h" based on the sinc function
// finite impulse response low pass filter
void impulseResponseLPF(float Fc, float Fs, unsigned short int N_taps, std::vector<float> &h)
{
	h.clear(); h.resize(N_taps, 0.0);
	float Normco= Fc/(Fs/2);

	for (int i = 0; i < N_taps; i++){
		if (i == ((N_taps-1)/2)){
			h[i] = Normco;
		}
		else{
			h[i] = Normco*(std::sin(PI* Normco*(i-(N_taps-1)/2))/(PI*Normco*(i-(N_taps-1)/2)));
		}

		h[i] = h[i]*(std::pow(std::sin(i*PI/N_taps),2));
	}
}

void bpf(float Fb,float Fe, float Fs, unsigned short int N_taps,  std::vector<float> &h)
{
	h.clear(); h.resize(N_taps, 0.0);
	float Norm_c= ((Fe+Fb)/2)/(Fs/2); //norm centre
	float Norm_p= (Fe-Fb)/(Fs/2); //norm pass

	for (int i = 0; i < N_taps; i++){
		if (i == ((N_taps-1)/2)){
			h[i] = Norm_p;
		}
		else{
			h[i] = Norm_p*(std::sin(PI*(Norm_p/2)*(i-(N_taps-1)/2))/(PI*(Norm_p/2)*(i-(N_taps-1)/2)));
		}
		h[i] = h[i]*std::cos(i*PI*Norm_c);
		h[i] = h[i]*std::pow(std::sin(i*PI/N_taps),2);
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// bring your own functionality
  // allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);
	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for(unsigned int n=0; n<y.size();n++){
		for(unsigned int k=0; k<h.size();k++){
			if((n-k)>=0 && (n-k)<x.size()){
				y[n]+=h[k]*x[n-k]; //from lecture weel2-1 pg15
			}
		}
	}
}

void fastConvWithState(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h,
						std::vector<float> &state, int U, int D, int gain)
{
	y.clear(); y.resize(x.size()*U/D);
	for(int n = 0; n<int(y.size()); n++){
		for (int s=(n*D)%U; s<int(h.size()); s+=U){
			int x_idx = int(n*D-s)/U;
			if (x_idx >= 0 && x_idx <int(x.size())*U){
				y[n] += h[s] * x[x_idx] * gain;
			}
			else{
				y[n] += h[s] * state[int(state.size())+x_idx] * gain;
			}
    	}
	}

	std::copy(x.end()-h.size()+1,x.end(),state.begin());
}
