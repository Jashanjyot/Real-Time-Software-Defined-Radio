/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"

// just DFT function (no FFT yet)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	Xf.clear(); Xf.resize(x.size(), static_cast<std::complex<float>>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
	// only the positive frequencies
	Xmag.clear(); Xmag.resize(Xf.size(), static_cast<float>(0));
	for (int i = 0; i < (int)Xf.size(); i++) {
		Xmag[i] = std::abs(Xf[i])/Xf.size();
	}
}

// add your own code to estimate the PSD

//////////////////////////////////////////////////////

// added IDFT

void IDFT(const std::vector<std::complex<float>> &Xf, std::vector<std::complex<float>> &x) {
	x.clear(); x.resize(Xf.size(), static_cast<std::complex<float>>(0));
	for (int k = 0; k < (int)x.size(); k++) {
		for (int m = 0; m < (int)x.size(); m++) {
			std::complex<float> expval(0, 2*PI*(k*m) / Xf.size());
			x[k] += Xf[m] * std::exp(expval);
		}
		x[k] /= Xf.size();
	}
}

// added FFT

// this is a support function used for bit reversal
unsigned int swap_bits(unsigned int x, unsigned char i, unsigned char j) {

	unsigned char bit_i = (x >> i) & 0x1L;
	unsigned char bit_j = (x >> j) & 0x1L;

	unsigned int val = x;
	val = bit_i ? (val | (0x1L << j)) : (val & ~(0x1L << j));
	val = bit_j ? (val | (0x1L << i)) : (val & ~(0x1L << i));

	return val;
}

// bit reversal is needed for the non-recursive (forward traversal) version of FFT
unsigned int bit_reversal(unsigned int x, unsigned char bit_size) {

	unsigned int val = x;

	for (int i=0; i < int(bit_size/2); i++)
		val = swap_bits(val, i, bit_size-1-i);

	return val;
}

// precomputation of twiddle factors (needed for faster versions of FFT)
void compute_twiddles(std::vector<std::complex<float>> &twiddles) {
	for (int k=0; k<(int)twiddles.size(); k++) {
			std::complex<float> expval(0.0, -2*PI*float(k)/ NFFT);
			twiddles[k] = std::exp(expval);
	}
}

// recursive version of FFT
// based on first principles from the mathematical derivation
void FFT_recursive(const std::vector<std::complex<float>> &x, \
	std::vector<std::complex<float>> &Xf) {

	if (x.size() > 1) {
		// declare vectors and allocate space for the even and odd halves
		std::vector<std::complex<float>> xe(int(x.size()/2)), xo(int(x.size()/2));
		std::vector<std::complex<float>> Xfe(int(x.size()/2)), Xfo(int(x.size()/2));

		// split into even and odd halves
		for (int k=0; k<(int)x.size(); k++)
			if ((k%2) == 0) xe[k/2] = x[k];
			else xo[k/2] = x[k];

		// call recursively FFT of half size for even and odd halves respectively
		FFT_recursive(xe, Xfe);
		FFT_recursive(xo, Xfo);

		// merge the results from the odd/even FFTs (each of half the size)
		for (int k=0; k<(int)xe.size(); k++) {
				std::complex<float> expval(0.0, -2*PI*float(k)/ x.size());
				std::complex<float> twiddle = std::exp(expval);
				Xf[k]           = Xfe[k] + twiddle * Xfo[k];
				Xf[k+xe.size()] = Xfe[k] - twiddle * Xfo[k];
		}
	} else {
		// end of recursion - copy time domain samples to frequency bins (default values)
		Xf[0] = x[0];
	}
}

// improved version of FFT
// the speed-up is mainly from avoiding the recomputation of twiddle factors
// the math library functions do not need to be recalled at each level of recursion
void FFT_improved(const std::vector<std::complex<float>> &x, \
	std::vector<std::complex<float>> &Xf, \
	const std::vector<std::complex<float>> &twiddles, \
	const unsigned char recursion_level) {

	if (x.size() > 1) {
		int half_size = int(x.size()/2);
		std::vector<std::complex<float>> xe(half_size), xo(half_size);
		std::vector<std::complex<float>> Xfe(half_size), Xfo(half_size);

		for (int k=0; k<half_size; k++) {
			xe[k] = x[k*2];
			xo[k] = x[k*2+1];
		}

		FFT_improved(xe, Xfe, twiddles, recursion_level+1);
		FFT_improved(xo, Xfo, twiddles, recursion_level+1);

		for (int k=0; k<half_size; k++) {
				// the precomputed twiddle factors are reused at each recursion level
				Xf[k]           = Xfe[k] + twiddles[k*(1<<(recursion_level-1))] * Xfo[k];
				Xf[k+half_size] = Xfe[k] - twiddles[k*(1<<(recursion_level-1))] * Xfo[k];
		}
	} else {
		Xf[0] = x[0];
	}
}

// optimized version of FFT is based on non-recursive forward traversal
// speed-up comes (in-part) from avoiding loading/unloading the stack with recursive calls
// however, in addition to pre-computation of twiddle factors, the central reason for
// speed-up is the in-place computation of the intermediate vectors at each FFT level
// the core idea is that ONLY two input values are used to compute two output levels
void FFT_optimized(const std::vector<std::complex<float>> &x, \
	std::vector<std::complex<float>> &Xf, \
	const std::vector<std::complex<float>> &twiddles) {

	unsigned char no_levels = (unsigned char)std::log2((float)x.size());
	for (int i=0; i<(int)x.size(); i++) {
		// time-domain values indexed in a bit-reversed manner before
		// copied to the leftmost (first) level from where forward traversal starts
		Xf[i] = x[bit_reversal(i, no_levels)];
	}

	unsigned int step_size = 1;

	std::complex<float> tmp;
	for (unsigned char l=0; l<no_levels; l++) {
		for (unsigned int p=0; p<x.size(); p+=2*step_size) {
			for (unsigned int k=p; k<p+step_size; k++) {
				// in-place computation - only two output values depend on two input values
				// indexing of twiddle factors exploits the fact that
				// index k for FFT of size n is the same as index 2k for FFT of size 2n
				tmp             = Xf[k] + twiddles[(k-p)*(1<<(no_levels-1-l))] * Xf[k+step_size];
				Xf[k+step_size] = Xf[k] - twiddles[(k-p)*(1<<(no_levels-1-l))] * Xf[k+step_size];
				Xf[k]           = tmp;
			}
		}
		step_size *= 2;
	}
}

// add your own code to estimate the PSD
void estimatePSD(std::vector<float> &samples, const int NFFTval, float Fs, std::vector<float> &freq, std::vector<float> &psd_est){

	int freq_bins = NFFTval; // rename the NFFT argument

	float df = Fs/freq_bins; // frequency increment (or resolution of the frequency bins)

	freq.resize(Fs/2); // come back

	// fill in the freq array
	for(unsigned int i=0; i< freq.size(); i++){
		freq[i] = i*df; // starting at 0, every index of freq will increment by df
	}

	// design the Hann window used to smooten the discrete data in order
	// to reduce the spectral leakage after the Fourier Transform
	std::vector<float> hann;
	hann.clear(); hann.resize(freq_bins);

	for(unsigned int i=0; i<hann.size(); i++){
		hann[i] = pow(sin(i*PI/freq_bins),2);
	}

	// create an empty vector where the PSD for each segment is computed
	std::vector<float> psd_list;
	psd_list.clear();

	// samples should be a multiple of frequency bins, so the
	// number of segments used for estimation is an integer
	int no_segments = (int)floor(samples.size()/(float)freq_bins);

	// iterate through all segments
	for(int k=0; k<no_segments; k++){
		// apply the hann window (using pointwise multiplication)
		// before computing the Fourier transform on a segment
		std::vector<float> windowed_samples;
		windowed_samples.clear(); windowed_samples.resize(freq_bins);

		// need to fill in each index of windowed_samples
		for(int j=0; j<freq_bins; j++){
			windowed_samples[j] = samples[k*freq_bins+j] * hann[j];
		}

		std::vector<std::complex<float>> Xf(freq_bins);

		DFT(windowed_samples, Xf);

		// since input is real, we keep only the positive half of the spectrum
		// however, we will also add the signal energy of negative frequencies
		// to have a better a more accurate PSD estimate when plotting

		std::vector<std::complex<float>> positive_Xf(freq_bins/2); // define new vector to kep only positive freq bins

		for (int i=0; i<freq_bins/2; i++){
			positive_Xf[i] = Xf[i]; // keep only positive freq bins
		}

		std::vector<float> psd_seg;
		psd_seg.clear(); psd_seg.resize(freq_bins/2);
		std::vector<float> psd_list;
		psd_list.clear(); psd_seg.resize(freq_bins/2);

		// fill in values for psd_seg vector
		for(int i=0; i<freq_bins/2; i++){
			psd_seg[i] = (1/(Fs*freq_bins/2))*pow(abs(positive_Xf[i]),2); // compute signal power
			psd_seg[i] = 2*psd_seg[i]; // add the energy from the negative freq bins
		}

		// translate to the decibel (dB) scale
		for(unsigned int i=0; i<psd_seg.size(); i++){
			psd_seg[i] = 10*log10(psd_seg[i]);
			// append to the list where PSD for each segment is stored
			// in sequential order
			psd_list.insert(psd_list.end(), psd_seg.begin(), psd_seg.end());
		}

		// compute the estimate to be returned by the function through averaging
		psd_est.resize(freq_bins/2);

		// iterate through all the frequency bins (positive freq only)
		// from all segments and average them (one bin at a time ...)
		for (int k=0; k<freq_bins/2; k++){
			// iterate through all the segments
			for (int i=0; i<no_segments; i++){
				psd_est[k] += psd_list[k + i*(int)(freq_bins/2)];
			}
			// compute the estimate for each bin
			psd_est[k] = psd_est[k] / no_segments;
		}
	}

}