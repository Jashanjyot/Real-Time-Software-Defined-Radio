/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_IOFUNC_H
#define DY4_IOFUNC_H

// add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>

// declaration of a function prototypes
void printRealVector(const std::vector<float> &);

void printComplexVector(const std::vector<std::complex<float>> &);

void readBinData(const std::string, std::vector<float> &);

void read8BitBinData(const std::string, std::vector<float> &);

void writeBinData(const std::string, const std::vector<float> &);

void readStdinBlockData(unsigned int num_samples,std::vector<char> &raw_data, std::vector<float> &block_data);

void writeStdoutBlockData(std::vector<float> &processed_data,std::vector<short int> &audio_data);

int getMode(int argc, char* argv[]);

void writeStdoutStereoAudioBlock(const std::vector<float> &audio_left, const std::vector<float> &audio_right);


#endif // DY4_IOFUNC_H
