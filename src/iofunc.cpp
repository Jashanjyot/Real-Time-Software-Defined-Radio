/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "iofunc.h"

// some basic functions for printing information from vectors
// or to read from/write to binary files in 32-bit float format
void printRealVector(const std::vector<float> &x)
{
	std::cout << "Printing float vector of size " << x.size() << "\n";
	for (int i = 0; i < (int)x.size(); i++)
		std::cout << x[i] << " ";
	std::cout << "\n";
}

void printComplexVector(const std::vector<std::complex<float>> &X)
{
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (int i = 0; i < (int)X.size(); i++)
		std::cout << X[i] << " ";
	std::cout << "\n";
}

// assumes data in the raw binary file is in 32-bit float format
void readBinData(const std::string in_fname, std::vector<float> &bin_data)
{
	std::ifstream fdin(in_fname, std::ios::binary);
	if(!fdin) {
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cout << "Reading raw binary from \"" << in_fname << "\"\n";
	}
	fdin.seekg(0, std::ios::end);
	const unsigned int num_samples = fdin.tellg() / sizeof(float);

	bin_data.resize(num_samples);
	fdin.seekg(0, std::ios::beg);
	fdin.read(reinterpret_cast<char*>(&bin_data[0]), num_samples*sizeof(float));
	fdin.close();
}

void read8BitBinData(const std::string in_fname, std::vector<unsigned char> &bin_data)
{
	std::ifstream fdin(in_fname, std::ios::binary);
	if(!fdin) {
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cout << "Reading raw binary from \"" << in_fname << "\"\n";
	}
	fdin.seekg(0, std::ios::end);
	const unsigned int num_samples = fdin.tellg() / sizeof(unsigned char);

	bin_data.resize(num_samples);
	fdin.seekg(0, std::ios::beg);
	fdin.read(reinterpret_cast<char*>(&bin_data[0]), num_samples*sizeof(unsigned char));
	fdin.close();
}

// assumes data in the raw binary file is 32-bit float format
void writeBinData(const std::string out_fname, const std::vector<float> &bin_data)
{
	std::cout << "Writing raw binary to \"" << out_fname << "\"\n";
	std::ofstream fdout(out_fname, std::ios::binary);
	for (int i=0; i<(int)bin_data.size(); i++) {
		fdout.write(reinterpret_cast<const char*>(&bin_data[i]),\
								sizeof(bin_data[i]));
	}
	fdout.close();
}

void write16BitBinData(const std::string out_fname, const std::vector<short> &bin_data)
{
	std::cout << "Writing raw binary to \"" << out_fname << "\"\n";
	std::ofstream fdout(out_fname, std::ios::binary);
	for (int i=0; i<(int)bin_data.size(); i++) {
		fdout.write(reinterpret_cast<const char*>(&bin_data[i]),\
								sizeof(bin_data[i]));
	}
	fdout.close();
}

// data in the redirected input is interpreted as bytes
void readStdinBlockData(unsigned int num_samples,std::vector<char> &raw_data, std::vector<float> &block_data)
{
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));

	for(unsigned int k=0; k<num_samples; k++){
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0); //automatically normalizes the data to the range -1 to 1
	}
}

//Write stdout as a batch by block
void writeStdoutBlockData(std::vector<float> &processed_data,std::vector<short int> &audio_data)
{
	audio_data.resize(processed_data.size(),0);

	for (unsigned int k = 0; k<processed_data.size();k++){
		if(std::isnan(processed_data[k])) audio_data[k] = 0;
		//prepare a block of audio data to be redirected to stdout at once
		else audio_data[k] = static_cast<short int> (processed_data[k] * 32767);	//2^15 normalized floats between -1 and 1
	}
	// a block write is approx an order of magnitude faster than writing each sample
	//(assuming input streams are processed in blocks of hundreds of kilobytes)
	fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
}


void writeStdoutStereoAudioBlock(const std::vector<float> &audio_left, const std::vector<float> &audio_right)
{
	// file descriptor for the output to be written
	if (audio_left.size() != audio_right.size()) {
		std::cout << "Something got messed up with audio channels\n";
		std::cout << "They must have the same size ... exiting\n";
		exit(1);
	} else {
		std::cout << "Writing raw audio to stdout\n";
	}
	
	std::vector<short int> audio_interlaced(audio_left.size()*2,0);

	for (unsigned int k = 0; k<audio_interlaced.size();k++){
		if(k%2==0)
			if (!std::isnan(audio_left[k/2]))	audio_interlaced[k] = static_cast<short int> (audio_left[k/2]* 32767);
			else audio_interlaced[k] = 0;
		else
			if (!std::isnan(audio_right[(k-1)/2])) audio_interlaced[k] = static_cast<short int> (audio_right[(k-1)/2]* 32767);
			else audio_interlaced[k] = 0;
		
	}
	fwrite(&audio_interlaced[0], sizeof(short int), audio_interlaced.size(), stdout); //write L+R in one block, aplay can not take separate writes!
}

int getMode(int argc, char* argv[]){
	int mode = 0;

	if (argc < 2) {
		std::cerr << "Operating in default mode " << std::endl;
	} else if (argc == 2) {
		mode = atoi(argv[1]);
		if (mode> 3) {
		    std::cerr << "Wrong mode " << mode << std::endl; 
		    exit(1);
		}
	} else {
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or" << std::endl;
		std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl; 
		std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl; 
		exit(1);
	}
	std::cerr << "Operating in mode" << mode << std::endl;
	return mode;
}