#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

void readStdinBlockData(unsigned int num_samples, std::vector<float> &block_data)
{
	std::vector<char> raw_data(num_samples);

	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for(int k=0; k<(int)num_samples; k++){
		//automatically normalizes the data to the range -1 to 1
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}
}

int main(){
    std::vector<float> audio_data;
    unsigned int BLOCK_SIZE  = 1024;
	readStdinBlockData (BLOCK_SIZE, audio_data);
    return 0;
}