/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "settings.h"
#include "project.h"
#include "mono.h"
#include "rfFrontend.h"
#include "stereoProcessing.h"

#include <atomic>
#include <functional>
#include <thread>
#include <vector>
#include <iostream>
#include <algorithm>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <utility>

#define QUEUE_ELEMS 10 // need to change

// Producer thread: rf_thread
void rf_thread(std::vector<float> &my_queue, std::atomic<int>& write_offset, std::atomic<int>& audio_read_offset, Settings& params, int ELEM_SIZE){
	
	rfFrontend rfFe(params);

	for (int block_id=0; ; block_id++) {	//loop through stdin
		readStdinBlockData(params.blockSize, rfFe.raw_data8, rfFe.raw_data32);
		if ((std::cin.rdstate()) != 0) {
			std::cerr << "End of input stream reached" << std::endl;
			break;
		}
		// std::cerr << "Read block " << block_id << std::endl;
		rfFe.getIQSamples();
		// *** RF Front-End Processing ***
		rfFe.processing(params);
		
		// check if the queue is full
        while(write_offset.load()>= (audio_read_offset.load()+QUEUE_ELEMS)){
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
    	}
		// compute the address offset from beginning of the queue
		std::vector<float>::difference_type address_offset = (write_offset.load()% QUEUE_ELEMS)*ELEM_SIZE;
		
		//COPY DATA TO QUEUE FROM rfFrontend fm_Demod:
		std::copy_n(rfFe.fm_Demod.begin(), rfFe.fm_Demod.size(), my_queue.begin()+address_offset);
		
		//increment the write offset
		write_offset.fetch_add(1);
		
		// std::cerr << "write_offset: " << write_offset.load()<< std::endl;

	}
	params.exit = 1;
}

// Consumer thread: audio_thread
void audio_thread(std::vector<float> &my_queue, std::atomic<int>& write_offset,std::atomic<int>& audio_read_offset, Settings& params, int ELEM_SIZE){
	
	//init variables
	std::vector<float> elem(ELEM_SIZE);
	std::vector<short int> audio_block; 
	Mono mono(params);
	Stereo stereo(params);

    while(1){
		// check if the queue is empty
        while(write_offset.load() <= audio_read_offset.load()){
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
		// std::cerr << "write_offset: " << write_offset.load()<< std::endl;
		// std::cerr << "read offset: " << audio_read_offset.load()<< std::endl;
        }
		// compute the address offset from the beginning of the queue
        std::vector<float>::difference_type address_offset = (audio_read_offset.load()% QUEUE_ELEMS)*ELEM_SIZE;
        
		// copy the data FROM the queue to the extracted element
		std::copy_n(my_queue.begin()+address_offset, elem.size(), elem.begin());
        
		// increment the read offset
		audio_read_offset.fetch_add(1);
		
		// std::cerr << "read offset: " << audio_read_offset.load()<< std::endl;

		// *** Mono Processing ***
		mono.processing(params, elem);
		if(params.path == 0)//Mono Audio
			writeStdoutBlockData(mono.processed_block,audio_block); 
		
		else if(params.path == 1){	//Stereo Audio
			// *** Stereo Processing ***
			stereo.processing(params, elem, mono);

			// // Write to file for stereo audio
			writeStdoutStereoAudioBlock(stereo.left_channel, stereo.right_channel);
		}
		if(params.exit) break;
	}
}

int main(int argc, char* argv[])
{
	int mode = 0;
	int path = 0;
	if(argc < 3){
		std::cerr << "Invalid input. The format is ./project <mode> <path>" << std::endl;
		std::cerr << "Possible modes range from 0 to 3" << std::endl;
		std::cerr << "Possible paths are: 0:Mono, 1:Stereo" << std::endl;
		std::cerr << "Using default mode." << std::endl;
	} else {
		mode = atoi(argv[1]);		//get mode (0, 1, 2, 3)
		path = atoi(argv[2]);		//get processing path eg) 0-mono, 1-stereo
	}
	Settings params(mode,path);

	int ELEM_SIZE = params.blockSize/2*params.rf_U/params.rf_D;	//size of fmDemod

	std::vector<float> my_queue(QUEUE_ELEMS*ELEM_SIZE);

	std::atomic<int> write_offset (0);
	std::atomic<int> audio_read_offset (0);
	//std::atomic<int> rds_read_offset = 0; 

	// RF thread - producer
	std::thread rf_prod_thread = std::thread (rf_thread, std::ref(my_queue), std::ref(write_offset), std::ref(audio_read_offset), std::ref(params), std::ref(ELEM_SIZE));

	// audio thread (mono & stereo) - consumer
	std::thread audio_cons_thread = std::thread (audio_thread, std::ref(my_queue),std::ref(write_offset), std::ref(audio_read_offset), std::ref(params), std::ref(ELEM_SIZE));

	// RDS thread - consumer
	//std::thread RDS_cons_thread = std::thread(RDS_thread, std::ref(fm_Demod), std::ref(write_offset), std::ref(rds_read_offset));

	rf_prod_thread.join();
	audio_cons_thread.join();
	//RDS_cons_thread.join();

	return 0;
}
