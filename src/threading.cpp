#include <vector>
#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>

#include "threading.h"

// for random value generation:
// #define RANDOM_BOUND 1000
// #define RANDOM_SEED 0 

//define variable to track when to terminate the threads
#define TOTAL_ELEMS 10
int generated_elems =0; //for illustrative purposed
//change to our generated outputs?
//size of queue:
#define QUEUE_ELEMS 3 //can be changed, its up to us

//producer thread: RF_thread ->  RF_front_end_thread
void rf_thread(std::vector<int> &my_queue, std::atomic<int>& write_offset,std::atomic<int>& read_offset){
    //thread runs until a custom end condition
    //using example 4 from lecture 10-2:
    impulseResponseLPF(lpf_Fc, params.rf_Fs*lpf_gain, lpf_taps, lpf_coeff); //
    while(generated_elems < TOTAL_ELEMS){ //change while loop
        //produce new elements
        std::vector<int>elem;
        my_generate(elem);
        generated_elems++;

        rfFrontend()

        while(write_offset.load()>= (read_offset.load()+QUEUE_ELEMS)){
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        std::vector<int>::difference_type address_offset = (write_offset.load()% QUEUE_ELEMS)*ELEM_SIZE;
        //COPY DATA TO QUEUE FROM PROCESSED ELEMENT:
        std::copy_n(elem.begin(), elem.size(), my_queue.begin()+address_offset);
        //increment the write offset
        write_offset.fetch_add(1);

        }

    }
}
//END


//consumers week 10-2 example 4:
void audio_thread(std::vector<int> &my_queue, std::atomic<int>& write_offset,std::atomic<int>& read_offset){
    std::vector<int> elem(ELEM_SIZE);
    while(1){
        while(write_ofset.load() <= read_offset.load()){
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        std::vector<int>::difference_type address_offset = (read_offset.load()% QUEUE_ELEMS)*ELEM_SIZE;
        std::copy_n(my_queue.begin()+address_offset, elem.size(), elem.begin());
        read_offset.fetch_add(1);

        my_compute(elem);
        if((generated_elems == TOTAL_ELEMS) && (wite_offset.load() <= read_offset.load())){
            break;
        }
    }
}

void RDS_thread(std::vector<int> &my_queue, std::atomic<int>& write_offset,std::atomic<int>& read_offset){
    std::vector<int> elem(ELEM_SIZE);
    while(1){
        while(write_ofset.load() <= read_offset.load()){
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        std::vector<int>::difference_type address_offset = (read_offset.load()% QUEUE_ELEMS)*ELEM_SIZE;
        std::copy_n(my_queue.begin()+address_offset, elem.size(), elem.begin());
        read_offset.fetch_add(1);

        my_compute(elem);
        if((generated_elems == TOTAL_ELEMS) && (wite_offset.load() <= read_offset.load())){
            break;
        }
    }
}

