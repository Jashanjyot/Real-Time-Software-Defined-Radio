#include <vector>
#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>

#include "settings.h"
#include "rfFrontend.h"
#include "mono.h"

class Threading{
    public:
        Threading(Settings &params);

        void rf_thread(std::vector<int> &my_queue, std::atomic<int>& write_offset,std::atomic<int>& read_offset);    
        void audio_thread(std::vector<int> &my_queue, std::atomic<int>& write_offset,std::atomic<int>& read_offset);
        void RDS_thread(std::vector<int> &my_queue, std::atomic<int>& write_offset,std::atomic<int>& read_offset);

};