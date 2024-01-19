#include "settings.h"
#include <iostream>

Settings::Settings(int mode, int processing_path){
    path = processing_path;
    exit = 0;
    switch(mode){
        case 0:
            rf_Fs = rf_Fs_mode0;
            if_Fs = if_Fs_mode0;
            audio_Fs = audio_Fs_mode0;
            break;
        case 1:
            rf_Fs = rf_Fs_mode1;
            if_Fs = if_Fs_mode1;
            audio_Fs = audio_Fs_mode1;
            break;
        case 2:
            rf_Fs = rf_Fs_mode2;
            if_Fs = if_Fs_mode2;
            audio_Fs = audio_Fs_mode2;
            break;
        case 3:
            rf_Fs = rf_Fs_mode3;
            if_Fs = if_Fs_mode3;
            audio_Fs = audio_Fs_mode3;
            break;
        default: break;
    }
    // RF-to-IF Resampling
    int factor = gcd(rf_Fs,if_Fs);
    rf_U = if_Fs/factor;
    rf_D = rf_Fs/factor;
    // IF-to-Audio Resampling
    factor = gcd(if_Fs,audio_Fs);
    if_U = audio_Fs/factor;
    if_D = if_Fs/factor;
    //select a block size that is a multiple of KB and a multiple of resample factors
    blockSize = 1024 * int(rf_D/rf_U) * int(if_D/if_U) * 2;//since 2 is a factor, it will always be even split between of IQ samples ^_^
};

int Settings::gcd(int a,int b) { //from https://stackoverflow.com/questions/41066488/how-the-following-function-of-finding-gcd-of-two-numbers-using-bitwise-operation
            while (b) {         // while b is non-zero
                int temp = b; // cache current value of b
                b = (a % b);  // set b equal to (a modulo b)
                a = temp;     // set a equal to the original value of b
            }
            return a; // once b is zero, return a
}