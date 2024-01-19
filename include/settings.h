/*
    From constraints document, the settings for group 18
						    Mode 0 	Mode 1	Mode 2	Mode 3
    RF Fs (Ksamples/s)	    2400	1152	2400	1920	#radio-frequency (RF) sampling rate
    IF Fs (Ksamples/s)	    240		288		240		320
    Audio Fs (Ksamples/s)	48		48		44.1	44.1
*/

#ifndef DY4_MODE_H
#define DY4_MODE_H

// Mode 0
const int rf_Fs_mode0 = 2400e3;
const int if_Fs_mode0 = 240e3;
const int audio_Fs_mode0 = 48e3;
// Mode 1
const int rf_Fs_mode1 = 1152e3;
const int if_Fs_mode1 = 288e3;
const int audio_Fs_mode1 = 48e3;
// Mode 2
const int rf_Fs_mode2 = 2400e3;
const int if_Fs_mode2 = 240e3;
const int audio_Fs_mode2 = 44.1e3;
// Mode 3
const int rf_Fs_mode3 = 1920e3;
const int if_Fs_mode3 = 320e3;
const int audio_Fs_mode3 = 44.1e3;

class Settings{
    public:
        int path;
        int rf_Fs;
        int if_Fs;
        int audio_Fs;
        int rf_U;
        int rf_D;
        int if_U;
        int if_D;
        float blockSize;
        int exit;

        Settings(int mode, int path);
        int gcd(int a,int b);
};

#endif