/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_DY4_H
#define DY4_DY4_H

// some general and reusable stuff
// our beloved PI constant
#define PI 3.14159265358979323846

// although we use DFT (no FFT ... yet), the number of points for a
// Fourier transform is defined as NFFT (same as matplotlib)
#define NFFT 512

#include <mutex>
#include <condition_variable>
#include <queue>

//Struct for PLL
struct PllState {
    float integrator = 0.0f;
    float phaseEst = 0.0f;
    float feedbackI = 1.0f;
    float feedbackQ = 0.0f;
    int trigOffset = 0;
};

//Struct for MUlti-Threading
struct SharedQueue {
    std::queue<std::vector<float>> queue;
    std::mutex mtx;
    std::condition_variable cv_not_full;
    std::condition_variable cv_not_empty;
    bool finished = false;
    const size_t max_queue_size = 10; // Maximum number of elements in the queue
    size_t block_size; // Expected size of each vector

    SharedQueue(int blockSize) : block_size(blockSize) {}
};

struct PllConfig {
    float frequency;      // Carrier frequency for the PLL to lock onto
    float Fs;             // Sampling frequency of the input signal
    float ncoScale;       // Scale factor for the NCO (Numerically Controlled Oscillator)
    float phaseAdjust;    // Phase adjustment for the NCO output
    float normBandwidth;  // Normalized bandwidth for the PLL loop filter

    // Constructor with default values matching the function defaults
    PllConfig(float freq, float samplingFreq, float scale = 1.0, float adjust = 0.0, float bandwidth = 0.001)
        : frequency(freq), Fs(samplingFreq), ncoScale(scale), phaseAdjust(adjust), normBandwidth(bandwidth) {}
};

struct PllStateRDS {
    float integrator;    // Integrator state for the loop filter
    float phaseEst;      // Current phase estimate
    float feedbackI;     // In-phase component of the feedback signal
    float feedbackQ;     // Quadrature component of the feedback signal
    float trigOffset;    // Offset for the internal oscillator
    float nextIOut;      // Next in-phase component of the NCO output
    float nextQOut;      // Next quadrature component of the NCO output

    // Constructor for initializing the state
    PllStateRDS() : integrator(0.0), phaseEst(0.0), feedbackI(1.0), feedbackQ(0.0), 
                 trigOffset(0.0), nextIOut(1.0), nextQOut(0.0) {}
};

#endif // DY4_DY4_H
