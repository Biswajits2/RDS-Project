/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
std::pair<std::vector<float>, std::vector<float>> convo(const std::vector<float>& input_sig, const std::vector<float>& filter_coeff, std::vector<float> initial_state);
void convoDS(std::vector<float> &iOut, const std::vector<float> &iData,
                           std::vector<float> &qOut, const std::vector<float> &qData,
                           const std::vector<float> &filter_coeff, const int unsigned rf_decim,
                           std::vector<float> &iState, std::vector<float> &qState);
void DSconvo(std::vector<float> &output_sig, const std::vector<float> &input_sig, const std::vector<float> &filter_coeff, std::vector<float> &initial_state, int audio_decim);
void convoResample(std::vector<float> &output_sig, const std::vector<float> &input_sig, const std::vector<float> &filter_coeff, std::vector<float> &initial_state, int audio_up, int audio_decim);
void bandPassFilter(float FB, float FE, float Fs, unsigned short int num_taps, std::vector<float> &h);
void fmPLL(std::vector<float> pllIn, const float freq, const float Fs, 
        const float phaseAdjust, const float normBandwidth, float &integrator, 
        float &phaseEst, float &feedbackI, float &feedbackQ, 
        std::vector<float> &ncoOut, float &trigOffset, float &state);
std::vector<float> fmPll(const std::vector<float>& pllIn, float freq, float Fs,
                         PllState& state, float ncoScale = 2.0f,
                         float phaseAdjust = 0.0f, float normBandwidth = 0.01f);
void delayBlock(std::vector<float> &output, const std::vector<float> &input, std::vector<float> &state);
void RDS_fmPll(PllStateRDS &state, std::vector<float> &ncoOut_I, std::vector<float> &ncoOut_Q, 
               const std::vector<float> &pllIn, const PllConfig &config);
void impulseResponseRootRaisedCosine(float Fs, int N_taps, std::vector<float> &impulseResponseRRC);
#endif // DY4_FILTER_H
