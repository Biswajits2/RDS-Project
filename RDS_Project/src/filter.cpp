/*
Comp Eng 3DY4 (Computer Systems unsigned integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
    h.clear();
	h.resize(num_taps, 0.0);
    float Norm_cutoff = Fc / (Fs / 2.0);
    for (unsigned int i = 0; i < num_taps; i++) {
        if (i == (num_taps - 1) / 2) {
            h[i] = Norm_cutoff;
        } else {
            float numerator = Norm_cutoff * sin(PI * Norm_cutoff * (i - (num_taps - 1) / 2.0));
            float denominator = Norm_cutoff * PI * (i - (num_taps - 1) / 2.0);
            h[i] = numerator / denominator;
        }
        h[i] *= pow(sin((i * PI) / num_taps), 2);  // Window function
    }
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
    y.clear();
    y.resize(x.size() + h.size() - 1, 0.0);

    // Perform the convolution
    for (size_t n = 0; n < y.size(); ++n) {
        for (size_t k = 0; k < h.size(); ++k) {
            if (n >= k && n - k < x.size()) {
                y[n] += h[k] * x[n - k];
            }
        }
    }
}

std::vector<float> convolveFIR_State(const std::vector<float> &filter, const std::vector<float> &inputBlock, std::vector<float> &state) {
    std::vector<float> outputBlock(inputBlock.size(), 0.0);
    outputBlock.clear();

    for (unsigned int outputIndex = 0; outputIndex < (unsigned int) outputBlock.size(); outputIndex++) {
        for (unsigned int filterIndex = 0; filterIndex < (unsigned int) filter.size(); filterIndex++) {
            unsigned int inputIndex = outputIndex - filterIndex;
			if (inputIndex >= 0) {
                outputBlock[outputIndex] += filter[filterIndex] * inputBlock[inputIndex];
            }
            else {
                outputBlock[outputIndex] += filter[filterIndex] * state[state.size() + inputIndex];
            }
        }
    }

    for (unsigned int stateIndex = 0; stateIndex < filter.size(); stateIndex++) {
		 state[stateIndex] = inputBlock[inputBlock.size() - state.size() + stateIndex];
    }
    return outputBlock;
}

void processAudioBlocks(const std::vector<float>& audio_left, const std::vector<float>& audio_right,
                         std::vector<float>& filteredLeft, std::vector<float>& filteredRight,
                         unsigned int blockSize, const std::vector<float>& h,
                         std::vector<float>& stateLeft, std::vector<float>& stateRight) {
    unsigned int position = 0;

    while (position < audio_left.size()) {
        // Extract a block of audio samples from the left and right channels
        std::vector<float> blockLeft(audio_left.begin() + position, audio_left.begin() + position + blockSize);
        std::vector<float> blockRight(audio_right.begin() + position, audio_right.begin() + position + blockSize);

        // Process the left and right blocks using a function 'convolveFIR_State' with parameters 'h', 'blockLeft', and 'stateLeft'
        std::vector<float> sumLeft = convolveFIR_State(h, blockLeft, stateLeft);
        std::vector<float> sumRight = convolveFIR_State(h, blockRight, stateRight);

        // Append the processed samples to the respective filteredLeft and filteredRight vectors
        filteredLeft.insert(filteredLeft.end(), sumLeft.begin(), sumLeft.end());
        filteredRight.insert(filteredRight.end(), sumRight.begin(), sumRight.end());

        // Move the position pounsigned inter to the next block
        position += blockSize;
    }
}

//Block Processing from Lab 3 refactored unsigned into CPP
std::pair<std::vector<float>, std::vector<float>> convo(const std::vector<float>& input_sig, const std::vector<float>& filter_coeff, std::vector<float> initial_state) {
    unsigned int input_len = input_sig.size();
    unsigned int filter_len = filter_coeff.size();
    unsigned int output_len = input_len;
    std::vector<float> output_sig(output_len, 0.0);

   //std::cout << "Size of the input signal is " << input_len << std::endl;
    //std::cout << "Size of the output is " << output_len() << std::endl;

    // Flip filter coefficient for convolution
    std::vector<float> flipped_filter_coeff = filter_coeff;
    std::reverse(flipped_filter_coeff.begin(), flipped_filter_coeff.end());

    // Convolution
    for (unsigned int i = 0; i < output_len; i++) {
        for (unsigned int j = 0; j < filter_len; j++) {
            if (i - j >= 0) {
                output_sig[i] += input_sig[i - j] * flipped_filter_coeff[j];
            } else {
                // Handling indices outside the input signal range using initial_state
                output_sig[i] += initial_state[initial_state.size() + (i - j)] * flipped_filter_coeff[j];
            }
        }
        //printf("Output Signal: %f\n", output_sig[i]);
    }
    //print the size of the output signal
    //std::cout << "Size of the output signal is " << output_sig.size() << std::endl;

    // Update initial_state based on output_len and input_len
    if (output_len < initial_state.size()) {
        std::copy(initial_state.begin() + input_len, initial_state.end(), initial_state.begin());
        std::copy(input_sig.begin(), input_sig.end(), initial_state.end() - input_len);
    } else {
        initial_state.assign(input_sig.end() - initial_state.size(), input_sig.end());
    }

    return {output_sig, initial_state};
}

void convoDS(std::vector<float> &iOut, const std::vector<float> &iData,
                           std::vector<float> &qOut, const std::vector<float> &qData,
                           const std::vector<float> &filter_coeff, const unsigned int rf_decim,
                           std::vector<float> &iState, std::vector<float> &qState) {
    
    //std::cerr << "I Data Size: " << iData.size() << std::endl;
	//std::cerr << "Q Data Size: " << qData.size() << std::endl;
    //std::cerr << "RF Decim: " << rf_decim << std::endl;
    // Clear and resize the output vectors based on the downsampling factor
    iOut.clear();
    iOut.resize(iData.size() / rf_decim, 0.0);
    qOut.clear();
    qOut.resize(qData.size() / rf_decim, 0.0);

    //std::cerr << "Size of iOut: " << iOut.size() << std::endl;
    //std::cerr << "Size of qOut: " << qOut.size() << std::endl;

    // Convolution loop with downsampling
    for (int i = 0; i < (int)iOut.size(); i++) {
        for (int j = 0; j < (int)filter_coeff.size(); j++) {
            int inputIndex = rf_decim * i - j;
            //std::cerr << "Input Index: " << inputIndex << std::endl;
            if (inputIndex >= 0) {
                iOut[i] += filter_coeff[j] * iData[inputIndex];
                qOut[i] += filter_coeff[j] * qData[inputIndex];
            } else {
                // Use state for initial samples
                iOut[i] += filter_coeff[j] * iState[iState.size() - abs(inputIndex)];
                qOut[i] += filter_coeff[j] * qState[qState.size() - abs(inputIndex)];
            }
        }
    }


    // Update state with the last samples of the current block
    std::vector<float> tempiState(&iData[iData.size() - iState.size()], &iData[iData.size()]);
    iState = tempiState;
    std::vector<float> tempqState(&qData[qData.size() - qState.size()], &qData[qData.size()]);
    qState = tempqState;
}

void DSconvo(std::vector<float> &output_sig, const std::vector<float> &input_sig, const std::vector<float> &filter_coeff, std::vector<float> &initial_state, int audio_decim)
{
    output_sig.clear();
    output_sig.resize(input_sig.size() / audio_decim, 0.0);
    //std::cout << "Size of the output signal is " << output_sig.size() << std::endl;
    //std::cout << "Size of the input signal is " << input_sig.size() << std::endl;
    //std::cout << "Size of the filter coefficient is " << filter_coeff.size() << std::endl;
    //std::cout << "Size of the Intial State is " << initial_state.size() << std::endl;
    for(int i = 0; i < (int)output_sig.size(); i++) {
        for(int j = 0; j < (int)filter_coeff.size(); j++) {
            if(audio_decim * i - j >= 0) {
                output_sig[i] += filter_coeff[j] * input_sig[audio_decim * i - j];
            }
            else {
                output_sig[i] += filter_coeff[j] * initial_state[initial_state.size() + (audio_decim * i - j)];
            }
        }
    }

    for(int n = 0; n < (int)initial_state.size(); n++) {
        initial_state[n] = input_sig[input_sig.size() - initial_state.size() + n];
    }
}

void convoResample(std::vector<float> &output_sig, const std::vector<float> &input_sig, const std::vector<float> &filter_coeff, std::vector<float> &initial_state, int audio_up, int audio_decim) {
    output_sig.clear(); 
    output_sig.resize(input_sig.size()*audio_up/audio_decim, 0.0);
    for (int i = 0; i < (int)output_sig.size(); i++) {
        for (int j = (i * audio_decim) % audio_up; j < (int)filter_coeff.size(); j += audio_up) {
            int n = ((i * audio_decim) - j) / audio_up;
            if(n >= 0) {
                output_sig[i] += filter_coeff[j]*input_sig[n];
            }
            else {
                output_sig[i] += filter_coeff[j] * initial_state[initial_state.size() + n];
            }
        }
        output_sig[i] *= audio_up;
    }
    std::vector<float> temp(&input_sig[input_sig.size() - initial_state.size()], &input_sig[input_sig.size()]);
    initial_state = temp;
}

//bandpass filter (referenced using LPF filter)

void bandPassFilter(float FB, float FE, float Fs, unsigned short int num_taps, std::vector<float> &h)   {
  h.clear(); 
  h.resize(num_taps, 0.0);

  float Nc = ((FE + FB)/2)/(Fs/2);
  float Np = (FE - FB)/(Fs/2);

    for (int i = 0; i < num_taps - 1; i++) {
        if (i == (num_taps - 1) / 2) {
            h[i] = Np;
        }
        else {
            h[i] = Np * sin(Np/2 * PI * (i - (num_taps - 1) / 2)) / (Np/2 * PI * (i - (num_taps - 1) / 2));
        }
        h[i] = h[i] * (cos(i*PI*Nc) * pow(sin(i * PI/num_taps), 2));
        //h[i] = h[i]* pow(sin(i * PI/num_taps), 2);
    }
}

void fmPLL(std::vector<float> pllIn, const float freq, const float Fs, 
        const float phaseAdjust, const float normBandwidth, 
        float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, 
        std::vector<float> &ncoOut, float &trigOffset, float &state) {

    float ncoScale = 2.0;

    float Cp = 2.666;
    float Ci = 3.555;

    float errorI;
    float errorQ;
    float errorD;

    float Kp = normBandwidth*Cp;
    float Ki = (normBandwidth*normBandwidth)*Ci;

    ncoOut.clear(); 
    ncoOut.resize((int)pllIn.size(), 0.0);

    float trigArg;
    ncoOut[0] = state;
    for(int i = 0; i < (int)pllIn.size(); i++){
        //phase detector
        errorI = pllIn[i] * (feedbackI);
        errorQ = pllIn[i] * (-1 * feedbackQ);
        //four-quadrant arctangent discriminator for phase error detection
        errorD = std::atan2(errorQ,errorI);
        //loop filter
        integrator = integrator + Ki*errorD;
        //update phase estimate
        phaseEst = phaseEst + Kp*errorD + integrator;
        //internal oscillator
        trigOffset += 1;
        trigArg = 2*PI*(freq/Fs)*(trigOffset) + phaseEst;
        feedbackI = std::cos(trigArg);
        feedbackQ = std::sin(trigArg);
        ncoOut[i+1] = std::cos(trigArg * ncoScale + phaseAdjust);
    }
    state = ncoOut[ncoOut.size()-1];
}

std::vector<float> fmPll(const std::vector<float>& pllIn, float freq, float Fs,
                         PllState& state, float ncoScale,
                         float phaseAdjust, float normBandwidth) {
    const float Cp = 2.666f;
    const float Ci = 3.555f;

    float Kp = normBandwidth * Cp;
    float Ki = (normBandwidth * normBandwidth) * Ci;

    std::vector<float> ncoOut(pllIn.size() + 1, 0);

    for (size_t k = 0; k < pllIn.size(); ++k) {
        float errorI = pllIn[k] * state.feedbackI;
        float errorQ = pllIn[k] * (-state.feedbackQ);

        float errorD = std::atan2(errorQ, errorI);

        state.integrator += Ki * errorD;
        state.phaseEst += Kp * errorD + state.integrator;

        state.trigOffset += 1;
        float trigArg = 2 * PI * (freq / Fs) * state.trigOffset + state.phaseEst;
        state.feedbackI = std::cos(trigArg);
        state.feedbackQ = std::sin(trigArg);

        ncoOut[k + 1] = std::cos(trigArg * ncoScale + phaseAdjust);
    }

    return ncoOut;
}

void delayBlock(std::vector<float> &output, const std::vector<float> &input, std::vector<float> &state) {
    output.clear();
    // Resize the output vector to match the size of the input vector.
    output.resize(input.size(), 0.0);

    // Loop through the output vector.
    for(int i = 0; i < (int)output.size(); i++) {
        // For the first part, use the state vector's contents.
        if(i < (int)state.size()) {
            output[i] = state[i];
        }
        // For the remaining part, copy from the input vector.
        else {
            output[i] = input[i - state.size()];
        }
    }

    // Update the state vector with the last 'state.size()' elements of the input vector.
    for(int j = 0; j < (int)state.size(); j++) {
        state[j] = input[input.size() - state.size() + j];
    }
}

void RDS_fmPll(PllStateRDS &state, std::vector<float> &ncoOut_I, std::vector<float> &ncoOut_Q, 
                const std::vector<float> &pllIn, const PllConfig &config) {
    const float Cp = 2.666;
    const float Ci = 3.555;

    float Kp = config.normBandwidth * Cp;
    float Ki = (config.normBandwidth * config.normBandwidth) * Ci;

    ncoOut_I.clear(); 
    ncoOut_I.resize(pllIn.size(), 0.0);

    ncoOut_Q.clear(); 
    ncoOut_Q.resize(pllIn.size(), 0.0);
    
    //std::cerr << "Size of pllIn: " << pllIn.size() << std::endl;
    //std::cerr << "Size of ncoOut_I: " << ncoOut_I.size() << std::endl;
    //std::cerr << "Size of ncoOut_Q: " << ncoOut_Q.size() << std::endl;


    float errorI, errorQ, errorD, trigArg;

    for(unsigned int k = 0; k < pllIn.size(); k++) {
        errorI = pllIn[k] * state.feedbackI;  // In-phase error
        errorQ = pllIn[k] * -state.feedbackQ; // Quadrature-phase error

        // Four-quadrant arctangent discriminator for phase error detection
        errorD = std::atan2(errorQ, errorI);

        // Loop filter
        state.integrator += Ki * errorD;
        // Update phase estimate
        state.phaseEst += Kp * errorD + state.integrator;

        // Internal oscillator
        state.trigOffset += 1;
        trigArg = 2 * M_PI * (config.frequency / config.Fs) * state.trigOffset + state.phaseEst;
        state.feedbackI = std::cos(trigArg);
        state.feedbackQ = std::sin(trigArg);

        ncoOut_I[k] = std::cos(trigArg * config.ncoScale + config.phaseAdjust);
        ncoOut_Q[k] = std::sin(trigArg * config.ncoScale + config.phaseAdjust);
    }

    // Save the last NCO output state for block continuity
    state.nextIOut = ncoOut_I.back();
    state.nextQOut = ncoOut_Q.back();
}

//RRC conversion from python to C++    [RDS]
void impulseResponseRootRaisedCosine(float Fs, int N_taps, std::vector<float> &impulseResponseRRC) {
    float T_symbol = 1/2375.0;
    float beta = 0.90;

    impulseResponseRRC.clear();
    impulseResponseRRC.resize(N_taps, 0.0);

    for(int k=0; k < N_taps; k++){
        float t  = float(k-N_taps/2)/Fs;
        if( t == 0.0){
            impulseResponseRRC[k] = 1.0 + beta*((4/PI)-1);
        }
        else if(t == (-1 * T_symbol)/(4*beta) || t == (1*T_symbol)/(4*beta)){
            impulseResponseRRC[k] = (beta/sqrt(2))*(((1 + 2/PI)* (sin(PI/(4*beta)))) + ((1-2/PI)*(cos(PI/(4*beta)))));
        }
        else{
            impulseResponseRRC[k] = (sin(PI*t*(1-beta)/T_symbol) + 4*beta*(t/T_symbol)*cos(PI*t*(1+beta)/T_symbol))/ (PI*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol);
        }
    }
}
 