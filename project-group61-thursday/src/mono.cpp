#include "dy4.h"
#include "filter.h"
#include "mono.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <iofunc.h>

unsigned int BLOCK_SIZE;
//Assume we are give the data when calling the function
//i.e. the front end will be responsible for reading the data from the file

//downsampling using lecture concepts
void downSampling(std::vector<float> &input, std::vector<float> &output, unsigned int decimate_factor) {
  output.clear();
  output.resize(input.size()/decimate_factor);    // the size might be different

  for(unsigned int i = 0; i<(unsigned int)output.size(); i++) {
    output[i] = input[i*decimate_factor];
  }
}

void downsample(const std::vector<float>& input, std::vector<float>& output, unsigned int decim_factor) {
    output.clear(); // Ensure the output vector is empty before starting
    for (size_t i = 0; i < input.size(); i += decim_factor) {
        output.push_back(input[i]);
    }
}

void expand(const std::vector<float>& input, std::vector<float>& output, unsigned int expand_factor) {
    output.clear(); // Ensure the output vector is empty before starting
    for (size_t i = 0; i < input.size(); ++i) {
        for (unsigned int j = 0; j < expand_factor; ++j) {
            output.push_back(input[i]); // Insert each sample expand_factor times
        }
    }
}
std::tuple<std::vector<float>, float, float> CustomDemod(const std::vector<float>& I, const std::vector<float>& Q, float prev_I = 0.0, float prev_Q = 0.0) {
    std::vector<float> fm_demod(I.size());
    
    for (size_t k = 0; k < (int)I.size(); k++) {
        float denominator = (I[k]*I[k]) + (Q[k]*Q[k]);
        
        if (I[k] == 0.0 || Q[k] == 0.0) {
            fm_demod[k] = 0.0;
        } else {
            fm_demod[k] = ((I[k] * (Q[k] - prev_Q)) - (Q[k] * (I[k] - prev_I))) / denominator;
        }
        
        prev_I = I[k];
        prev_Q = Q[k];
    }

    return {fm_demod, prev_I, prev_Q};
}