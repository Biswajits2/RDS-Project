
#ifndef DY4_MONO_H
#define DY4_MONO_H

// add headers as needed
#include <iostream>
#include <vector>
#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include "iofunc.h"

// declaration of a function prototypes
void downSampling(std::vector<float> &input, std::vector<float> &output, int decimate_factor);
void downsample(const std::vector<float>& input, std::vector<float>& output, unsigned int decim_factor);
void expand(const std::vector<float>& input, std::vector<float>& output, unsigned int expand_factor);
std::tuple<std::vector<float>, float, float> CustomDemod(const std::vector<float>& I, const std::vector<float>& Q, float prev_I, float prev_Q);

#endif // DY4_MONO_H
