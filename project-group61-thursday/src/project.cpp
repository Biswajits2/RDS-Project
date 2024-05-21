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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include "mono.h"
#include "logfunc.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>


void popFirstElement(SharedQueue* sharedQueue, std::vector<float>& signal) {
	std::cerr << "Entering Pop First Element " <<  std::endl;
	//std::cerr << "Queue Size Before Pop: " << sharedQueue->queue.size() << std::endl;

    std::unique_lock<std::mutex> lock(sharedQueue->mtx);	
    sharedQueue->cv_not_empty.wait(lock, [sharedQueue] {
        return !sharedQueue->queue.empty() || sharedQueue->finished;
    });

    if (!sharedQueue->queue.empty()) {
        signal = std::move(sharedQueue->queue.front());
        sharedQueue->queue.pop();
		sharedQueue->cv_not_full.notify_one();
    }
	//std::cerr << "Queue Size After Pop: " << sharedQueue->queue.size() << std::endl;
}

void RDS(const std::vector<float>& rds_coeff_64,std::vector<float>& rds_coeff_114, std::vector<float>& rds_state_64, 
		std::vector<float>& rds_state_114, std::vector<float>& rds_state_delay, std::vector<float>& lpf_coeff_rds,
		std::vector<float>& rds_state_lpf, std::vector<float>& rds_rrc_coeff, std::vector<float>& rds_state_rrc, int SPS, bool block_track,
		float rds_ManLastbit, int rds_lastbit, SharedQueue* sharedQueue) {

	while (1) {

		std::cerr << "Entering RDS" << std::endl;
		std::vector<float> rds_data;

		// Demodulate the downsampled I/Q signals
		std::vector<float> fm_demod;
		popFirstElement(sharedQueue, fm_demod);

		//Bandpass Filter from 54kHz to 60kHz
		DSconvo(rds_data, fm_demod, rds_coeff_64, rds_state_64, 1);

		//Delay
		std::vector<float> rds_data_delayed;
		delayBlock(rds_data_delayed, rds_data, rds_state_delay);

		//Squaring Non-Linearity
		std::vector<float> rds_data_squared(rds_data.size());
		for (int i = 0; i < rds_data.size(); i++) {
			rds_data_squared[i] = rds_data[i] * rds_data[i];
		}

		//Bandpass Filter from 113.5kHz to 114.5kHz
		std::vector<float> rds_data_filtered;
		DSconvo(rds_data_filtered, rds_data_squared, rds_coeff_114, rds_state_114, 1);

		//PLLRDS
		PllConfig pllConfig (114e3, 240e3);
		PllStateRDS pllState;
		std::vector<float> ncoOut_I;
		std::vector<float> ncoOut_Q;
		RDS_fmPll(pllState, ncoOut_I, ncoOut_Q, rds_data_filtered, pllConfig);

		/*std::vector<float> dex;
		genIndexVector(dex, (int)ncoOut_I.size());
		logVector("file1", dex, ncoOut_I);

		//put the output of the PLL into the a csv file
		std::ofstream rds_out("rds_out.csv");
		for (int i = 0; i < ncoOut_I.size(); i++) {
			rds_out << ncoOut_I[i] << "," << ncoOut_Q[i] << std::endl;
		}
		rds_out.close();*/

		//Mixer
		std::vector<float> rds_data_mixed(rds_data_filtered.size());
		std::vector<float> rds_data_mixed_Q(rds_data_filtered.size());
		for (int i = 0; i < rds_data_filtered.size(); i++) {
			rds_data_mixed[i] = 2 * ncoOut_I[i] * rds_data_filtered[i];
			rds_data_mixed_Q[i] = 2 * ncoOut_Q[i] * rds_data_filtered[i];
		}

		//Lowpass Filter
		std::vector<float> rds_data_lpf_out;
		DSconvo(rds_data_lpf_out, rds_data_mixed, lpf_coeff_rds, rds_state_lpf, 1);

		//Rastional Resampler
		std::vector<float> rds_data_resampled;
		convoResample(rds_data_resampled, rds_data_lpf_out, lpf_coeff_rds, rds_state_lpf, 1, 1);

		//Root Raised Cosine Filter
		std::vector<float> rds_data_rrc_out;
		DSconvo(rds_data_rrc_out, rds_data_resampled, rds_rrc_coeff, rds_state_rrc, 1);

		// std::ofstream rds_rrc("rds_out_rrc.csv");
		// for (int i = 0; i < rds_data_rrc_out.size(); i++) {
		// 	rds_rrc << rds_data_rrc_out[i] << std::endl;
		// }
		// rds_rrc.close();

		int max_index = 0;
		// If it's the first block, identify the peak value and its index.
		if(block_track) {
			int max = -1; // Initialize max to hold the peak's magnitude (set to -1 initially).
			for (int i = 0; i < SPS; i++) {
				if (std::abs(rds_data_rrc_out[i + 2 * SPS]) > max) {
					max_index = i; // Update max_index to current index if a new peak is found.
					max = std::abs(rds_data_rrc_out[i + 2 * SPS]); // Update max to the new peak's magnitude.
				}
			}
			//block_track = false; // Ensure this block runs only once by setting block_track to false.
		}

		int index = 0; // Initialize index to track the current position in CDR_I and CDR_Q arrays.
		std::vector<float> CDR_I(rds_data_rrc_out.size() / SPS); // Prepare CDR_I to hold in-phase components.
		std::vector<float> CDR_Q(rds_data_rrc_out.size() / SPS); // Prepare CDR_Q to hold quadrature components.

		// Start from max_index, iterate over rds_data_rrc_out with step size SPS to extract CDR samples.
		for (int i = max_index; i < (int)rds_data_rrc_out.size(); i += SPS) {
			CDR_I[index] = rds_data_rrc_out[i]; // Assign the in-phase component to CDR_I.
			index++; // Move to the next position in CDR_I (and CDR_Q).
		}

		// std::ofstream rds_cdr("rds_out_cdr.csv");
		// for (int i = 0; i < CDR_I.size(); i++) {
		// 	rds_cdr << CDR_I[i] << std::endl;
		// }
		// rds_cdr.close();

		//Manchester
		std::vector<int> rds_data_manchester(CDR_I.size());
		for (int i = 0; i < CDR_I.size(); i++) {
			if (i = CDR_I.size() - 1) {
				rds_ManLastbit = CDR_I[i];
			} else if (block_track) {
				if (rds_ManLastbit > 0) {
					rds_data_manchester[i] = 1;
				} else {
					rds_data_manchester[i] = 0;
				}
			} else {
				if (CDR_I[i] > 0) {
					rds_data_manchester[i] = 1;
				} else {
					rds_data_manchester[i] = 0;
				}
			}
		}

		// std::ofstream rds_man("rds_out_man.csv");
		// for (int i = 0; i < rds_data_manchester.size(); i++) {
		// 	rds_man << rds_data_manchester[i] << std::endl;
		// }
		// rds_man.close();

		//Differeintial Decoding
		std::vector<int> rds_data_diff(CDR_I.size());
		rds_data_diff[0] = rds_data_manchester[0];
		for (int i = 1; i < CDR_I.size(); i++) {
			if (!block_track) {
				if (i == (CDR_I.size() - 1)) {
					rds_lastbit = rds_data_manchester[i];
				} else {
					rds_data_diff[i] = rds_data_manchester[i] ^ rds_data_manchester[i - 1];
				}
			}
			else {
				rds_data_diff[i] = rds_data_manchester[0] ^ rds_lastbit;
			}
		}
		block_track = false; // Ensure this block runs only once by setting block_track to false.

		// std::ofstream rds_diff("rds_out_diff.csv");
		// for (int i = 0; i < rds_data_diff.size(); i++) {
		// 	rds_diff << rds_data_diff[i] << std::endl;
		// }
		// rds_diff.close();
	}
}

void stereo(std::vector<float> &audio_block, const std::vector<float> &demod, const std::vector<float> &mono_data, 
			const std::vector<float> &audio_coeff, std::vector<float> &stereo_state, std::vector<float> &carrier_coeffs, 
			std::vector<float> &channel_coeffs, std::vector<float> &carrier_state, std::vector<float> &channel_state, 
			float audio_FS, int audio_expand, int audio_decim, bool needs_resampler)
{
	std::cerr << "Entering Stereo Block" << std::endl;
	std::vector<float> carrier_data;
	auto start_carrier = std::chrono::high_resolution_clock::now();
	DSconvo(carrier_data, demod, carrier_coeffs, carrier_state, 1);
	auto stop_carrier = std::chrono::high_resolution_clock::now();
	auto duration_carrier = std::chrono::duration_cast<std::chrono::microseconds>(stop_carrier - start_carrier);
	std::cerr << "Time taken by Bandpass carrier in Stereo Processing: "
			<< duration_carrier.count() << " microseconds" << std::endl;

	std::vector<float> channel_data;
	auto start_channel = std::chrono::high_resolution_clock::now();
	DSconvo(channel_data, demod, channel_coeffs, channel_state, 1);
	auto stop_channel = std::chrono::high_resolution_clock::now();
	auto duration_channel = std::chrono::duration_cast<std::chrono::microseconds>(stop_channel - start_channel);
	std::cerr << "Time taken by Bandpass carrier in Stereo Processing: "
			<< duration_channel.count() << " microseconds" << std::endl;
	//std::cerr << "Channel Data Size: " << channel_data.size() << std::endl;

	std::vector<float> nco_out;
	PllState pllState;
	auto start_pll = std::chrono::high_resolution_clock::now();
	nco_out = fmPll(carrier_data, 19e3, audio_FS, pllState);
	auto stop_pll = std::chrono::high_resolution_clock::now();
	auto duration_pll = std::chrono::duration_cast<std::chrono::microseconds>(stop_pll - start_pll);
	std::cerr << "Time taken by PLL in Stereo Processing: "
			<< duration_pll.count() << " microseconds" << std::endl;
	//std::cerr << "NCO Data Size: " << nco_out.size() << std::endl;

	std::vector<float> stereoMixed;
	stereoMixed.clear();
	stereoMixed.resize(nco_out.size());
	//Mixer
	auto start_mixer = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < channel_data.size(); i++) {
		stereoMixed[i] = channel_data[i] * nco_out[i] * 2;
	}
	auto stop_mixer = std::chrono::high_resolution_clock::now();
	auto duration_mixer = std::chrono::duration_cast<std::chrono::microseconds>(stop_mixer - start_mixer);
	std::cerr << "Time taken by Mixer in Stereo Processing: "
			<< duration_mixer.count() << " microseconds" << std::endl;

	std::vector<float> processed_stereo;
	if(needs_resampler) {
		auto start_resample_DF = std::chrono::high_resolution_clock::now();
		convoResample(processed_stereo, stereoMixed, audio_coeff, stereo_state, audio_expand, audio_decim);
		auto stop_resample_DF = std::chrono::high_resolution_clock::now();
		auto duration_resample_DF = std::chrono::duration_cast<std::chrono::microseconds>(stop_resample_DF - start_resample_DF);
		std::cerr << "Time taken by convoResample in Stereo Processing for Digital Filtering: "
				<< duration_resample_DF.count() << " microseconds" << std::endl;
	} else {
		auto start_DSconvo_DF = std::chrono::high_resolution_clock::now();
		DSconvo(processed_stereo, stereoMixed, audio_coeff, stereo_state, audio_decim);
		auto stop_DSconvo_DF = std::chrono::high_resolution_clock::now();
		auto duration_DSconvo_DF = std::chrono::duration_cast<std::chrono::microseconds>(stop_DSconvo_DF - start_DSconvo_DF);
		std::cerr << "Time taken by DSconvo in Stereo Processing for Digital Filtering: "
				<< duration_DSconvo_DF.count() << " microseconds" << std::endl;
	}

	//std::cerr << "Processed Stereo Size: " << processed_stereo.size() << std::endl;

	audio_block.clear();
	audio_block.resize(mono_data.size()*2,0);

	//Combiner
	auto start_combiner = std::chrono::high_resolution_clock::now();
	for (unsigned int i = 0; i < mono_data.size(); i++) {
		audio_block[2*i] = (mono_data[i] + processed_stereo[i])/2;
		audio_block[2*i+1] = (mono_data[i] - processed_stereo[i])/2;
	}
	auto stop_combiner = std::chrono::high_resolution_clock::now();
	auto duration_combiner = std::chrono::duration_cast<std::chrono::microseconds>(stop_combiner - start_combiner);
	std::cerr << "Time taken by Combiner in Stereo Processing: "
			<< duration_combiner.count() << " microseconds" << std::endl;
}

// Function to perform mono processing
void monoProcessing(const std::vector<float>& audio_coeff,
                    std::vector<float>& audio_block,
                    int num_taps, int audio_taps,
                    int audio_decim, int audio_expand,
					std::vector<float>& state_audio, bool resample,
					std::vector<float>& channel_state, 
					std::vector<float>& channel_coeffs,
					std::vector<float>& carrier_state,
					std::vector<float>& carrier_coeffs,
					std::vector<float>& delay_state,
					float audio_FS,
					std::vector<float>& stereo_state,
					char path, SharedQueue* sharedQueue) {
	while(1) {
		std::cerr << "Entering monoProcessing" << std::endl;
		//std::cerr << "Audio Coeff Size: " << audio_coeff.size() << std::endl;
		//std::cerr << "Audio Block Size Befroe convo: " << audio_block.size() << std::endl;
		std::vector<float> mono_data;
		std::vector<float> audio_delayed;
		
		// Demodulate the downsampled I/Q signals
		std::vector<float> fm_demod;
		popFirstElement(sharedQueue, fm_demod);
		
		//std::cerr << "FM Demod Size: " << fm_demod.size() << std::endl;
		//std::cerr << "Audio Taps: " << audio_taps << std::endl;

		if (path == 's') {
			std::cerr << "Stereo Processing" << std::endl;
			auto start_delay = std::chrono::high_resolution_clock::now();
			delayBlock(audio_delayed, fm_demod, delay_state);
			auto stop_delay = std::chrono::high_resolution_clock::now();
			auto duration_delay = std::chrono::duration_cast<std::chrono::microseconds>(stop_delay - start_delay);
			std::cerr << "Time taken by delayBlock in Stereo Processing: "
					<< duration_delay.count() << " microseconds" << std::endl;
			//std::cerr << "Delay Block Size " << audio_delayed.size() << std::endl;
			if (resample) {
				std::cerr << "Resampling in Stereo" << std::endl;
				auto start_resample_stereo = std::chrono::high_resolution_clock::now();
				convoResample(mono_data, audio_delayed, audio_coeff, state_audio, audio_expand, audio_decim);
				auto stop_resample_stereo = std::chrono::high_resolution_clock::now();
				auto duration_resample_stereo = std::chrono::duration_cast<std::chrono::microseconds>(stop_resample_stereo - start_resample_stereo);
				std::cerr << "Time taken by convoResample in Stereo Processing: "
						<< duration_resample_stereo.count() << " microseconds" << std::endl;
			} else {
				auto start_DSconvo_stereo = std::chrono::high_resolution_clock::now();
				DSconvo(mono_data, audio_delayed, audio_coeff, state_audio, audio_decim);
				auto stop_DSconvo_stereo = std::chrono::high_resolution_clock::now();
				auto duration_DSconvo_stereo = std::chrono::duration_cast<std::chrono::microseconds>(stop_DSconvo_stereo - start_DSconvo_stereo);
				std::cerr << "Time taken by DSconvo in Stereo Processing: "
						<< duration_DSconvo_stereo.count() << " microseconds" << std::endl;
			}
			//std::cerr << "Mono Data Size after convo: " << mono_data.size() << std::endl;
			stereo(audio_block, fm_demod, mono_data, audio_coeff, stereo_state, carrier_coeffs, channel_coeffs, 
					carrier_state, channel_state, audio_FS, audio_expand, audio_decim, resample);

		} else {
			if (resample) {
				std::cerr << "Resampling" << std::endl;
				auto start_resample_mono = std::chrono::high_resolution_clock::now();
				convoResample(audio_block, fm_demod, audio_coeff, state_audio, audio_expand, audio_decim);
				auto stop_resample_mono = std::chrono::high_resolution_clock::now();
				auto duration_resample_mono = std::chrono::duration_cast<std::chrono::microseconds>(stop_resample_mono - start_resample_mono);
				std::cerr << "Time taken by convoResample in Mono Processing: "
						<< duration_resample_mono.count() << " microseconds" << std::endl;
			} else {
				auto start_DSconvo_mono = std::chrono::high_resolution_clock::now();
				DSconvo(audio_block, fm_demod, audio_coeff, state_audio, audio_decim);
				auto stop_DSconvo_mono = std::chrono::high_resolution_clock::now();
				auto duration_DSconvo_mono = std::chrono::duration_cast<std::chrono::microseconds>(stop_DSconvo_mono - start_DSconvo_mono);
				std::cerr << "Time taken by DSconvo in Mono Processing: "
						<< duration_DSconvo_mono.count() << " microseconds" << std::endl;
			}
		}

		//Getting data ready to write
		//std::cerr << "Converting to short int" << std::endl;
		//std::cerr << "Audio Block Size: " << audio_block.size() << std::endl;
		std::vector<short int> audioBlockInt(audio_block.size());

		for (int i = 0; i < audio_block.size(); i++) {
			if (std::isnan(audio_block[i])) {
				audioBlockInt[i] = 0;
			} else {
				audioBlockInt[i] = static_cast<short int>(audio_block[i] * 16384);
			}
		}

		//std::cerr << "Audio Block Int Size: " << audioBlockInt.size() << std::endl;

		fwrite(&audioBlockInt[0], sizeof(short int), audioBlockInt.size(), stdout);
	}
}

void processBlock(const std::vector<float>& block_data,
                  const std::vector<float>& rf_coeff,
                  int rf_decim,
				  std::vector<float>& state_i, 
				  std::vector<float>& state_q,
				  float prev_I, float prev_Q,  
				  std::vector<float>& fm_demodulated_signal,
				  SharedQueue* sharedQueue) {

	std::cerr << "Entering processBlock" << std::endl;

	//auto sharedQueue = (SharedQueue*)arg;

    // Getting the I's and Q's from the block data
    std::vector<float> iData, qData;
	iData.clear();
	qData.clear();
    iData.reserve(block_data.size() / 2);
    qData.reserve(block_data.size() / 2);

	auto start_split = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < block_data.size(); i++) {
        if (i % 2 > 0) {
            qData.push_back(block_data[i]);
        } else {
            iData.push_back(block_data[i]);
        }
    }
	auto stop_split = std::chrono::high_resolution_clock::now();
	auto duration_split = std::chrono::duration_cast<std::chrono::microseconds>(stop_split - start_split);
	std::cerr << "Time taken by Splitting the data: "
			<< duration_split.count() << " microseconds" << std::endl;

	std::vector<float> i_filter, q_filter;
	

	auto start_convo = std::chrono::high_resolution_clock::now();
    convoDS(i_filter, iData, q_filter, qData, rf_coeff, rf_decim, state_i, state_q);
	auto stop_convo = std::chrono::high_resolution_clock::now();
	auto duration_convo = std::chrono::duration_cast<std::chrono::microseconds>(stop_convo - start_convo);
	std::cerr << "Time taken by convoDS in RF Front end: "
			<< duration_convo.count() << " microseconds" << std::endl;


	fm_demodulated_signal.clear();

	auto start_custom = std::chrono::high_resolution_clock::now();
    std::tie(fm_demodulated_signal, prev_I, prev_Q) = CustomDemod(i_filter, q_filter, prev_I, prev_Q);
	auto stop_custom = std::chrono::high_resolution_clock::now();
	auto duration_custom = std::chrono::duration_cast<std::chrono::microseconds>(stop_custom - start_custom);
	std::cerr << "Time taken by CustomDemod in RF Front end: "
			<< duration_custom.count() << " microseconds" << std::endl;

	//std::cerr << "Queue Size Before Push: " << sharedQueue->queue.size() << std::endl;

	std::unique_lock<std::mutex> lock(sharedQueue->mtx);

	while (sharedQueue->queue.size() >= sharedQueue->max_queue_size) {
            // Wait on the condition variable until notified by the consumer
            sharedQueue->cv_not_full.wait(lock);
	}
	sharedQueue->queue.push(std::move(fm_demodulated_signal));
	//std::cerr << "Queue Size After Push: " << sharedQueue->queue.size() << std::endl;
    lock.unlock();

	sharedQueue->cv_not_empty.notify_one();

	//std::cerr << "Queue Size After Push: " << sharedQueue->queue.size() << std::endl;

}

// Function to process blocks from stdin
void processBlocksFromStdin(int block_size, float rf_FS, int rf_Fc, int rf_decim,
                             int num_taps, SharedQueue* sharedQueue) {

	//Variables and States reuired for IQ processing and audio processing
	std::vector<float> state_i(num_taps - 1, 0.0);
	std::vector<float> state_q(num_taps - 1, 0.0);
	float prev_I = 0.0, prev_Q = 0.0;


	//Finding the coefficients for IQ processing                
	std::vector<float> rf_coeff(num_taps);
	//have a time stap
	auto start_rfimp = std::chrono::high_resolution_clock::now();
	impulseResponseLPF (rf_FS, rf_Fc, num_taps, rf_coeff);
	auto stop_rfimp = std::chrono::high_resolution_clock::now();
	auto duration_rfimp = std::chrono::duration_cast<std::chrono::microseconds>(stop_rfimp - start_rfimp);
	std::cerr << "Time taken by impulseResponseLPF in RF Front end: "
			<< duration_rfimp.count() << " microseconds" << std::endl;


	// Demodulate the downsampled I/Q signals
    std::vector<float> fm_demodulated_signal;

	std::cerr << "Entering processBlocksFromStdin" << std::endl;
    for (unsigned int block_id = 0; ; block_id++) {
        std::vector<float> block_data(block_size);
        readStdinBlockData(block_size, block_id, block_data);
        if ((std::cin.rdstate()) != 0) {
            std::cerr << "End of input stream reached" << std::endl;

			std::lock_guard<std::mutex> lock(sharedQueue->mtx);
        	sharedQueue->finished = true;
			sharedQueue->cv_not_empty.notify_all();

            exit(1);
        }

        processBlock(block_data, rf_coeff, rf_decim, state_i, 
					state_q, prev_I, prev_Q, fm_demodulated_signal, sharedQueue);
    }
}

int main(int argument, char * argument_v[]) {
	//Declaration of variables
	//Integers
	int mode  = 0;
	char path;
	int rf_decim = 0;
	int rf_Fc = 100e3;
	int audio_decim = 0;
	int block_size = 0;
	int num_taps = 101;
	int audio_expand = 1;
	int audio_taps = 101;
	int rf_FS = 0;
	int audio_FS = 0;
	int audio_Fc = 16e3;
	int SPS = 0;
	int RDS_Fs_LPF = SPS * 2375;
	int RDS_fs = 4.56e6;
	int RDS_decim = 0;
	int RDS_expand = 0;

	//Boolean
	bool resample = false;
	bool block_track = true;

	if (argument < 2) {
		//Execute Mode 0 of the argument 1
		std::cerr << "Not enough Arguments" << std::endl;
	} else if (argument == 3) {
		mode = atoi(argument_v[2]);			//converting the argument to integer, this looks at the second argument
		path = argument_v[1][0];			//getting the first character of the first argument
		if (mode < 0 || mode > 3) {
			std::cerr << "Invalid mode: " << mode << std::endl;
			exit(1);
		}
	} else {
		std::cerr << "Invalid entry" << std::endl;
		exit(1);
	}
	//std::cerr << "Argument Vector: " << argument_v[1] << std::endl;
	//std::cerr << "Mode: " << mode << std::endl;
	//std::cerr << "Path: " << path << std::endl;
	//m - mono, s - stereo, r - RDS
	if ((argument_v[1][0] == 'm')) {
		if (mode == 0) {		
			std::cerr << "Executing Mono with mode: 0"<< std::endl;

			//Setting Parameters
			rf_FS = 2.4e6;
			rf_decim = 10;
			audio_FS = 240e3;
			audio_decim = 5;
			block_size = 1024 * rf_decim * audio_decim * 2;

		} else if (mode == 1) {
			std::cerr << "Executing Mono with mode: 1"<< std::endl;

			//Setting Parameters
			rf_FS = 1.152e6;
			rf_decim = 8;
			audio_FS = 144e3;
			audio_decim = 4;
			block_size = 1024 * rf_decim * audio_decim * 2;

		} else if (mode == 2) {
			std::cerr << "Executing Mono with mode: 2"<< std::endl;

			//Setting Parameters
			rf_FS = 2.4e6;
			rf_decim = 10;
			audio_FS = 240e3;
			audio_expand = 147;
			audio_decim = 800;
			audio_taps *= audio_expand;
			block_size = 8 * rf_decim * audio_decim * 2;
			resample = true;

		}else if (mode == 3) {
			std::cerr << "Executing Mono with mode: 3"<< std::endl;
			
			//Setting Parameters
			rf_FS = 1.92e6;
			rf_decim = 6;
			audio_FS = 360e3;
			audio_expand = 441;
			audio_decim = 3200;
			audio_taps *= audio_expand;
			block_size = 4 * rf_decim * audio_decim * 2;
			resample = true;

		}

	} else if (argument_v[1][0] == 's') {
		//std::cerr << "Mode: " << mode << std::endl;
		if (mode == 0) {
			std::cerr << "Executing Stereo with mode: 0"<< std::endl;

			//Setting Parameters
			rf_FS = 2.4e6;
			rf_decim = 10;
			audio_FS = 240e3;
			audio_decim = 5;
			block_size = 1024 * rf_decim * audio_decim * 2;
		} else if (mode == 1) {
			std::cerr << "Executing Stereo with mode: 1"<< std::endl;

			//Setting Parameters
			rf_FS = 1.152e6;
			rf_decim = 8;
			audio_FS = 144e3;
			audio_decim = 4;
			block_size = 1024 * rf_decim * audio_decim * 2;
		} else if (mode == 2) {
			std::cerr << "Executing Mono with mode: 2"<< std::endl;

			//Setting Parameters
			rf_FS = 2.4e6;
			rf_decim = 10;
			audio_FS = 240e3;
			audio_expand = 147;
			audio_decim = 800;
			audio_taps = 5000;
			block_size = 8 * rf_decim * audio_decim * 2;
			resample = true;
			
		} else if (mode == 3) {
			std::cerr << "Executing Mono with mode: 3"<< std::endl;

			//Setting Parameters
			rf_FS = 1.92e6;
			rf_decim = 6;
			audio_FS = 360e3;
			audio_expand = 441;
			audio_decim = 3200;
			audio_taps = 12800;
			block_size = 4 * rf_decim * audio_decim * 2;
			resample = true;
		}

	} else if (argument_v[1][0] == 'r') {
		if (mode == 0) {
			std::cerr << "Executing RDS with mode: 0"<< std::endl;

			//Setting Parameters
			rf_FS = 2.4e6;
			rf_decim = 10;
			audio_FS = 240e3;
			RDS_fs = 186.96e6;
			audio_decim = 5;
			SPS = 41;
			RDS_Fs_LPF = SPS * 2375;
			RDS_decim = 1920;
			RDS_expand = 779;
			block_size = 1024 * rf_decim * audio_decim * 2;
		} else if (mode == 2) {
			std::cerr << "Executing RDS with mode: 2"<< std::endl;

			//Setting Parameters
			rf_FS = 2.4e6;
			rf_decim = 10;
			audio_FS = 240e3;
			audio_expand = 147;
			audio_decim = 800;
			audio_taps = 5000;
			SPS = 12;
			RDS_Fs_LPF = SPS * 2375;
			RDS_fs = 4.56e6;
			RDS_decim = 160;
			RDS_expand = 19;
			block_size = 8 * rf_decim * audio_decim * 2;
			resample = true;
		}
	} else {
		std::cerr << (int)argument_v[1] << std::endl;
		std::cerr << "Invalid Process - Mode Block" << std::endl;
		exit(1);
	}

	//vectors for the coefficients used for Mono processing and Stereo processing
	std::vector<float> audio_coeff(audio_taps);

	//vectors for the coefficients used for Stereo processing
	std::vector<float> channel_coeffs;
	std::vector<float> carrier_coeffs;

	//Vectors for the coefficients used for RDS processing
	std::vector<float> rds_coeff_64;
	std::vector<float> rds_coeff_114;
	std::vector<float> lpf_coeff_rds;
	std::vector<float> rds_rrc_coeff;

	//State for Mono processing
	std::vector<float> state_audio(audio_taps - 1, 0.0);

	//States for Stereo processing
	std::vector<float> channel_state(audio_taps - 1, 0.0);
	std::vector<float> carrier_state(audio_taps - 1, 0.0);
	std::vector<float> stereo_state(audio_taps - 1, 0.0);
	std::vector<float> delay_state((num_taps - 1) / 2, 0.0);

	//States for RDS processing
	std::vector<float> rds_state_64(audio_taps - 1, 0.0);
	std::vector<float> rds_state_114(audio_taps - 1, 0.0);
	std::vector<float> rds_state_delay((num_taps - 1) / 2, 0.0);
	std::vector<float> rds_state_lpf((num_taps - 1) / 2, 0.0);
	std::vector<float> rds_state_rrc((num_taps - 1) / 2, 0.0);
	int rds_lastbit = 0;
	float rds_ManLastbit = 0;

	//Audio Block for Mono Processing
	std::vector<float> audio_block;

	//Finding the coefficients for Mono processing
	impulseResponseLPF (audio_FS * audio_expand, audio_Fc, audio_taps, audio_coeff);

	//Finding the coefficients for Stereo processing
	bandPassFilter (22e3, 54e3, audio_FS, num_taps, channel_coeffs);
	bandPassFilter (18.5e3 , 19.5e3, audio_FS, num_taps, carrier_coeffs);

	//Finding the coefficients for RDS processing
	bandPassFilter (54e3, 60e3, audio_FS, num_taps, rds_coeff_64);
	bandPassFilter(113.5e3, 114.5e3, audio_FS, num_taps, rds_coeff_114);
	impulseResponseLPF (RDS_Fs_LPF, 3e3, num_taps, lpf_coeff_rds);
	impulseResponseRootRaisedCosine(RDS_Fs_LPF, num_taps, rds_rrc_coeff);


	//Creating the Shared Queue for Multi-threading
	SharedQueue* sharedQueue = new SharedQueue{(block_size/2/rf_decim)};

	//Creating and calling the thread for RF Front End Processing
	std::thread processBlocksThread(processBlocksFromStdin, block_size, rf_FS, rf_Fc, rf_decim, num_taps, sharedQueue);
	
	//Creating and calling the threads for Mono and Stereo Processing
	std::thread monoProcessingThread(monoProcessing, std::ref(audio_coeff), std::ref(audio_block), num_taps, audio_taps, 
								audio_decim, audio_expand, std::ref(state_audio), resample, std::ref(channel_state), 
								std::ref(channel_coeffs), std::ref(carrier_state), std::ref(carrier_coeffs), 
								std::ref(delay_state), audio_FS, std::ref(stereo_state), path,sharedQueue);

	/*std::thread rdsThread(RDS, std::ref(rds_coeff_64), std::ref(rds_coeff_114), std::ref(rds_state_64),
						std::ref(rds_state_114), std::ref(rds_state_delay), std::ref(lpf_coeff_rds),
						std::ref(rds_state_lpf), std::ref(rds_rrc_coeff), std::ref(rds_state_rrc), SPS, block_track,
						rds_ManLastbit, rds_lastbit,sharedQueue);*/
				

	processBlocksThread.join();
	monoProcessingThread.join();
	//rdsThread.join();
	
	/*processBlocksFromStdin(block_size, rf_FS, rf_Fc, rf_decim, 
							num_taps, sharedQueue);

	// Mono Processing
    monoProcessing(audio_coeff, audio_block, num_taps, audio_taps, 
					audio_decim, audio_expand, state_audio, resample, channel_state, channel_coeffs,
					carrier_state, carrier_coeffs, delay_state, audio_FS, stereo_state, path, sharedQueue);*/
	return 1;
}

