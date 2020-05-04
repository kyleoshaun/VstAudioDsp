//------------------------------------------------------------------------
// Project     : VST SDK
//
// Category    : Examples
// Filename    : plugprocessor.cpp
// Created by  : Steinberg, 01/2018
// Description : HelloWorld Example for VST 3
//
//-----------------------------------------------------------------------------
// LICENSE
// (c) 2019, Steinberg Media Technologies GmbH, All Rights Reserved
//-----------------------------------------------------------------------------
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// 
//   * Redistributions of source code must retain the above copyright notice, 
//     this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation 
//     and/or other materials provided with the distribution.
//   * Neither the name of the Steinberg Media Technologies nor the names of its
//     contributors may be used to endorse or promote products derived from this 
//     software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
// IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.
//-----------------------------------------------------------------------------

#include "../include/plugprocessor.h"
#include "../include/plugids.h"
#include "../include/firFilterDemoCoeffs.h"

#include <cmath>
#include <math.h>

#include "base/source/fstreamer.h"
#include "pluginterfaces/base/ibstream.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"
#include "public.sdk/source/vst/vstaudioprocessoralgo.h"

namespace Steinberg {
namespace FilterDemo {

//-----------------------------------------------------------------------------
FilterDemoProcessor::FilterDemoProcessor ()
{
	// Init members
	outputDelayBuf = 0;
	mBypass = false;
	fbFilterCoeffsA = 0;
	resonatorPoleRadius = 0;
	resonatorFreq = 0;
	ffwFilterCoeffsB = 0;
	inputDelayBuf = 0;
	noOfZeroPairs = 1;
	noOfPolePairs = 0;
	noOfSecOrderFilters = std::max(noOfZeroPairs, noOfPolePairs);
	cutoffFreq = 0;
	zeroRadius = 1;
	// register its editor class
	setControllerClass (MyControllerUID);
}

//-----------------------------------------------------------------------------
tresult PLUGIN_API FilterDemoProcessor::initialize (FUnknown* context)
{
	//---always initialize the parent-------
	tresult result = AudioEffect::initialize (context);
	if (result != kResultTrue)
		return kResultFalse;

	//---create Audio In/Out buses------
	// we want a stereo Input and a Stereo Output
	addAudioInput (STR16 ("AudioInput"), Vst::SpeakerArr::kStereo);
	addAudioOutput (STR16 ("AudioOutput"), Vst::SpeakerArr::kStereo);

	return kResultTrue;
}

//-----------------------------------------------------------------------------
tresult PLUGIN_API FilterDemoProcessor::setBusArrangements (Vst::SpeakerArrangement* inputs,
                                                            int32 numIns,
                                                            Vst::SpeakerArrangement* outputs,
                                                            int32 numOuts)
{
	// we only support one in and output bus and these buses must have the same number of channels
	if (numIns == 1 && numOuts == 1 && inputs[0] == outputs[0])
	{
		return AudioEffect::setBusArrangements (inputs, numIns, outputs, numOuts);
	}
	return kResultFalse;
}

//-----------------------------------------------------------------------------
tresult PLUGIN_API FilterDemoProcessor::setupProcessing (Vst::ProcessSetup& setup)
{
	// here you get, with setup, information about:
	// sampleRate, processMode, maximum number of samples per audio block
	return AudioEffect::setupProcessing (setup);
}

//-----------------------------------------------------------------------------
tresult PLUGIN_API FilterDemoProcessor::setActive (TBool state)
{
	Vst::SpeakerArrangement arr;
	if (getBusArrangement(Vst::kOutput, 0, arr) != kResultTrue)
		return kResultFalse;
	int32 numChannels = Vst::SpeakerArr::getChannelCount(arr);
	size_t outputDelayBufSize = sizeof(Vst::Sample64) * BIQUAD_NO_OF_A_COEFFS + 0.5; //max 1 sec delay (round up)
	size_t inputDelayBufSize  = sizeof(Vst::Sample64) * BIQUAD_NO_OF_B_COEFFS + 0.5;
	size_t secOrderFilterCoeffSize = 2 * sizeof(double);

	if (state) // Initialize
	{
		// Delay Buffers
		outputDelayBuf = (double***)std::malloc(sizeof(double**) * numChannels);
		inputDelayBuf  = (double***)std::malloc(sizeof(double**) * numChannels);
		for (int channelIdx = 0; channelIdx < numChannels; channelIdx++)
		{
			inputDelayBuf[channelIdx]  = (double**)std::malloc(sizeof(double*) * MAX_NO_OF_2ND_ORDER_FILTERS);
			outputDelayBuf[channelIdx] = (double**)std::malloc(sizeof(double*) * MAX_NO_OF_2ND_ORDER_FILTERS);
			double** channelInputDelayBuf  = inputDelayBuf[channelIdx];
			double** channelOutputDelayBuf = outputDelayBuf[channelIdx];
			for (int filterIdx = 0; filterIdx < MAX_NO_OF_2ND_ORDER_FILTERS; filterIdx++)
			{
				channelInputDelayBuf[filterIdx]  = (double*)std::malloc(inputDelayBufSize);
				channelOutputDelayBuf[filterIdx] = (double*)std::malloc(outputDelayBufSize);

				memset(channelInputDelayBuf[filterIdx],  0, inputDelayBufSize);
				memset(channelOutputDelayBuf[filterIdx], 0, outputDelayBufSize);
			}
		}


		// Filter Coefficients
		ffwFilterCoeffsB       = (double**)std::malloc(sizeof(double*) * MAX_NO_OF_2ND_ORDER_FILTERS);
		fbFilterCoeffsA = (double**)std::malloc(sizeof(double*) * MAX_NO_OF_2ND_ORDER_FILTERS);
		for (int filterIdx = 0; filterIdx < MAX_NO_OF_2ND_ORDER_FILTERS; filterIdx++)
		{
			ffwFilterCoeffsB[filterIdx] = (double*)std::malloc(secOrderFilterCoeffSize);
			fbFilterCoeffsA[filterIdx]  = (double*)std::malloc(secOrderFilterCoeffSize);

			memset(ffwFilterCoeffsB[filterIdx], 0, secOrderFilterCoeffSize);
			memset(fbFilterCoeffsA[filterIdx], 0, secOrderFilterCoeffSize);

			//Initialize as bypass (Set b1 = 1 for each 2nd order component)
			double* filterFfwCoeffsA = ffwFilterCoeffsB[filterIdx];
			filterFfwCoeffsA[0] = 1;
		}
	}
	else // Release
	{
		if (outputDelayBuf)
		{
			//Delay Buffers
			for (int channelIdx = 0; channelIdx < numChannels; channelIdx++)
			{
				double** channelInputDelayBuf = inputDelayBuf[channelIdx];
				double** channelOutputDelayBuf = outputDelayBuf[channelIdx];
				for (int filterIdx = 0; filterIdx < MAX_NO_OF_2ND_ORDER_FILTERS; filterIdx++)
				{
					std::free(channelInputDelayBuf[filterIdx]);
					std::free(channelOutputDelayBuf[filterIdx]);
				}
				std::free(outputDelayBuf[channelIdx]);
				std::free(inputDelayBuf[channelIdx]);
			}
			std::free(outputDelayBuf);
			outputDelayBuf = 0;
			std::free(inputDelayBuf);
			inputDelayBuf = 0;

			//Filter Coefficients
			for (int filterIdx = 0; filterIdx < MAX_NO_OF_2ND_ORDER_FILTERS; filterIdx++)
			{
				std::free(ffwFilterCoeffsB[filterIdx]);
				std::free(fbFilterCoeffsA[filterIdx]);
			}
			std::free(fbFilterCoeffsA);
			fbFilterCoeffsA = 0;
		}
	}
	return AudioEffect::setActive (state);
}

//-----------------------------------------------------------------------------
template <typename Sample>
tresult FilterDemoProcessor::processAudio(Sample** in, Sample** out, int32 numSamples, int32 numChannels)
{

	for (int channelIdx = 0; channelIdx < numChannels; channelIdx++)
	{
		Sample* channelInputBuffer  = in[channelIdx];
		Sample* channelOutputBuffer = out[channelIdx];
		double** channelOutputDelayBuffer  = outputDelayBuf[channelIdx];
		double** channelInputDelayBuffer   = inputDelayBuf[channelIdx];
		double x_n = 0; // Current Input  Sample
		double y_n = 0; // Current Output Sample

		//Filter Implementation (DFI Cascaded-Biquad IIR Filter)
		for (int sampleIdx = 0; sampleIdx < numSamples; sampleIdx++)
		{
			x_n = channelInputBuffer[sampleIdx];
			y_n = 0;
			for (int filterIdx = 0; filterIdx < noOfSecOrderFilters; filterIdx++)
			{
				// Get buffers/coeffs for current filter component (Current Biquad)
				double* biquadFbCoeffs_a = fbFilterCoeffsA[filterIdx];
				double* biquadFfCoeffs_b = ffwFilterCoeffsB[filterIdx];
				double* outputDelayBuffer = channelOutputDelayBuffer[filterIdx];
				double* inputDelayBuffer = channelInputDelayBuffer[filterIdx];

				//Compute output of current filter component
				y_n = biquadFfCoeffs_b[0] * x_n + biquadFfCoeffs_b[1] * inputDelayBuffer[1]  + biquadFfCoeffs_b[2] * inputDelayBuffer[2] -
					                              biquadFbCoeffs_a[1] * outputDelayBuffer[1] - biquadFbCoeffs_a[2] * outputDelayBuffer[2];

				//Update delay buffers for current filter component
				outputDelayBuffer[2] = outputDelayBuffer[1];
				outputDelayBuffer[1] = y_n;
				inputDelayBuffer[2] = inputDelayBuffer[1];
				inputDelayBuffer[1] = x_n;

				//Set output to input of next filter component
				x_n = y_n;
			}

			//Set output of last filter component to output of system
			channelOutputBuffer[sampleIdx] = y_n;
		}
	}

	return kResultOk;
}


//-----------------------------------------------------------------------------
tresult PLUGIN_API FilterDemoProcessor::process (Vst::ProcessData& data)
{
	//--- Read inputs parameter changes-----------
	if (data.inputParameterChanges)
	{
		int32 numParamsChanged = data.inputParameterChanges->getParameterCount ();
		for (int32 index = 0; index < numParamsChanged; index++)
		{
			Vst::IParamValueQueue* paramQueue =
			    data.inputParameterChanges->getParameterData (index);
			if (paramQueue)
			{
				double exponenetialFreqFactor = 1;
				double minCutoffFreq = 0.00306;
				Vst::ParamValue value;
				int32 sampleOffset;
				int32 numPoints = paramQueue->getPointCount ();
				switch (paramQueue->getParameterId ())
				{
					case FilterDemoParams::kNoZeroPairs:
						if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) ==
							kResultTrue)
							noOfZeroPairs = 10 * value; // Range = 0 - 10
						break;
					case FilterDemoParams::kNoPolePairs:
						if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) ==
							kResultTrue)
							noOfPolePairs = 10 * value; // Range = 0 - 10
						break;
					case FilterDemoParams::kResonatorQ:
						if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) ==
							kResultTrue)
							resonatorPoleRadius = value*0.2 + 0.8; // Range = 0.80 -  1.00
						break;
					case FilterDemoParams::kCutoffFreq:
						if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) ==
							kResultTrue)
							exponenetialFreqFactor = pow(2.0, (10.0 * value)); // 2 ^ 10*val -- Range = 1 - 1024 (exp scale for audio frequency)
							cutoffFreq = exponenetialFreqFactor * minCutoffFreq; // Range = minCutoffFreq - M_PI
							resonatorFreq = cutoffFreq - 0.09;// offset slightly
							zeroRadius = 1.0; // True Zero (On Unit Circle)
						break;
					case FilterDemoParams::kBypassId:
						if (paramQueue->getPoint (numPoints - 1, sampleOffset, value) ==
						    kResultTrue)
							mBypass = (value > 0.5f);
						break;
				}
				noOfSecOrderFilters = std::max(noOfPolePairs, noOfZeroPairs);
			}
		}

		//Feedback (B)
		//Resonator Poles
		for (int filterIdx = 0; filterIdx < noOfPolePairs; filterIdx++)
		{
			double* fbCoeffsi = fbFilterCoeffsA[filterIdx];
			fbCoeffsi[0] = 1;
			fbCoeffsi[1] = -2 * resonatorPoleRadius * cos(resonatorFreq);
			fbCoeffsi[2] = resonatorPoleRadius * resonatorPoleRadius;
		}


		//Feedforward (A)
		// Evenly distributed zeros f=cutoff to f=fs/2 (Nyquist)
		double interZeroFreqDelta = (M_PI - cutoffFreq) / noOfZeroPairs;
		for (int filterIdx = 0; filterIdx < noOfZeroPairs; filterIdx++)
		{
			double* ffwCoeffsi = ffwFilterCoeffsB[filterIdx];
			double zeroFreq = cutoffFreq + interZeroFreqDelta * filterIdx;
			ffwCoeffsi[0] = 1;
			ffwCoeffsi[1] = -2 * zeroRadius * cos(zeroFreq);
			ffwCoeffsi[2] = zeroRadius * zeroRadius;
		}

	}

	//--- Process Audio---------------------
	//--- ----------------------------------
	if (data.numInputs == 0 || data.numOutputs == 0)
	{
		// nothing to do
		return kResultOk;
	}

	if (data.numSamples > 0)
	{
		// Process Algorithm
		// Ex: algo.process (data.inputs[0].channelBuffers32, data.outputs[0].channelBuffers32,
		// data.numSamples);
		void** in = getChannelBuffersPointer(processSetup, data.inputs[0]);
		void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);

		if (Vst::kSample32 == processSetup.symbolicSampleSize)
		{
			processAudio<Vst::Sample32>((Vst::Sample32**)in, (Vst::Sample32**)out, data.numSamples, data.inputs[0].numChannels);
		}
		else
		{
			processAudio<Vst::Sample64>((Vst::Sample64**)in, (Vst::Sample64**)out, data.numSamples, data.inputs[0].numChannels);
		}
	}
	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API FilterDemoProcessor::setState (IBStream* state)
{
	if (!state)
		return kResultFalse;

	// called when we load a preset or project, the model has to be reloaded

	IBStreamer streamer (state, kLittleEndian);

	float savedParam1 = 0.f;
	if (streamer.readFloat (savedParam1) == false)
		return kResultFalse;

	int32 savedParam2 = 0;
	if (streamer.readInt32 (savedParam2) == false)
		return kResultFalse;

	int32 savedBypass = 0;
	if (streamer.readInt32 (savedBypass) == false)
		return kResultFalse;

	mBypass = savedBypass > 0;

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API FilterDemoProcessor::getState (IBStream* state)
{
	// here we need to save the model (preset or project)

	float toSaveParam1 = 0;
	int32 toSaveParam2 = 0;
	int32 toSaveBypass = mBypass ? 1 : 0;

	IBStreamer streamer (state, kLittleEndian);
	streamer.writeFloat (toSaveParam1);
	streamer.writeInt32 (toSaveParam2);
	streamer.writeInt32 (toSaveBypass);

	return kResultOk;
}

//------------------------------------------------------------------------
} // namespace
} // namespace Steinberg
