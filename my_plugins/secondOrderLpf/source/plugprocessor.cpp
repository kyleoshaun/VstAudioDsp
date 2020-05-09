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

#include <cmath>
#include <math.h>

#include "base/source/fstreamer.h"
#include "pluginterfaces/base/ibstream.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"
#include "public.sdk/source/vst/vstaudioprocessoralgo.h"

namespace Steinberg {
namespace SecondOrderLpf {

//-----------------------------------------------------------------------------
SecondOrderLpfProcessor::SecondOrderLpfProcessor ()
{
	// Init members
	inputDelayBuf = 0;
	outputDelayBuf = 0;
	mBypass = false;
	fbFilterCoeffsA = 0;
	ffFilterCoeffsB = 0;
	AL = BIQUAD_NO_OF_FB_COEFFS;
	BL = BIQUAD_NO_OF_FF_COEFFS;
	cutoffFreq = 0;

	// register its editor class
	setControllerClass (MyControllerUID);
}

//-----------------------------------------------------------------------------
tresult PLUGIN_API SecondOrderLpfProcessor::initialize (FUnknown* context)
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
tresult PLUGIN_API SecondOrderLpfProcessor::setBusArrangements (Vst::SpeakerArrangement* inputs,
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
tresult PLUGIN_API SecondOrderLpfProcessor::setupProcessing (Vst::ProcessSetup& setup)
{
	// here you get, with setup, information about:
	// sampleRate, processMode, maximum number of samples per audio block
	return AudioEffect::setupProcessing (setup);
}

//-----------------------------------------------------------------------------
tresult PLUGIN_API SecondOrderLpfProcessor::setActive (TBool state)
{
	Vst::SpeakerArrangement arr;
	if (getBusArrangement(Vst::kOutput, 0, arr) != kResultTrue)
		return kResultFalse;
	int32 numChannels = Vst::SpeakerArr::getChannelCount(arr);
	//TODO: This should be BIQUAD_NO_OF_FB_COEFFS
	size_t delayBufSize = processSetup.sampleRate * sizeof(Vst::Sample64) + 0.5; //max 1 sec delay (round up)

	if (state) // Initialize
	{
		inputDelayBuf  = (double**)std::malloc(sizeof(double*) * numChannels);
		outputDelayBuf = (double**) std::malloc(sizeof(double*) * numChannels);
		for (int channelIdx = 0; channelIdx < numChannels; channelIdx++)
		{
			inputDelayBuf[channelIdx]  = (double*)std::malloc(delayBufSize);
			outputDelayBuf[channelIdx] = (double*) std::malloc(delayBufSize);
			memset(outputDelayBuf[channelIdx], 0, delayBufSize);
		}

		//Initialize Filter as bypass
		fbFilterCoeffsA = (double*) std::malloc(sizeof(double) * BL);
		ffFilterCoeffsB = (double*)std::malloc(sizeof(double) * AL);
		memset(fbFilterCoeffsA, 0, sizeof(double) * BL);
		memset(ffFilterCoeffsB, 0, sizeof(double) * AL);
		ffFilterCoeffsB[0] = 1; //Bypass x[n]
		fbFilterCoeffsA[0] = 1; // a0 is always 1 in DSP theory
	}
	else // Release
	{
		if (outputDelayBuf)
		{
			for (int channelIdx = 0; channelIdx < numChannels; channelIdx++)
			{
				std::free(outputDelayBuf[channelIdx]);
			}
			std::free(outputDelayBuf);
			outputDelayBuf = 0;
		}
		if (inputDelayBuf)
		{
			for (int channelIdx = 0; channelIdx < numChannels; channelIdx++)
			{
				std::free(inputDelayBuf[channelIdx]);
			}
			std::free(inputDelayBuf);
			inputDelayBuf = 0;
		}
		if (fbFilterCoeffsA)
		{
			std::free(fbFilterCoeffsA);
			fbFilterCoeffsA = 0;
		}
		if (ffFilterCoeffsB)
		{
			std::free(ffFilterCoeffsB);
			ffFilterCoeffsB = 0;
		}

	}
	return AudioEffect::setActive (state);
}

//-----------------------------------------------------------------------------
template <typename Sample>
tresult SecondOrderLpfProcessor::processAudio(Sample** in, Sample** out, int32 numSamples, int32 numChannels)
{

	for (int channelIdx = 0; channelIdx < numChannels; channelIdx++)
	{
		Sample* channelInputBuffer  = in[channelIdx];
		Sample* channelOutputBuffer = out[channelIdx];
		double* channelinputDelayBuffer  = inputDelayBuf[channelIdx];
		double* channelOutputDelayBuffer = outputDelayBuf[channelIdx];
		double x_n = 0; // Current Input  Sample
		double y_n = 0; // Current Output Sample

		//Resonator Implementation (2nd order Feedback IIR Filter)
		for (int sampleIdx = 0; sampleIdx < numSamples; sampleIdx++)
		{
			x_n = channelInputBuffer[sampleIdx];
			y_n = 0;
			
			//Compute Output sample
			y_n = ffFilterCoeffsB[0] * x_n + ffFilterCoeffsB[1] * channelinputDelayBuffer[1]  + ffFilterCoeffsB[2] * channelinputDelayBuffer[2] -
				                             fbFilterCoeffsA[1] * channelOutputDelayBuffer[1] - fbFilterCoeffsA[2] * channelOutputDelayBuffer[2];
			
			//Update Delay Buffers
			channelinputDelayBuffer[2]  = channelinputDelayBuffer[1];
			channelinputDelayBuffer[1]  = x_n;
			channelOutputDelayBuffer[2] = channelOutputDelayBuffer[1];
			channelOutputDelayBuffer[1] = y_n;

			channelOutputBuffer[sampleIdx] = y_n;
		}
	}

	return kResultOk;
}


//-----------------------------------------------------------------------------
tresult PLUGIN_API SecondOrderLpfProcessor::process (Vst::ProcessData& data)
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
				Vst::ParamValue value;
				int32 sampleOffset;
				int32 numPoints = paramQueue->getPointCount ();
				switch (paramQueue->getParameterId ())
				{
					case SecondOrderLpfParams::kCutoffFreq: //resonatorFreq
						if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) ==
							kResultTrue)
							//cutoffFreq = M_PI * value; // value is input as 0-1 -> Map to F=0-Fs/2, f=0-1/2, w=0-pi
							cutoffFreq = 0.01 * pow(2.0, 8.295 * value); // Convert linear 0-1 input to exponential (base 2) 0.01-PI cutoff
						break;
					case SecondOrderLpfParams::kResonanceQ://resonatorQ
						if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) ==
							kResultTrue)
							resonanceQFactor = MIN_RESONANCE_Q_FACTOR + value * (MAX_RESONANCE_Q_FACTOR - MIN_RESONANCE_Q_FACTOR);
						break;
					case SecondOrderLpfParams::kBypassId:
						if (paramQueue->getPoint (numPoints - 1, sampleOffset, value) ==
						    kResultTrue)
							mBypass = (value > 0.5f);
						break;
				}
			}
		}

		//Coefficient parameter calculation
		double beta = 0.5 * (1.0 - sin(cutoffFreq) / (2.0 * resonanceQFactor)) / (1.0 + sin(cutoffFreq) / (2.0 * resonanceQFactor));
		double gamma = (0.5 + beta) * cos(cutoffFreq);

		//Feedback Coefficients (A)
		fbFilterCoeffsA[1] = -2.0 * gamma;
		fbFilterCoeffsA[2] = 2.0 * beta;

		//Feed-forward Coefficients (B)
		ffFilterCoeffsB[0] = (0.5 + beta - gamma) / 2.0;
		ffFilterCoeffsB[1] = 0.5 + beta - gamma;
		ffFilterCoeffsB[2] = ffFilterCoeffsB[0];

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
tresult PLUGIN_API SecondOrderLpfProcessor::setState (IBStream* state)
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
tresult PLUGIN_API SecondOrderLpfProcessor::getState (IBStream* state)
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
