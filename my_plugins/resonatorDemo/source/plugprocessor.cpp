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
namespace ResonatorDemo {

//-----------------------------------------------------------------------------
ResonatorDemoProcessor::ResonatorDemoProcessor ()
{
	// Init members
	delayBuf = 0;
	delayBufIdx = 0;
	mBypass = false;
	resonatorCoeffsB = 0;
	resonatorBL = 3;
	resonatorA0 = 1;
	resonatorPoleRadius = 0;
	resonatorFreq = 0;

	// register its editor class
	setControllerClass (MyControllerUID);
}

//-----------------------------------------------------------------------------
tresult PLUGIN_API ResonatorDemoProcessor::initialize (FUnknown* context)
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
tresult PLUGIN_API ResonatorDemoProcessor::setBusArrangements (Vst::SpeakerArrangement* inputs,
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
tresult PLUGIN_API ResonatorDemoProcessor::setupProcessing (Vst::ProcessSetup& setup)
{
	// here you get, with setup, information about:
	// sampleRate, processMode, maximum number of samples per audio block
	return AudioEffect::setupProcessing (setup);
}

//-----------------------------------------------------------------------------
tresult PLUGIN_API ResonatorDemoProcessor::setActive (TBool state)
{
	Vst::SpeakerArrangement arr;
	if (getBusArrangement(Vst::kOutput, 0, arr) != kResultTrue)
		return kResultFalse;
	int32 numChannels = Vst::SpeakerArr::getChannelCount(arr);
	size_t delayBufSize = processSetup.sampleRate * sizeof(Vst::Sample64) + 0.5; //max 1 sec delay (round up)

	if (state) // Initialize
	{
		delayBuf = (double**) std::malloc(sizeof(double*) * numChannels);
		for (int channelIdx = 0; channelIdx < numChannels; channelIdx++)
		{
			delayBuf[channelIdx] = (double*) std::malloc(delayBufSize);
			memset(delayBuf[channelIdx], 0, delayBufSize);
		}
		delayBufIdx = 0;

		//Initialize resonator as bypass
		resonatorA0 = 1;
		resonatorCoeffsB = (double*) std::malloc(sizeof(double) * resonatorBL);
		memset(resonatorCoeffsB, 0, sizeof(double) * resonatorBL);
	}
	else // Release
	{
		if (delayBuf)
		{
			for (int channelIdx = 0; channelIdx < numChannels; channelIdx++)
			{
				std::free(delayBuf[channelIdx]);
			}
			std::free(delayBuf);
			delayBuf = 0;
			std::free(resonatorCoeffsB);
			resonatorCoeffsB = 0;
		}
	}
	return AudioEffect::setActive (state);
}

//-----------------------------------------------------------------------------
template <typename Sample>
tresult ResonatorDemoProcessor::processAudio(Sample** in, Sample** out, int32 numSamples, int32 numChannels)
{

	for (int channelIdx = 0; channelIdx < numChannels; channelIdx++)
	{
		Sample* channelInputBuffer  = in[channelIdx];
		Sample* channelOutputBuffer = out[channelIdx];
		double* channelDelayBuffer  = delayBuf[channelIdx];
		double x_n = 0; // Current Input  Sample
		double y_n = 0; // Current Output Sample

		//Resonator Implementation (2nd order Feedback IIR Filter)
		for (int sampleIdx = 0; sampleIdx < numSamples; sampleIdx++)
		{
			x_n = channelInputBuffer[sampleIdx];
			y_n = 0;
			y_n = resonatorA0 * x_n - resonatorCoeffsB[1] * channelDelayBuffer[1] - resonatorCoeffsB[2] * channelDelayBuffer[2];
			channelDelayBuffer[2] = channelDelayBuffer[1];
			channelDelayBuffer[1] = y_n;

			channelOutputBuffer[sampleIdx] = y_n;
		}
	}

	return kResultOk;
}


//-----------------------------------------------------------------------------
tresult PLUGIN_API ResonatorDemoProcessor::process (Vst::ProcessData& data)
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
					case ResonatorDemoParams::kResonatorFreq: //resonatorFreq
						if (paramQueue->getPoint (numPoints - 1, sampleOffset, value) ==
						    kResultTrue)
							resonatorFreq = M_PI * value; // value is input as 0-1 -> Map to F=0-Fs/2, f=0-1/2, w=0-pi
						break;
					case ResonatorDemoParams::kResonatorQ://resonatorQ
						if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) ==
							kResultTrue)
							resonatorPoleRadius = value*1.0 + 0; // Range = 0.00 -  1.00
						break;
					case ResonatorDemoParams::kBypassId:
						if (paramQueue->getPoint (numPoints - 1, sampleOffset, value) ==
						    kResultTrue)
							mBypass = (value > 0.5f);
						break;
				}
			}
		}
		resonatorCoeffsB[0] = 1; //meaningless (b0 is always 1 in  DSP theory)
		resonatorCoeffsB[1] = -2 * resonatorPoleRadius * cos(resonatorFreq); 
		resonatorCoeffsB[2] = resonatorPoleRadius * resonatorPoleRadius;
		resonatorA0 = 1.0;

		//Algorithm to normalize peak magnitude to 1:
		/*(resonatorCoeffsB[2] == 0) ?
		1.0 :
		(1.0 - resonatorCoeffsB[2]) * sqrt(1.0 - (resonatorCoeffsB[1] * resonatorCoeffsB[1]) / (4.0 * resonatorCoeffsB[2]));*/

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
tresult PLUGIN_API ResonatorDemoProcessor::setState (IBStream* state)
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
tresult PLUGIN_API ResonatorDemoProcessor::getState (IBStream* state)
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
