// Copyright 2019 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Simple check that the spike generator behaves.

#include <algorithm>
#include <random>
#include <vector>

#include "Eigen/Core"
#include "eidos/audition/gammatone_filterbank_config.pb.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "eidos/audition/utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "third_party/audition/models/amt/meddis1986_hair_cell_synapse_model.h"
#include "third_party/audition/models/bruce_carney/bruce2018_spike_generator.h"
#include "third_party/audition/models/gammatone/gammatone_filterbank_model.h"

using testing::AllOf;
using testing::Each;
using testing::Gt;
using testing::Le;

namespace eidos {
namespace audition {
namespace {

constexpr double kSpontRate = 60.0;  // Spontaneous firing rate (spikes/sec).
constexpr double kAbsoluteRefractoryPeriod = 0.7E-3;  // sec
constexpr double kRelativeRefractoryPeriod = 0.6E-3;  // sec

}  // namespace

TEST(Bruce2018SpikeGeneratorTest, SingleFiberRandomStimulus) {
  constexpr double kNumSamples = 100000;  // Number of samples to test.
  constexpr double kSampleRate = 100000;  // 100 kHz
  constexpr double kMinRate = 11.0;       // Minimal firing rate.
  constexpr double kMaxRate = 153.0;      // Maximum firing rate.
  constexpr double kMaxRepeats = 6;  // Increase the initial stimulii buffer
                                     // size by (maximum) this factor.

  for (int n = 1; n <= kMaxRepeats; ++n) {
    const std::vector<double> &rates =
        GetRandomUniformVector(kNumSamples * n, kMinRate, kMaxRate);
    const std::vector<double> &spike_times = Bruce2018SpikeGenerator(
        rates, 1.0 / kSampleRate /* sample_period */, kSpontRate,
        kAbsoluteRefractoryPeriod, kRelativeRefractoryPeriod);
    EXPECT_LT(2, spike_times.size());  // At least two spikes.
    EXPECT_TRUE(std::is_sorted(spike_times.begin(), spike_times.end()));

    // Check that all the spike times are sane. Please note, we don't expect a
    // spike at time 0.
    const double stimulus_length = rates.size() / kSampleRate;
    ASSERT_THAT(spike_times, Each(AllOf(Gt(0.0), Le(stimulus_length))));
  }
}

TEST(Bruce2018SpikeGeneratorTest, InputsFromMeddis1986HairCell) {
  // Generate firing rates using Meddis' hair cell model given the random
  // stimulus.
  constexpr int kNumChannels = 251;
  constexpr int kNumSamples = 10000;  // Long enough so we get a good average.
  const Eigen::ArrayXXd input =
      Eigen::ArrayXXd::Random(kNumChannels, kNumSamples);
  constexpr int kSampleRate = 16000;  // Hz.
  Eigen::ArrayXXd rates;
  Meddis1986HairCellSynapse(input, kSampleRate,
                            false /* subtract spontaneous firing rate */,
                            &rates);
  EXPECT_EQ(kNumChannels, rates.rows());
  EXPECT_EQ(kNumSamples, rates.cols());

  // Run through the spike generator.
  Eigen::ArrayXXd spikes;
  Bruce2018SpikeGenerator(rates, 1.0 / kSampleRate /* sample_period */,
                          kSpontRate, kAbsoluteRefractoryPeriod,
                          kRelativeRefractoryPeriod, &spikes);
  EXPECT_EQ(kNumChannels, spikes.rows());
  EXPECT_EQ(kNumSamples, spikes.cols());
  EXPECT_FALSE(spikes.isZero());
  EXPECT_EQ(0.0, spikes.minCoeff());
  EXPECT_EQ(1.0, spikes.maxCoeff());

  // Expect these lower bounds on the number of spikes.
  EXPECT_LT(7000, (spikes.array() != 0).count());
}

TEST(Bruce2018SpikeGeneratorTest, InputsFromMaGammatoneAndMeddis1986HairCell) {
  // Given the random stimulus, approximate basilar membrane displacements using
  // gammatone filters followed by the Meddis' model producing the firing rates.
  constexpr int kNumChannels = 251;
  constexpr double kSampleRate = 16000.0;  // Hz.
  constexpr int kNumSamples = 20000;  // Long enough so we get a good average.
  constexpr int kNumEars = 1;
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);

  // Estimate basilar membrane displacements.
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(kSampleRate);
  GammatoneFilterbankConfig gammatone_config;
  gammatone_config.set_filter_type(GAMMATONE_FILTERBANK_COOKE_AND_MA);
  GammatoneFilterbankModel model(gammatone_config);
  model.Init(stimulus_config);
  Eigen::ArrayXXd bm_output;
  model.ProcessSegment(input, &bm_output);
  const ChannelProperties &channel_properties = model.GetChannelProperties();
  EXPECT_EQ(channel_properties.NumChannels(), bm_output.rows());
  EXPECT_EQ(kNumSamples, bm_output.cols());

  // Estimate firing rates.
  Eigen::ArrayXXd rates;
  Meddis1986HairCellSynapse(bm_output, kSampleRate,
                            false /* subtract spontaneous firing rate */,
                            &rates);
  EXPECT_EQ(kNumChannels, rates.rows());
  EXPECT_EQ(kNumSamples, rates.cols());

  // Run through the spike generator.
  Eigen::ArrayXXd spikes;
  Bruce2018SpikeGenerator(rates, 1.0 / kSampleRate /* sample_period */,
                          kSpontRate, kAbsoluteRefractoryPeriod,
                          kRelativeRefractoryPeriod, &spikes);
  EXPECT_EQ(kNumChannels, spikes.rows());
  EXPECT_EQ(kNumSamples, spikes.cols());
  EXPECT_FALSE(spikes.isZero());
  EXPECT_EQ(0.0, spikes.minCoeff());
  EXPECT_EQ(1.0, spikes.maxCoeff());

  // Expect these lower bounds on the number of spikes.
  EXPECT_LT(2400, (spikes.array() != 0).count());
}

TEST(Bruce2018SpikeGeneratorTest, SimpleModelInterfaceRandomInputs) {
  // Generate random firing rate inputs in the range (kMinRate, kMaxRate).
  constexpr double kNumChannels = 50;
  constexpr double kNumSamples = 200000;  // Number of samples to test.
  constexpr double kMinRate = 11.0;       // Minimal firing rate.
  constexpr double kMaxRate = 153.0;      // Maximum firing rate.
  const Eigen::ArrayXXd input =
      kMinRate + (Eigen::ArrayXXd::Random(
          kNumChannels, kNumSamples) * 0.5 + 0.5) + (kMaxRate - kMinRate);

  // Run the model.
  constexpr double kSampleRate = 100000;  // 100 kHz
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(kSampleRate);
  Bruce2018SpikeGeneratorModel model;
  Eigen::ArrayXXd spikes;
  model.Init(stimulus_config);
  model.ProcessSegment(input, &spikes);

  EXPECT_EQ(kNumChannels, spikes.rows());
  EXPECT_EQ(kNumSamples, spikes.cols());
  EXPECT_FALSE(spikes.isZero());
  EXPECT_EQ(0.0, spikes.minCoeff());
  EXPECT_EQ(1.0, spikes.maxCoeff());

  // Expect these lower bounds on the number of spikes.
  constexpr int kMinNumSpikes = 8000;
  EXPECT_LT(kMinNumSpikes, (spikes.array() != 0).count());
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
