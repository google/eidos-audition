# Eidos (Audition)

Eidos (Audition) TODO list.

## Version: pre-alpha
TODO for the pre-alpha version:

### Miscellaneous:
  - [X] Implement basic cochlear feature visualizer.
  - [X] Implement time-smoothing the output for NAP and IHC excitations. For
        Baumgarte's model we can achieve something similar with heavy
        downsampling.
  - [X] Visualizer end-to-end test.
  - [X] Fix the bug with visualizer: specifying pyplot imshow's extent
        screws up orientation of the y axis.
  - [X] Implement Jackson's Poisson process-based spike generator.
  - [X] Remove ears from the auditory model interface (monaural only).
  - [X] Use doubles in the output of auditory model interface instead of floats.
  - [X] Implement resampler (C++ version of Matlab/Octave polyphase resample()).
  - [X] Migrate Jackson's spike generator to use Eigen throughout: Decided not
        to do this for now, because test mocks need to be extended to Eigen.
  - [X] Support .npz output in feature extractor. This can contain the outputs
        from all the stages in the processing.
  - [X] Visualizer: support for multiple plots in one figure.
  - [X] Visualizer: Show additional labels with the actual values of frequency
        channels.
  - [X] Check Sumner's synaptic adaptation of inner hair cells known as Meddis'2000
        in DSAM. It is a stochastic model, and with the exception of a probabilistic
        output option, is the same as the DSAM Meddis 95 module.
  - [X] Support parallel processing of the input file lists.
  - [X] Decide where framer (downsampler) is needed on the C++ side in the
        cochlear feature extractor. Existing visualizer implementation
        may be enough.
  - [X] Implement downsample step functionality for all the models (move this out
        out of Baumgarte).
  - [ ] Investigate reimplementing parallelization using absl::Barrier API.
  - [ ] Resampler (Matlab-like API) is currently broken. Even downsampling by a
        factor of 2 does not produce a coherent signal. Investigate.

### API:
  - [X] Remove dependency on feature extractor config in each model.
  - [X] Add support for composite auditory models consisting of cascade of
        several models.
  - [X] Support running the extractor from a proto file.
  - [X] Always sort the center frequencies in ascending order and re-arrange
        the outputs correspondingly.
  - [X] Compute and pass (between the filters in the cascade) the center
        frequency parameters.
  - [X] Figure out how to pass model-specific configuration around. Migrate
        via constructors or via Init() (in which case remove abstract Init()
        from the base class).
  - [X] Support passing in the pipeline proto as textual string to the feature
        extractor.
  - [X] Support configuring the number of channels.
  - [X] Support configuring the channel ranges.
  - [X] Support passing stimulus config proto as a parameter in feature
        extractor tool.
  - [X] Implement pipeline proto-specific end-to-end test.
  - [ ] Add more pipeline configurations under "configs".

### CARFAC:
  - [X] Implement very initial version of IHC activations using CARFAC.
  - [X] Implement end-to-end shell test for CARFAC feature extractor.
  - [X] Add configuration.
  - [X] Extend tests to check for various configurations.
  - [ ] Fix CARFAC includes. In client code the headers can only be included
        as "cpp/*.h". Including under "carfac/*.h" will be ideal.
  - [ ] Add the BM and OHC signals to the output of CARFAC extractor.
  - [ ] Add check against the golden data.
  - [ ] Check if the filter produces more detailed results for higher sampling
        rates.

### Zilany, et. al. (2014):
  - [X] Possibly integrate the Zilany et. al. (2014) IHC model.
  - [X] IHC: Need different scaling of the interval for characteristic
        frequencies, e.g. similar to CARFAC's ERB scale.
  - [X] Possibly change the Nepali waveform for a smaller speech file --
        at the moment it contains a full sentence which makes tests slow.
  - [X] Configurable number of channels/partitions.
  - [ ] Configurable min/max frequency limits.
  - [ ] IHC: Very slow. Optimize Eigen operations.
  - [ ] IHC: Implement checks against the output of the reference model.
  - [ ] Check Zilany 2014 synapse model.
  - [ ] Fix configuration: spikes_zilany2014_bruce2018 (MacOS X failure).
  - [ ] Fix configuration: spikes_zilany2014_jackson (Linux failure).
  - [ ] Fix configuration: spikes_zilany2014_zhang2001 (MacOS X failure).

### Ray Meddis's IHC synapse model (1986):
  - [X] Implement Ray Meddis' inner hair cell synapse model with tests.
  - [X] Integrate Gammatone and/or BM portion of Baumgarte's model with the
        test for Meddis synapse.
  - [X] Extend Meddis' test using the real BM inputs.
  - [X] Integrate with Jackson's spike generator and visualizer.
  - [X] Possibly rename this model to Meddis'86 (or Meddis'95?).

### Bruce, et. al. (2018) IHC and synapse model:
  - [X] Implement spike generator from Bruce, et. al. (2018).
  - [X] Implement inner hair cell synapse model from Bruce, et. al. (2018).
  - [X] Add "real" BM inputs (Gammatone) to check the synapse model.
  - [ ] The spike generator should take number of trials into account. Currently
        it doesn't.
  - [ ] The spike generator should have channel-specific spontaneous firing
        rate and refractory period parameters.
  - [ ] Make sure the IHC synapse model is CF-specific.
  - [ ] Check against the Matlab-generated reference data.
  - [ ] Implement FFGN version from Yingchun Zhou and Stilian Stoev.

### Zhang and Carney (2001):
  - [X] Investigate simple spike generator.
  - [ ] Check how to turn number of trials into the number of fibers, hence
        turning the output from PSTH into an estimate of total AN discharge
        from a large population of fibers.

### Baumgarte (2000):
  - [X] Implement peripheral ear model from Frank Baumgarte's thesis.
  - [X] Configurable number of channels/partitions.
  - [ ] Configurable min/max frequency limits.

### Gammatone filters:
  - [X] Implement gammatone filter (simulating BM) with tests.
  - [X] Make particular Gammatone filter selection (Slaney/Ma) configurable.
  - [X] Configurable number of channels/partitions.
  - [X] Configurable min/max frequency limits.
  - [X] Slaney: The filter bank contains N channels that extend from half the
        sampling rate to "low_freq". Make sure that this mode is enabled *only*
	if the "high_freq" is not in the stimulus configuration.
  - [X] Create two simple configuration pipelines (for Slaney's and Ma's
        filters).
