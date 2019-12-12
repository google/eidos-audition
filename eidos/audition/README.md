# Models of Auditory Periphery

This library represents a simple pipeline integrating various models
of auditory periphery developed over the years in academia. The pipeline
can be used to combine models representing various stages of auditory
perception for extracting features that can be used by machine learning
algorithms.

## Directory Structure

Directory layout for the library looks as follows

```
|-- WORKSPACE
|-- build
|-- eidos
    |-- audition
    |-- stubs
    |-- utils
|-- third_party
    |-- audition
        |-- models
```

*  The [WORKSPACE](../../WORKSPACE) file contains Bazel workspace definitions
   for the external dependencies, such as [protocol buffers](https://developers.google.com/protocol-buffers).
*  The [build](../../build) directory contains the Bazel build files for
   external dependencies as well as simple unit tests for making sure the
   dependencies are built correctly.
*  Directory [eidos/audition](../../eidos/audition) is the main directory that
   contains the pipeline and the corresponding tools.
*  Directory [third_party/audition/models](../../third_party/audition/models)
   contains the implementations of the auditory models derived from the
   third-party code. Each subdirectory refers either to an individual model or
   to a group of models released by a certain lab or an open-source initiative.
   The original models have been modified to operate together in a single
   pipeline. Please see the documentation in individual directories for model
   information as well as licensing.

## Models

The library includes miscellaneous models of various stages of the processing
by the auditory periphery:

*  The [model](../../third_party/audition/models/baumgarte/baumgarte_model.h) of
   auditory periphery originally implemented by Frank Baumgarte for his PhD
   thesis [8]. The peripheral ear model is based on the structure of Eberhard
   Zwicker's "Analogmodell" consisting of analog electrical elements and makes
   use of Wave Digital Filters (WDF) algorithm originally described by Wolfgang
   Peisl in his PhD thesis (see [9]). The model provides estimates for basilar
   membrane displacements and inner and outer hair cell action potentials.
*  Cascade of Asymmetric Resonators with Fast-Acting Compression (CARFAC)
   [model](carfac_model.h) by Richard F. Lyon [6,7]. This model provides
   estimates for basilar membrane displacements and transmembrane voltages
   across inner hair cells and fine tuning of nonlinearities by adaptive gain
   control (AGC) filters.
*  Gammatone filterbank
   [model](../../third_party/audition/models/gammatone/gammatone_filterbank_model.h)
   based on the original filterbank by Malcolm Slaney [1] and a faster
   version developed by Ning Ma from University of Sheffield. This model can
   be used for providing estimates for basilar membrane displacements.
*  The
   [model](../../third_party/audition/models/amt/meddis1986_hair_cell_synapse_model.h)
   of inner hair synapse originally developed by Ray Meddis [4,5].  Multiple
   implementations of this model exist in various toolkits, such as Auditory
   Modeling Toolkit (AMT), Malcolm Slaney's Auditory Toolbox and Development System
   for Auditory Modelling (DSAM).
*  Synaptic inner hair cell
   [model](../../third_party/audition/models/dsam/sumner2002_hair_cell_synapse_model.h)
   by Sumner, C. J, Lopez-Poveda, E. A., O'Mard, L. P. and Meddis, R. [17].
*  A
   [model](../../third_party/audition/models/bruce_carney/zilany2014_ihc_model.h)
   of auditory periphery by Zilany, M.S.A., Bruce, I.C., Nelson, P.C. and
   Carney, L.H. [10,11,12]. The derived implementation provided by this library
   supports estimates for hair cell transmembrane potentials, excluding the full
   synaptic model available in the original implementation.
*  Hair cell synaptic
   [model](../../third_party/audition/models/bruce_carney/bruce2018_hair_cell_synapse_model.h)
   and the accompanying spike generation
   [model](../../third_party/audition/models/bruce_carney/bruce2018_spike_generator.h)
   developed by Ian Bruce, et. al. [14].
*  The spike generation [model](../../third_party/audition/models/bruce_carney/zhang2001_spike_generator.h)
   originally developed by Xuedong Zhang, et. al. [16].
*  The spike generation [model](../../third_party/audition/models/bruce_carney/jackson_spike_generator.h)
   from Jackson, B. S. and Carney L. H. [15].

## Building and Testing

The main prerequisite required for building the pipeline is the
[Bazel](https://bazel.build/) build system. Please follow the
documentation to install Bazel on your system. The code has been
developed and tested on Linux and Mac OS X and should work out of the
box on these operating systems.

To build all the targets, run the following build command from the
root directory

```shell
> bazel build -c opt ...
```

Following will run all the tests

```shell
> bazel test -c opt ...
```

## Running

### Auditory Feature Extractor

Auditory feature [extractor](auditory_feature_extractor_main.cc) is a simple
tool for extracting auditory representations from the audio. To build it, run

```shell
> bazel build -c opt eidos/audition:auditory_feature_extractor
```

The binary can be found in
`bazel-bin/eidos/audition/auditory_feature_extractor`. To get the list of
available command-line parameters please run this tool with `--helpshort`
option.

Following example modes of operation are available:

*  To process the input waveform using a list of model names supplied on the
   command-line:

   ```shell
   > bazel-bin/eidos/audition/auditory_feature_extractor \
      --models GAMMATONE_SLANEY,SUMNER_SYNAPSE_2002 \
      --waveform_file input.wav --output_file output.npy
   ```

*  To specify the processing pipeline configuration as a string protocol buffer
   run:

    ```shell
    > bazel-bin/eidos/audition/auditory_feature_extractor \
       --config_proto_contents "pipeline { models { model_type: MODEL_BAUMGARTE } }" \
       --waveform_file input.wav --output_file output.npy
    ```

*  To specify the processing pipeline configuration as a text protocol buffer
   file run:

   ```shell
   > bazel-bin/eidos/audition/auditory_feature_extractor \
      --config_proto_file pipeline.textproto \
      --waveform_file input.wav \
      --output_file_format npz --output_file output.npz
   ```

*  To process multiple audio files specified as a file list run:

   ```shell
   > bazel-bin/eidos/audition/auditory_feature_extractor \
      --models BAUMGARTE,SUMNER_SYNAPSE_2002 \
      --waveform_file_list file_list.txt --output_dir /tmp \
      --output_file_format npz
   ```

   Use `--num_tasks` flag to specify the number of workers if parallel
   processing is required.

### Visualizer

[Visualizer](visualize_auditory_signals.py) is a simple Python tool for
visualing the auditory responses extracted using auditory feature
[extractor](auditory_feature_extractor_main.cc).  Build it using

```shell
> bazel build -c opt eidos/audition:visualize_auditory_signals
```

The binary can be found in
`bazel-bin/eidos/audition/visualize_auditory_signals`. Please use the
`--helpshort` flag to get the list of available options.

The visualizer requires an input file (provided by `--input_signal_file` flag)
produced by the auditory feature extractor in either Numpy `.npy` or `.npz`
format. Example usage:

```shell
> bazel-bin/eidos/audition/visualize_auditory_signals \
   --input_signal_file output.npz
```

---

## References

1.  Slaney, M. (1993): ["An Efficient Implementation of the Patterson-Holdsworth
    Auditory Filter Bank"](https://engineering.purdue.edu/~malcolm/apple/tr35/PattersonsEar.pdf), Apple Technical Report #35.
2.  Moore, B. C. J. and Glasberg, B. R (1983): ["Suggested formulae for
    calculating auditory-filter bandwidths and excitation patterns."](https://asa.scitation.org/doi/abs/10.1121/1.389861),
    J. Acoust. Soc. Am. 74: 750-753.
3.  Glasberg, B. R. and Moore, B. C. J. (1990): ["Derivation of auditory filter
    shapes from notched-noise data."](http://www.academia.edu/download/42800388/Derivation_of_auditory_filter_shapes_fro20160218-28794-178jcg8.pdf), Hearing Research, 47: 103-138.
4.  Ray Meddis, Michael J. Hewitt, and Trevor M. Shackleton (1990): ["Implementation
    details of a computation model of the inner hair cell auditory nerve
    synapse."](https://asa.scitation.org/doi/pdf/10.1121/1.399379), The Journal of the Acoustical Society of America 87, 1813.
5.  Ray Meddis (1986): ["Simulation of mechanical to neural transduction in the
    auditory receptor."](https://www.researchgate.net/profile/Ray_Meddis/publication/20286351_Simulation_of_mechanical_to_neural_transduction_in_the_auditory_receptor/links/0c96052a1db13b03e0000000.pdf), Journal of the Acoustical Society of America 79(3),
    702--711.
6.  Lyon, R. (2017): "The CARFAC Digital Cochlear Model.", In ["Human
    and Machine Hearing: Extracting Meaning from Sound."](http://dicklyon.com/Lyon_Hearing_book_companion_color.pdf),
    Chapter 15, (pp. 293-298), Cambridge: Cambridge University Press.
7.  Lyon, R. (2011): ["Using a Cascade of Asymmetric Resonators with
    Fast-Acting Compression as a Cochlear Model for Machine-Hearing
    Applications."](https://ai.google/research/pubs/pub37215.pdf), Autumn Meeting of the Acoustical Society of Japan
    (2011), pp. 509-512.
8.  Baumgarte, F. (2000): ["Ein psychophysiologisches Gehoermodell zur
    Nachbildung von Wahrnehmungsschwellen fuer die Audiocodierung."](ftp://ftp.tnt.uni-hannover.de/pub/papers/2000/Diss-FB.pdf), PhD Dissertation, University of Hannover, February.
9.  Zwicker, E. and Peisl, W. (1990): ["Cochlear preprocessing in analog models,
    in digital models and in human inner ear"](https://www.sciencedirect.com/science/article/pii/037859559090081Y),
    Hearing Research, vol. 44 (no. 2-3): March, pp. 209-216.
10.  Zilany, M.S.A., Bruce, I.C., Ibrahim, R.A., and Carney,
    L.H. (2013): ["Improved parameters and expanded simulation
    options for a model of the auditory periphery."](http://www.ece.mcmaster.ca/~ibruce/papers/zilany_aro2013.pdf), in Abstracts of the 36th ARO Midwinter Research Meeting.
11.  Zilany, M.S.A., Bruce, I.C., Nelson, P.C., and Carney, L.H. (2009): ["A
     phenomenological model of the synapse between the inner hair cell and
     auditory nerve: Long-term adaptation with power-law dynamics."](http://www.ece.mcmaster.ca/~ibruce/papers/Bruce_ARO2017_poster.pdf), Journal of the Acoustical Society of America 126(5): 2390-2412.
12.  Zilany, M.S.A., Bruce, I.C., and Carney, L.H. (2014): ["Updated parameters
     and expanded simulation options for a model of the auditory periphery."](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3985897), Journal of the Acoustical Society of America.
13.  Ibrahim, R. A., and Bruce, I. C. (2010): ["Effects of peripheral tuning
    on the auditory nerve's representation of speech envelope and
    temporal fine structure cues."](https://mimm.mcmaster.ca/publications/pdfs/IbrahimBruce_ISH2010.pdf) in The Neurophysiological Bases of Auditory Perception, eds. E. A. Lopez-Poveda and A. R. Palmer and R. Meddis, Springer, NY, pp. 429-438.
14.  Bruce, I.C., Erfani, Y., and Zilany, M.S.A. (2018): ["A
     phenomenological model of the synapse between the inner hair cell
     and auditory nerve: Implications of limited neurotransmitter
     release sites."](http://www.ece.mcmaster.ca/~ibruce/papers/Bruce_ARO2017_poster.pdf), Hearing research, 360, 40--54,      (Special Issue on "Computational Models in Hearing").
15.  Jackson, B.S. and Carney L. H. (2005): ["The spontaneous-rate histogram of the
     auditory nerve can be explained by only two or three spontaneous rates
     and long-range dependence."](https://link.springer.com/article/10.1007/s10162-005-5045-6), J. Assoc. Res. Otolaryngol. 6:148-159.
16.  Zhang X., Heinz M.G., Bruce I.C., Carney L.H. (2001): ["A phenomenological
     model for the responses of auditorynerve fibers: I. Nonlinear tuning
     with compression and suppression."](https://mimm.mcmaster.ca/publications/pdfs/1.1336503.pdf), J. Acoust. Soc. Am. 109:648-670.
17.  Sumner, C. J, Lopez-Poveda, E. A., O'Mard, L. P. and Meddis, R. (2002): ["A revised model of the inner-hair      cell and auditory-nerve complex."](https://web.stanford.edu/group/boahen/proj/Sumner_Meddis_2002.pdf), The Journal of the Acoustical Society of America (JASA), vol.111, no.5, pp. 2178--2188.

---

## See also

The following lists are not nearly exhaustive, but should provide a useful and,
in certain cases, crucial references.

### Books and mailing lists

*  Lyon, R. (2017): [“Human and Machine Hearing: Extracting Meaning from Sound.”](https://www.amazon.co.uk/Human-Machine-Hearing-Extracting-Meaning/dp/1107007534), Cambridge: Cambridge University Press.
*  Meddis R., Lopez-Poveda E., Fay R. R. and Popper, A. N. (eds) 2010: [“Computational Models of the Auditory System”](https://www.amazon.com/Computational-Auditory-Springer-Handbook-Research/dp/144191370X/ref=sr_1_14), Series: "Springer Handbook of Auditory Research", Book 35, Springer.
*  Manley G. A., Gummer A. W., Popper A. N. and Fay, R. R. (eds) 2017: [“Understanding the Cochlea”](https://www.amazon.com/Understanding-Springer-Handbook-Auditory-Research/dp/3319520717/), Series: "Springer Handbook of Auditory Research", Book 62, Springer.
*  [AUDITORY](http://www.auditory.org/) list.

   An email list for the discussion of organizational aspects of auditory
   perception. The list was created in 1992 by Professor Albert S. Bregman of the
   McGill University Department of Psychology. It is currently maintained by [Dan
   Ellis](http://www.ee.columbia.edu/~dpwe/), who also runs that web site.

### Software

*  [DSAM](http://dsam.org.uk/): The Development System for Auditory Modelling (DSAM)

   The Development System for Auditory Modelling (DSAM) created in University of Essex
   is a computational library designed specifically for producing time-sampled auditory
   system simulations.

*  [AMT](http://amtoolbox.sourceforge.net/): The Auditory Modeling Toolbox

   The Auditory Modeling Toolbox (AMT) from Acoustics Research Institute (ARI) of
   Austrian Academy of Sciences is a Matlab/Octave toolbox intended to
   serve as a common ground for reproducible research in auditory modeling.

*  [Auditory Toolbox](https://engineering.purdue.edu/~malcolm/interval/1998-010) by Malcolm Slaney

   Collection of tools that implement several popular auditory models in Matlab.

*  [Auditory Models](https://www.urmc.rochester.edu/labs/carney/publications-code/auditory-models.aspx)
   from Laurel Carney's Lab.
