Documentation of peripheral ear model software

(c) Frank Baumgarte, 1999, 2000
Email: baumgart@tnt.uni-hannover.de
WWW:   http://www.tnt.uni-hannover.de/~baumgart

Acknowledgements:
  The peripheral ear model is based on the structure of
  Eberhard Zwicker's "Analogmodell". That model consists
  of analog electrical elements.

  Wolfgang Peisl described the basic WDF algorithm in his
  PhD thesis:
  W. Peisl: "Beschreibung aktiver nichtlinearer Effekte der
  peripheren Schallverarbeitung des Gehoers durch ein
  Rechnermodell", (in German), Dissertation,
  Technical University of Munich

This software is originally available from
http://www.tnt.uni-hannover.de/~baumgart/peri_ear_v0.1.tgz

A comprehensive documentation of the software provided here
is contained in the dissertation:
F. Baumgarte: "Ein psychophysiologisches Gehoermodell zur
Nachbildung von Wahrnehmungsschwellen fuer die Audiocodierung",
(in German), dissertation, University of Hannover, 2000.
Available online via:
ftp://ftp.tnt.uni-hannover.de/pub/papers/2000/Diss-FB.pdf

This software is written for UNIX platforms using ANSI-C.
MS-Windows operating systems are currently not supported.
However, you should be able to compile and link the C source
files within an appropriate Windows software development tool.

UNIX Installation:

1) Move the zipped tar-file "peri_ear_v0.1.tgz" to the directory
   in which you want to install the software.
2) Excute the following UNIX shell command to extract all files:
	tar -xzvf peri_ear_v0.1.tgz
3) Go to subdirectory "peri_ear_v0.1/" and execute
        make
   The executable will be written to "bin/peri_ear"

Running the program:
The filenames of an input and output file must be specified on
the command line when starting the program in a UNIX shell:
	peri_ear INPUT_FILENAME OUTPUT_FILENAME

INPUT_FILENAME is the name of an audio signal file
OUTPUT_FILENAME is the name of the model output signal filename

The standard input file is a raw (*.raw) audio file without
header. The default parameter setting of the program assumes
a sampling rate of 100 kHz and a 16 bit integer representation
of each sample. (It is recommended that the sampling rate is
it least 4 times the highest signal frequency).

The sound pressure level of the input audio signal is assumed to
be SPL = 10 * log10(sum(x*x)). Where sum(x*x) means summation
of all squared signal samples x.

A more convenient way to deal with audio files of different
formats is to use a header which provides the program with the
necessary parameters. Examples of such formats are the WAVE (*.wav)
format or AU (*.au) format. These formats can be processed by this
software if you define AFSP in the header file "model.h". Additionally
you need another software package which provides the file reading
functions. It can be obtained from

current version:        AFsp-v4r2.tar.gz
original site:          ftp://ftp.TSP.ECE.McGill.CA/pub/AFsp
                        http://www.TSP.ECE.McGill.CA/software
mirror site:            ftp://ftp.tnt.uni-hannover.de/pub/audio/AFsp

Two shell environment variables should be defined:
AUDIO_INCLUDE_PATH path to "libtsp.a"
AUDIO_LIBRARY_PATH path to "libtsp.h"

The file "Makefile" of the peri_ear package must be modified to
incorporate these environment variables (see "Makefile" for details).

The output signals of the model are written to the output file
"OUTPUT_FILENAME". With the default parameter setting the modeled
inner hair cell excitation signal at 251 equidistant locations in
the cochlea are generated. The default sampling rate of the output
is equal to the input sampling rate divided by 10. The samples are
written in 32 bit binary floating point format in consecutive blocks
of 251 samples each starting with model section 1 with increasing
section numbers.

The output as well as the input byte ordering can be swapped, if
necessary. The swapping is controlled by defining OUTPUT_BYTESWAP
or INPUT_BYTESWAP in "model.h". This might be necessary when using
machines with different "native" byte ordering like PC as opposed
to SUN and SGI.


Parameters:
The current version does not support parameter setting in the
command line. However, parameters can be changed in the source code.
The main model parameters are described here:

Name		File		Description
---------------------------------------------------------------------
DOWNSAMPLE	param.h		downsampling factor for output signal
FS_DEFAULT	param.h		default input sampling rate
MAX_BARK	param.h		highest critical-band rate covered by
				the model
DELTA_Z		param.h		spectral (spatial) resolution of model
				as fraction of critical-band width

inp_sel		get_para.c	audio signal is input at outer ear or
				at oval window
feedback	get_para.c	active or passive model (without outer
				hair cell model)
nonlin		get_para.c	characteristic of nonlinearity
coupling	get_para.c	coupling of neighboring model section


Testing:
A verification of the installation can be done by comparison with
example output files. These outputs can be generated by running the
model using following file names:

	bin/peri_ear examples/T1kL60.raw examples/T1kL60_tst.exc
	bin/peri_ear examples/imp2.raw   examples/imp2_tst.exc

Check the identity of examples/T1kL60_tst.exc and examples/T1kL60.exc
or examples/imp2_tst.exc and examples/imp2.exc.

The input file examples/T1kL60.* contains a 1 kHz sinusoid with 60 dB SPL.
The input file examples/imp2.* contains an impuls of amplitude 1000 at
t=0 and an impulse of -1000 at t=0.05 s. Both files are given in RAW
format (headerless) and AU format.

The output signals can be visualized using MatLab. Two MatLab figures
of both output signal files are provided in TIF format. The figures show
the magnitude output with zero amplitude in deep blue and larger
amplitudes gradually becoming red. Only the first 10 ms are shown.
The figures were generated using the following MatLab commands:

	[f_ptr, message] = fopen(example/T1kL60.exc,'rb','ieee-be');
	exc = fread(f_ptr, [251 1000], 'float32');
	fclose (f_ptr);
	imagesc(abs(exc(:,1:100)));

-------------------------------------------------------------------------

An addendum to the "peri_ear_v0.1" software package.
http://www.tnt.uni-hannover.de/~baumgarte/peri_ear/
http://geocities.com/fbaumgarte/centralear1.txt

Author: Frank Baumgarte
Please send comments and suggestions to fbaumgarte@yahoo.com

1 Introduction --------------
This document describes an example
implementation of the central ear model. This model is intended to be
used as enhancement of the peripheral ear model of this software
package. It can be used to generate masked thresholds or to compare
the perceptual similarity of audio signals. A similar model is part of
my PhD thesis (see http://www.tnt.uni-hannover.de/~baumgart). However,
the thesis is only available in German. This description includes a
simple method to generate a perceptual distance measure between two
audio signals based on the simulated masked threshold. This measure is
NOT a single number but a time-frequency function. Currently there is
no algorithm available to integrate this function in order to generate
a single number for the distance. Thus, this distance measure is
suited for a more detailed analysis of distortions, e.g. to find the
time instance and frequency were distortions occur. The resolution of
your own ears is often not sufficient for that task.

ATTENTION: For a successful distance measurement, both audio files
must have perfect time and level alignment. If there are linear
distortions, e.g. band-pass filtering of one audio signal, the model
will predict much larger distances than actually audible.

2 Central Ear Model -------------------
The input to the central ear
model are the inner hair cell (IHC) model outputs from the peripheral
ear model. The corresponding variable is called "u_ihc". The central
ear model is applied independently to the output of the IHC of each
section. Hence, the model is described here for one section only.

2.1 Low-Pass Filtering
----------------------

Several 1st order low-pass filters are used in the model. The
implementation of such a filter has one paramater k according to

  y(n) = k y(n-1) + (1-k) x(n)              (1)

x: input signal
y: output signal

This parameter can be derived from a given time constant t by

  k = exp(-1/(t Fs))                        (2)

t: time constant [s]
Fs: sampling rate of x(n) and y(n) [Hz]

In the following description the time index "n" is suppressed. From
the context it should be obvious which variable is a time signal and
should be associated with a time index. Variables corresponding to y
in (1) need to be initialized - preferably by the long-term average,
if known.

2.2 Activity
------------
From the IHC output "u_ihc" the activity "act" is calculated:

  act = (1-k_rect) (u_ihc^2) + k_rect act    (3)

with the time constant t_rect=0.002s. The time constant t_rect is
applied as explained in Sect.1 Eq. (1). The index is a shortcut for
"RECTification".

The variable "DOWNSAMPLE" in the file "param.h" of the peri_ear
package can be set to a value between 1...10, e.g. no downsampling of
u_ihc is applied if DOWNSAMPLE=1. Downsampling reduces the time
resolution, but it speeds up the processing.

2.3 Time Spreading
------------------

The activity is smeared in the time domain by convolution (*) with a
spreading function "spread_func" and an offset, "internal noise
(inoise)" is added.

  spr = act * spread_func + inoise           (4)

The spreading function is equal to a gaussian probability density function with variance

  sigma = Fs 0.005sec / sqrt(pi)             (5)

The spreading function is truncated at the -40 dB point relative to
the maximum. The truncation limits the number of filter coefficients.

The sum of all filter coefficients is normalized to 1.  Fs is the
sampling rate of "act". The output is downsampled by a factor of 256,
i.e. every 256th sample is used. In other words, from every 256 output
samples the first sample is taken and the remaining 255 samples are
discarded.

The internal noise is always positive and can be adjusted at the
different center frequency, such that the hearing threshold in quiet
is properly simulated. However, here we assume a simple "flat" noise
floor of

  inoise = 1

This signal is referred to as specific loudness "sl".

2.4 Comparison With Reference
-----------------------------

The specific loudness is computed for a reference audio file and a
"test" audio file. These quantities are refered to as "sl_ref" and
"sl_tst". Usually the goal is to detect audible changes between test
and reference with the model. Additional information about how far the
change is above the masked threshold is often desired. The coparison
is done by the following steps:

  ratio = sl_tst/sl_ref                          (6)

  r_inc = Max(ratio, 1)                          (7)
  r_incf = k_fast  r_incf + (1-k_fast) r_inc     (8)
  r_incs = k_slow  r_incs + (1-k_slow) r_inc     (9)
  r_inclp = w_slow r_incs + (1-w_slow) r_incf   (10)

  r_dec = Max(1/ratio, 1)                       (11)
  r_decf = k_fast  r_decf + (1-k_fast) r_dec    (12)
  r_decs = k_slow  r_decs + (1-k_slow) r_dec    (13)
  r_declp = w_slow r_decs + (1-w_slow) r_decf   (14)

where the short notation means:
r_inc: ratio increment
r_incf: LP-filtered ratio increment with "fast" time const.
r_incs: LP-filtered ratio increment with "slow" time const.
r_inclp: LP-filtered ratio increment with combined time const.
r_dec ... "dec"=decrement otherwise same as "inc"

The parameters are:

  t_fast = 0.002 s
  t_slow = 0.050 s
  w_slow = 0.9

Initialization:
  r_incf = 1
  r_incs = 1
  r_decf = 1
  r_decs = 1

The distance measure of the loudness change with respect to the masked
threshold is calculated as:

  thr_dist = Max((r_inclp-1)/(r_thr-1), ((r_declp-1)/(r_thr-1))^2)   (15)

It is assumed that the change is at masked threshold if thr_dist=1. If
thr_dist>1 the change is assumed to be above masked threshold.

The threshold value "r_thr" is derived as follows:

  fl_noise = fln_lf + (fln_hf -fln_lf) nu / nsect    (17)
  fl_lim = Max(fl_tone, Min(fl_noise, fl_adapt))     (18)
  rtn = rtn_lf + (rtn_hf-rtn_lf) nu / nsect          (16)
  r_thr = rtt + (1-cos(pi(fl_lim-fl_tone)/(fl_noise-fl_tone))) 0.5 (rtn - rtt)    (19)

The short notation means:
  fl_noise: fluctuation for noise-like audio signal
  fl_lim:   limited (truncated) fluctuation
  rtn:      threshold for ratio for noise-like audio signal
  rtt:      threshold for ratio for tone-like audio signal

with the parameters:

  rtn_lf = 1.4   ratio threshold for noise-like input at low frequency limit of audio band
  rtn_hf = 1.4   ratio threshold at high frequency limit of audio band
  rtt    = 1.045 ratio threshold for tone-like input
  fln_lf = 0.5   fluctuation for noise-like input at low frequency limit of audio band
  fln_hf = 0.5   fluctuation for noise-like input at low frequency limit of audio band
  fl_tone = 0.1  fluctuation for tone-like input

  nsect: number of model sections, i.e. 241
  nu: index of current section

The fluctuation "fl_adapt" is estimated by the following algorithm:

  flmax(n) = Max(sl_ref(n), k_fl flmax(n-1) + (1-k_fl) sl_ref(n))               (20)
  flmin(n) = Min(sl_ref(n), k_fl flmin(n-1) + (1-k_fl) sl_ref(n))               (21)

  fl_dev(n) = 1.00001 - fl_min(n) / fl_max(n)                             (22)
  fl_adapt(n) = Max(fl_dev(n), k_adapt fl_adapt(n-1) + (1-k_adapt) fl_dev(n))  (23)


with the parameters:

  t_fl    = 0.050 s
  t_adapt = 0.020 s

Initialization:
  flmax = 1
  flmin = 1
  fl_adapt = fl_tone

fl_dev is the non-LP filtered fluctuation derived from the deviation of the smoothed "envelope" of the maxima (flmax) and the smoothed envelope of the minima (flmin).

Further comments
----------------

The input audio signal for the peripheral ear model should have a
sampling rate larger than about 4 times the highest signal
frequency. This is necessary to reduce aliasing arising from the
nonlinear processing. For instance if you have an audio file sampled
at 44.1 kHz with a frequency range up to 20 kHz you can upsample by a
factor of 2 to get 88.2 kHz sampling rate for the model
processing. The sampling rate used in the central ear model is always
a fraction of that. It is determined by the down-sampling factor given
above, e.g. 256.

For some applications, the regular 16 bit resolution is not sufficient
for the model, since the threshold in quiet is in the range of the
LSB. If you encounter problems related to the modeling of threshold in
quiet it is recommended to use a floating-point format, e.g. ".au".

One way to visually evaluate audible distortions is to plot "thr_dist"
in the time-frequency domain. For the color mapping a compressing
function can be applied to enhance the color resolution around the
value 1. Depending on the application, it can also make sense to
generate two separate plots - one for increased loudness and one for
decreased loudness in the test file with respect to the reference. The
center frequencies fRES of the 241 model sections are approximately
given by eqs. (A.38) to (A.40) in the Ph.D. thesis. In the formulas
DELTAz=0.1, NU is the section index. (This is approximate, since the
filters have slightly shifted maxima in frequency).
