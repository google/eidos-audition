# Auditory test data.

The IHC excitation files in this directory are the golden reference files from
the original implementation by Frank Baumgarte:

1.   `1kHz_60dB_SPL_tone`: 1 kHz sinusoid with 60 dB SPL.
2.   `impulse_1000at0s_and_-1000at0.05s`: An impulse of amplitude 1000 at t=0
     and an impulse of -1000 at t=0.05 sec.

Note: The `.npy` file was produced by eidos cochlear feature generator. It is
*not* part of golden files and is kept in this directory for convenience only.
